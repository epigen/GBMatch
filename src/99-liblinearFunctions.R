library(LiblineaR)
library(ROCR)
library(caret)


#####Functions#############################################################################

#multi_roc function
multi_roc=function(class,decisionValues,true_labels){
  true_labels[true_labels==class]="_sel_"
  true_labels[true_labels!="_sel_"]="_noSel_"
  true_labels[true_labels=="_sel_"]=1
  true_labels[true_labels=="_noSel_"]=0
  
  nas=is.na(as.data.frame(decisionValues)[,class])
  #print(class)
  if (length(unique(true_labels[!nas]))==2){
    pred <- prediction(as.data.frame(decisionValues)[,class],true_labels)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    auc <- performance(pred, measure = "auc")
    dt=data.table(fpr=perf@x.values[[1]],tpr=perf@y.values[[1]],class=class,auc=unlist(auc@y.values))}else{
      dt=data.table(fpr=0.5,tpr=0.5,class=class,auc=0)
    }  
  return(dt)
}

#binarize vectors
binarize=function(vector,low,high){
  ret=ifelse(vector<=quantile(vector,low,na.rm=TRUE),"low",ifelse(vector>=quantile(vector,high,na.rm=TRUE),"high",NA))
  return(as.character(ret))
}

#liblinear function

find_param=function(data,labels,type=NULL,cost=NULL,scaleAndCenter=FALSE){
  #parameter selection
  
  if(is.null(type)){
    tryTypes=c(0:7)}else{tryTypes=type}
  if(is.null(cost)){
    tryCosts=c(1000,1,0.001)}else{tryCosts=cost}
  
  bestCost=NA
  bestScore=-Inf
  bestType=NA
  
  set.seed(1234)
  train_test <- createDataPartition(labels, p = .5,list = FALSE, times = 1)
  if(scaleAndCenter==TRUE){
  train_train_data=scale(data[-train_test,],center=TRUE,scale=TRUE)
  }else{train_train_data=data[-train_test,]}
  
  train_train_annot=labels[-train_test]
  if(scaleAndCenter==TRUE){
  train_test_data=scale(data[train_test,],attr(train_train_data,"scaled:center"),attr(train_train_data,"scaled:scale"))
  }else{train_test_data=data[train_test,]}
  train_test_annot=labels[train_test]
  
  na_cols_train=apply(train_train_data,2,function(x)sum(is.na(x)))
  na_cols_test=apply(train_test_data,2,function(x)sum(is.na(x)))
  remove=c(names(na_cols_train[na_cols_train>0]),names(na_cols_test[na_cols_train>0]))
  
  train_train_data=train_train_data[,!colnames(train_train_data)%in%remove]
  train_test_data=train_test_data[,!colnames(train_test_data)%in%remove]
  
  for(ty in tryTypes){
    for(co in tryCosts){
      
      m=LiblineaR(data=train_train_data,target=train_train_annot,type=ty,cost=co,bias=TRUE,verbose=FALSE,svr_eps=0.01)
      p=predict(m,train_test_data,proba=FALSE,decisionValues=TRUE)
      roc_mat=unique(rbindlist(sapply(colnames(p$decisionValues),multi_roc,decisionValues=p$decisionValues,true_labels=as.character(train_test_annot),simplify=FALSE))[,c("class","auc"),with=FALSE])
      meanAUC=mean(roc_mat$auc,na.rm=TRUE)
      score=meanAUC
      
      cat("Type: ", ty,"\n")
      cat("Results for C=",co," : ",score," meanAUC.\n",sep="")
      if(score>bestScore){
        bestCost=co
        bestScore=score
        bestType=ty
      }
    }
  }
  cat("Best model: ",bestType,"--Best cost: ",bestCost,"--Best meanAUC: ",bestScore)
  
  return(list(type=bestType,cost=bestCost))
  
}




run_pred=function(data,labels,samples,cross=10,param,scaleAndCenter=FALSE){
  feature_mat_all=data.table()
  decision_mat_all=data.table()
  decision_list=list()
  i=0
  set.seed(1234)
  folds=createFolds(labels,k=cross)
  for (fold_name in names(folds)){
    i=i+1
    message(fold_name)
    fold=folds[[fold_name]]
    if(scaleAndCenter==TRUE){
      print("Attention: scaling and centering data.")
      train_data=scale(data[-fold,,drop=FALSE],center=TRUE,scale=TRUE)
    }else{train_data=data[-fold,,drop=FALSE]}
    train_annot=labels[-fold]
    
    if(scaleAndCenter==TRUE){
      test_data=scale(data[fold,,drop=FALSE],attr(train_data,"scaled:center"),attr(train_data,"scaled:scale"))
    }else{test_data=data[fold,,drop=FALSE]}
    
    
    test_annot=labels[fold]
    test_samples=samples[fold]
    
    na_cols_train=apply(train_data,2,function(x)sum(is.na(x)))
    na_cols_test=apply(test_data,2,function(x)sum(is.na(x)))
    remove=c(names(na_cols_train[na_cols_train>0]),names(na_cols_test[na_cols_train>0]))
    
    train_data=train_data[,!colnames(train_data)%in%remove,drop=FALSE]
    test_data=test_data[,!colnames(test_data)%in%remove,drop=FALSE]
      
    # Re-train best model with best cost value.
    m=LiblineaR(data=train_data,target=train_annot,type=param$type,cost=param$cost,bias=TRUE,verbose=FALSE,svr_eps=0.01)
    
    feature_mat=as.data.table(t(m$W))
    feature_names=colnames(m$W)
    feature_mat[,region:=feature_names,]
    feature_mat[,fold:=fold_name,]
    feature_mat_long=melt(feature_mat,id.vars=c("fold","region"))
    feature_mat_all=rbindlist(list(feature_mat_all,feature_mat_long))
    
    # Make prediction
    pr=FALSE
    dv=FALSE
    if(param$type==0 || param$type==7) pr=TRUE
    if(param$type%in%c(0:7)) dv=TRUE
    p=predict(m,test_data,proba=pr,decisionValues=dv)
    
    decision_mat=as.data.table(p$decisionValues)
    decision_mat[,true_label:=as.character(test_annot),]
    decision_mat[,pred_label:=as.character(p$predictions),]
    decision_mat[,fold:=fold_name,]
    decision_mat[,sample:=test_samples,]
    decision_mat_long=melt(decision_mat,id.vars=c("true_label","pred_label","fold","sample"))
    decision_mat_all=rbindlist(list(decision_mat_all,decision_mat_long))
    
    #use for sanity check (only works for 2 categories)
    #decision_list$predictions[[i]]=p$decisionValues[,1]
    #decision_list$labels[[i]]=as.character(test_annot)
    #decision_list$predlabels[[i]]=as.character(p$predictions)
    #pred=prediction(unlist(decision_list$predictions),unlist(decision_list$labels))
    #perf=performance(pred, 'tpr', 'fpr')
    #plot(perf)
  }
  
  
  #train model on ALL data for classification task of other data e.g. RRBS
  message("Training final model...")
  if(scaleAndCenter==TRUE){
  data_scaled=scale(data,scale=TRUE,center=TRUE)
  m_combined=LiblineaR(data=data_scaled,target=labels,type=param$type,cost=param$cost,bias=TRUE,verbose=FALSE,svr_eps=0.01)
  attr_center=attr(data_scaled,"scaled:center")
  attr_scale=attr(data_scaled,"scaled:scale")
  }else{
    m_combined=LiblineaR(data=data,target=labels,type=param$type,cost=param$cost,bias=TRUE,verbose=FALSE,svr_eps=0.01)
    attr_center=NA
    attr_scale=NA
  }
  return(list(features=feature_mat_all,decisions=decision_mat_all,model=m_combined,attr_center=attr_center,attr_scale=attr_scale)) 
}

#assess and plot prediction also with random labels
check_prediction=function(data,labels,samples,cross=10,nReps=10,type=NULL,cost=NULL,scaleAndCenter=FALSE){
  
  
  param=find_param(data,labels,type,cost,scaleAndCenter=FALSE)
  pred=run_pred(data=data,labels=labels,samples=samples,param=param,cross = cross,scaleAndCenter=scaleAndCenter)
  
  #random labels
  nReps=nReps
  set.seed(1234)
  rand_labs=lapply(seq_len(nReps), function(x) sample(labels))
  i=0
  decision_mat_all=data.table()
  for (rand_lab in rand_labs){
    i=i+1
    param_rand=find_param(data,rand_lab,type,cost,scaleAndCenter=FALSE)
    pred_rand=run_pred(data=data,labels=rand_lab,samples=samples,param=param_rand,cross = cross,scaleAndCenter=scaleAndCenter)
    decision_mat=pred_rand$decisions
    decision_mat[,rep:=i,]
    decision_mat_all=rbindlist(list(decision_mat_all,decision_mat))
  }
  
  decision_mat_all[,rand:=TRUE,]
  decision_mat_all=rbindlist(list(decision_mat_all,pred$decisions[,c("rand","rep"):=list(FALSE,1),]),use.names = TRUE)
  
  decision_mat_wide=reshape(decision_mat_all,idvar=c("true_label","pred_label","fold","sample","rand","rep"),timevar="variable",direction="wide")
  setnames(decision_mat_wide,names(decision_mat_wide),gsub("value.","",names(decision_mat_wide)))
  
  
  roc_mat_cv=decision_mat_wide[,rbindlist(sapply(names(decision_mat_wide)[-c(1:6)],multi_roc,decisionValues=eval(parse(text=paste0("data.frame(",paste0(sapply(names(decision_mat_wide)[-c(1:6)],function(x)paste0(x,"=",x)),collapse = ","),")"))),true_labels=as.character(true_label),simplify=FALSE)),by=c("rep","rand")]
  
  
  auc_annot=roc_mat_cv[,list(mean_auc=mean(auc)),by=c("class","rand")]
  auc_annot[,x:=0.9,]
  auc_annot[,y:=ifelse(rand==TRUE,0.03,0.1),]
  
  auc_annot[,label:=ifelse(rand==TRUE,paste0("Control=",round(as.numeric(mean_auc),2)),paste0("AUC=",round(as.numeric(mean_auc),2))),]
  
  classes=unique(roc_mat_cv$class)
  inernal_plot_list=list()
  for (plot_class in classes){
    pl=ggplot(roc_mat_cv[class==plot_class],aes(x=fpr,y=tpr,col=rand))+geom_line(aes(group=paste0(rep,rand),alpha=rand))+geom_text(data=auc_annot[class==plot_class],alpha=1,aes(x=x,y=y,label=label))+scale_color_manual(values=c("TRUE"="darkgrey","FALSE"="blue"))+scale_alpha_manual(values=c("TRUE"=0.4,"FALSE"=1))+annotate(geom="text",x=0,y=1,hjust=0,vjust=1,label=plot_class)
    inernal_plot_list[[plot_class]]=pl
  }
  
  return(list(plot=inernal_plot_list,decision_mat=decision_mat_wide,model=pred$model,attr_center=pred$attr_center,attr_scale=pred$attr_scale,auc=mean(auc_annot[rand==FALSE]$mean_auc),auc_rand=mean(auc_annot[rand==TRUE]$mean_auc)))
}
#############################################################################################

##################for predicting properties from methylation prediction######################

#########prepare the data#################
prepare_data=function(meth_data_dt,annotation_all,min_na_ratio=0.999){
  meth_data_mat=dcast(meth_data_dt,id~region,value.var="methyl")
  row.names(meth_data_mat)=meth_data_mat$id
  meth_data_mat$id=NULL
  
  region_na_ratios=apply(meth_data_mat,2,function(x){sum(!is.na(x))/length(x)})
  min_na_ratio=min_na_ratio
  
  #check for low coverage regions --> select regions to include
  meth_data_mat_sel=meth_data_mat[,region_na_ratios>min_na_ratio]
  #check for low coverage sample --> select samples to include
  sample_na_ratios=apply(meth_data_mat_sel,1,function(x){sum(!is.na(x))/length(x)})
  meth_data_mat_sel=meth_data_mat_sel[sample_na_ratios>0.8,]
  
  meth_data_imputed=t(impute.knn(t(meth_data_mat_sel),k=5)$data)
  region_na_ratios_red=region_na_ratios[region_na_ratios>min_na_ratio]
  
  #reduce annotation to selected samples in meth data
  annotation=annotation_all[N_number_seq%in%row.names(meth_data_imputed)]
  
  #make sure annotation is in the same order as the matrix
  sample_match=match(annotation$N_number_seq,row.names(meth_data_imputed))
  stopifnot(length(row.names(meth_data_imputed)[!row.names(meth_data_imputed)%in%annotation$N_number_seq])==0)
  stopifnot(length(annotation$N_number_seq[!annotation$N_number_seq%in%row.names(meth_data_imputed)])==0)
  stopifnot(nrow(annotation[N_number_seq%in%annotation[duplicated(N_number_seq)]$N_number_seq])==0)
  stopifnot(all(annotation$N_number_seq[order(sample_match)]==row.names(meth_data_imputed)))
  annotation=annotation[order(sample_match)]
  annotation[annotation==""]=NA
  stopifnot(all(annotation$N_number_seq==row.names(meth_data_imputed)))
  return(list(meth_data_imputed=meth_data_imputed,annotation=annotation))
}

#########prepare the data for predefined features#################
prepare_data_fixed=function(meth_data_dt,annotation_all,features){
  meth_data_mat=dcast(meth_data_dt[region%in%features],id~region,value.var="methyl")
  row.names(meth_data_mat)=meth_data_mat$id
  meth_data_mat$id=NULL
  
  #region_na_ratios=apply(meth_data_mat,2,function(x){sum(!is.na(x))/length(x)})
  #min_na_ratio=min_na_ratio
  
  #check for low coverage regions --> select regions to include
  #meth_data_mat_sel=meth_data_mat[,region_na_ratios>min_na_ratio]
  #check for low coverage sample --> select samples to include
  sample_na_ratios=apply(meth_data_mat,1,function(x){sum(!is.na(x))/length(x)})
  meth_data_mat=meth_data_mat[sample_na_ratios>0.8,]
  
  meth_data_imputed=t(impute.knn(t(meth_data_mat),k=5)$data)
  
  #reduce annotation to selected samples in meth data
  annotation=annotation_all[N_number_seq%in%row.names(meth_data_imputed)]
  
  #make sure annotation is in the same order as the matrix
  sample_match=match(annotation$N_number_seq,row.names(meth_data_imputed))
  stopifnot(length(row.names(meth_data_imputed)[!row.names(meth_data_imputed)%in%annotation$N_number_seq])==0)
  stopifnot(length(annotation$N_number_seq[!annotation$N_number_seq%in%row.names(meth_data_imputed)])==0)
  stopifnot(nrow(annotation[N_number_seq%in%annotation[duplicated(N_number_seq)]$N_number_seq])==0)
  stopifnot(all(annotation$N_number_seq[order(sample_match)]==row.names(meth_data_imputed)))
  annotation=annotation[order(sample_match)]
  annotation[annotation==""]=NA
  stopifnot(all(annotation$N_number_seq==row.names(meth_data_imputed)))
  return(list(meth_data_imputed=meth_data_imputed,annotation=annotation))
}






#############prediction##########################################
meth_pred_analysis=function(meth_data_imputed,annotation,column_annotation,set_targets=NULL,to_analyse,meth_sel,set_scale=FALSE,type=4,cost=1,cross=NULL,nReps=10){
  for(selected in to_analyse){
    
    if(file.exists(paste0("dat_",selected,meth_sel,".RData"))){
      load(file=paste0("dat_",selected,meth_sel,".RData"))
      
      pdf(paste0("roc_",selected,meth_sel,".pdf"),height=3.5,width=4.5)
      for (name in names(pl)){
        plot_list=pl[[name]]$plot
        for (i in c(1:length(plot_list))){
          print(pl[[name]]$plot[[i]]+ggtitle(name)) 
        }
      }
      dev.off()
      next
    }
    if (is.null(set_targets)){
      targets=column_annotation[[selected]]}else{
      targets=set_targets  
    }
    pl=list()
    for (target in targets){
      if(!target%in%names(annotation)){
        print(paste0(target,"not found!"))
        next
      }
      annotation_sel=annotation[!is.na(get(target))]
      ##Check if enough classes  
      nclasses=length(na.omit(unique(annotation_sel[,get(target)])))
      if(nclasses<2){
        print(paste0(target,": only one class. Skipping"))
        next}
      if(nclasses==2){
        tab=table(annotation_sel[,get(target)])
        remove=names(tab[tab<2])
        annotation_sel=annotation_sel[!get(target)%in%remove]
        nclasses=length(na.omit(unique(annotation_sel[,get(target)])))
        if(nclasses<2){
          print(paste0(target,": only one class. Skipping"))
          next
        }
      }
      #decide if categorical or continous and make categorical
      if (length(unique(annotation_sel[,get(target)]))>6&!"character"%in%is(annotation_sel[,get(target)])){
        data_type="cont"
      }else{data_type="cate"}
      if (data_type=="cate"&nclasses>20){
        print(paste0(target,": more than 20 classes. Skipping"))
        next
      }
      if (data_type=="cont"){
        annotation_sel[,eval(target):=binarize(get(target),0.2,0.8)]
        data_type="cate"
      }
      if(any(grepl("^[0-9]",annotation_sel[,get(target)],perl=TRUE))){
        annotation_sel[,eval(target):=paste0("x",get(target))]
      }
      ##now actually perform prediction   
      annotation_sel=annotation_sel[!is.na(get(target))]
      meth_data_sel=meth_data_imputed[row.names(meth_data_imputed)%in%annotation_sel$N_number_seq,]
      if(any(annotation_sel$N_number_seq!=row.names(meth_data_sel))){stop("features and targets don't match")}
      ############turn scaling on and off here#################################################################
      ##use scaling for m-values and no scaling for beta-values. Scaled m-values seem too work slightly better than unscaled beta values, but scaled beta values perform equally well. However logically, it doesn't really make sense to scale beta-values.
      
      if (is.null(cross)){cross=nrow(meth_data_sel)}
      pl[[target]]=check_prediction(data=meth_data_sel,labels=gsub("/|\\+","",gsub(" |-","_",annotation_sel[,get(target),])),samples=annotation_sel$N_number_seq,cross=cross,nReps=nReps,scaleAndCenter=set_scale,type=type,cost=cost)
      #########################################################################################################
    }
    save(pl,file=paste0("dat_",selected,meth_sel,".RData"))
    pdf(paste0("roc_",selected,meth_sel,".pdf"),height=3.5,width=4.5)
    for (name in names(pl)){
      plot_list=pl[[name]]$plot
      for (i in c(1:length(plot_list))){
        print(pl[[name]]$plot[[i]]+ggtitle(name)) 
      }
    }
    dev.off()
  }
}


