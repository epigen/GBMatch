#! /bin/bash


bam1=$1
bam2=$2
out_dir=$3
out=$4
min_reads=$5

samtools_path=/cm/shared/apps/samtools/1.3/bin/samtools
methclone_path=/data/groups/lab_bock/shared/resources/tools/methclone/bin/

tempdir=$out_dir/$out
mkdir -p $tempdir || exit 1

export TMPDIR=$tempdir

bam1_sorted=$(dirname $bam1)/$(basename $bam1 .bam)_sorted.bam
bam2_sorted=$(dirname $bam2)/$(basename $bam2 .bam)_sorted.bam

echo "bam1: " $bam1_sorted
echo "bam2: " $bam2_sorted
echo "out: " $out

if [ ! -f $bam1_sorted.bai ]; then
	echo "bam1: Sorting and indexing!"
	$samtools_path sort -o $bam1_sorted $bam1
	$samtools_path index $bam1_sorted
fi

if [ ! -f $bam2_sorted.bai ]; then
	echo "bam2: Sorting and indexing!"
	$samtools_path sort -o $bam2_sorted $bam2
	$samtools_path index $bam2_sorted
fi

if [ ! -f $tempdir/$out.output.txt.gz ] && [ ! -f $tempdir/$out.output.txt ]; then
	$methclone_path/methclone $bam1_sorted $bam2_sorted $tempdir/$out.output.txt.gz $out $min_reads
fi
gunzip $tempdir/$out.output.txt.gz
