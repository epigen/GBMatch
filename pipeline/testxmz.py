'''
Script originally written by Andreas, modified by Nathan.
it takes bismark alignment output file and calculates overall
methylation levels for spike-in controls.
'''

import sys
import os

def calc_methylation(bamfile, region):

	if os.path.isfile(bamfile):

		# extract SAM from BAM
		samfile = os.path.splitext(bamfile)[0] + "."+region+".sam"

		os.system("samtools view "+bamfile+" "+region+" > "+samfile)

		# parse SAM file
		file = open(samfile)

		# statistics helper variables
		count_unmethylated_lowercase = 0
		count_methylated_uppercase = 0

		count_meth_cpg = 0
		count_meth_chg = 0
		count_meth_chh = 0

		count_unmeth_cpg = 0
		count_unmeth_chg = 0
		count_unmeth_chh = 0

		while True:
			line = file.readline()
			if not line:
				break

			meth_calls = line.rstrip().split("\t")[13].split(":")[2]

			for char in meth_calls:
				if char != ".":
					if char in ["z","x","h"]:
						count_unmethylated_lowercase += 1

						if char == "z":
							count_unmeth_cpg += 1

						elif char == "x":
							count_unmeth_chg += 1

						elif char == "h":
							count_unmeth_chh += 1

					elif char in ["Z","X","H"]:
						count_methylated_uppercase += 1

						if char == "Z":
							count_meth_cpg += 1

						elif char == "X":
							count_meth_chg += 1

						elif char == "H":
							count_meth_chh += 1


		c_total = count_unmethylated_lowercase+count_methylated_uppercase
		message = ""

		if (c_total == 0):
			message += region + "_" + "count\t%d" % c_total
		else:
			cpg_total = count_meth_cpg+count_unmeth_cpg
			chg_total = count_meth_chg+count_unmeth_chg
			chh_total = count_meth_chh+count_unmeth_chh

			cpg_total_pct = float(cpg_total)/float(c_total)
			chg_total_pct = float(chg_total)/float(c_total)
			chh_total_pct = float(chh_total)/float(c_total)

			message += region + "_" + "count\t%d\n" % c_total
			message += region + "_" + "meth\t%f" % (float(count_methylated_uppercase)/float(c_total))
			#message += region + "_" + "CpG\t%f\n" % (float(count_meth_cpg)/float(cpg_total))
			#message += region + "_" + "CHG\t%f\n" % (float(count_meth_chg)/float(chg_total))
			#message += region + "_" + "CHH\t%f\n" % (float(count_meth_chh)/float(chh_total))
			#message += region + "_" + "CpG_count\t%d\n" % cpg_total
			#message += region + "_" + "CpG_pct\t%f\n" % cpg_total_pct
			#message += region + "_" + "CHG_count\t%d\n" % chg_total
			#message += region + "_" + "CHG_pct\t%d\n" % chg_total_pct
			#message += region + "_" + "CHH_count\t%d\n" % chh_total
			#message += region + "_" + "CHH_pct\t%f\n" % chh_total_pct

		print message

		#outfilename = os.path.splitext(bamfile)[0] + "." + region + ".txt"
		#out = open (outfilename, "w")
		#out.write(message+"\n")
		#out.close()

		# delete SAM file
		os.system("rm "+samfile)

	else:
		message = region + "_" + "count\tNo_file"
		print(message)
		exit(1)


if __name__ == '__main__':

	bamfile = sys.argv[1]
	region = sys.argv[2]

	calc_methylation(bamfile, region)
