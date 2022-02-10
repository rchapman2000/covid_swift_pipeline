# This file takes the vcf post-processed by bcftools csq and converts it to a standardized variant file format.

import subprocess 
import argparse
import shutil
import sys
import time
from datetime import datetime
import re
import os.path
import pandas as pd
import sys
import math

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-name', help="Provide sample name.")
	parser.add_argument('-file', help="Provide processed vcf from bcftools csq.")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	sample_name = args.name
	input_vcf = args.file

	output_variants_file = open(sample_name + "_bcftools_variants.csv","w+")
	output_variants_file.write("SAMPLE,gene,AAPOS,AAREF,AASUB,TCOV,VCOV,AAFREQ,NTPOS,SNPID,NSP,NSPPOS,NSPREF,NSPSUB\n")

	for line in open(input_vcf):
		# Separate out the fields into more manageable chunks
		csq_fields = line.split("BCSQ=")
		vcf_fields = csq_fields[0]

		# Grabbing allele frequencies based on DP4
		dp4=line.split("DP4=")[1].split(";")[0]
		ad = line.split("AD=")[1].split(";")[0]
		allele_ref = int(ad.split(",")[0])
		if (int(dp4.split(",")[2]) + int(dp4.split(",")[3])) == 0:
			allele_alt = 0
		else:
			allele_alt = int(ad.split(",")[1])
		
		fixed_depth = int(dp4.split(",")[0]) + int(dp4.split(",")[1]) + int(dp4.split(",")[2]) + int(dp4.split(",")[3])

		# These fields are supposedly formatted as Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change,
		# but in practice not so
		csq_fields = csq_fields[1].split("|")
		snpid = csq_fields[-1].split("\t")[0]

		# Checking this for sanity, also...
		# some lines don't have amino acid information field (that will have less # of elements in list, and/or no snpid information)
		if (allele_ref + allele_alt) > 0 and allele_alt > 0 and len(csq_fields)>=7 and ">" in snpid:
			# For indels, base allele frequency off IMF value
			if("IMF=" in line):
				af = float(line.split("IMF=")[1].split(";")[0])
			else:
				af = allele_alt / fixed_depth
			
			change_type = csq_fields[0]
			gene = csq_fields[1]
			
			# aa_change format is in 38T>38L, time to extract info we need
			aa_change = csq_fields[-2]
			aapos = ''.join([i for i in aa_change.split(">")[0] if i.isdigit()])
			if change_type == "stop_lost" or change_type == "*stop_lost":
				aaref = "*"
			else:
				aaref = ''.join([i for i in aa_change.split(">")[0] if not i.isdigit()])

			# Getting rid of the long sgRNAs...
			if len(aaref) <= 20 and af >=0.01:
				# Figure out aa substitution

				# Synonymous mutations will either be marked, or in the case of complex, will just have aachange of 60S, vs 60S>60T 
				if change_type == "synonymous" or change_type == "*synonymous" or ">" not in aa_change:
					aasub = aaref
				elif change_type == "stop_gain" or change_type == "*stop_gain":
					aasub = "*"
				else:
					aasub = ''.join([i for i in aa_change.split(">")[1] if not i.isdigit()])

				

				# Correct snpid, which is in format 29651G>C+29652T>C+29653A>T
				ntpos = int(''.join([i for i in snpid.split(">")[0] if i.isdigit()]))
				snp_list = snpid.split("+")
				nt_ref = ""
				nt_sub = ""
				for ind_snp in snp_list:
					nt_ref = nt_ref + ''.join([i for i in ind_snp.split(">")[0] if not i.isdigit()])
					nt_sub = nt_sub + ''.join([i for i in ind_snp.split(">")[1] if not i.isdigit()])
				
				snpid = nt_ref + str(ntpos) + nt_sub
				# Time to do mature peptides!
				nsp = "-"
				nsppos = "-"
				for mature_peptide in open("mat_peptides_additions.txt"):
					mat_name = mature_peptide.split(',')[0]
					mat_start = int(mature_peptide.split(',')[1])
					mat_end = int(mature_peptide.split(',')[2])
					mat_correction = int(mature_peptide.split(',')[3])

					# Check to see if mutation falls under this mature peptide
					if ntpos >= mat_start and ntpos <= mat_end:
						# Subtracts the difference between mature peptide and start of protein
						mat_nuc_num = ntpos - mat_start + 1
						mat_aa_num = math.ceil(mat_nuc_num/3)
						
						nsppos = mat_aa_num
						nsp = mat_name

				# Output file in format SAMPLE,gene,AAPOS,AAREF,AASUB,TCOV,VCOV,AAFREQ,NTPOS,SNPID,NSP,NSPPOS,NSPREF,NSPSUB
				output_variants_file.write(",".join([sample_name,gene,aapos,aaref,aasub,str(fixed_depth),str(allele_alt),str(af),str(ntpos),snpid,nsp,str(nsppos),aaref,aasub])+"\n")
