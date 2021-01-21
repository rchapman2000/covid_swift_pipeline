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

fixed_file = open("filtered_variants.txt", "w+")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-name', help="Provide sample name.")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	sample_name = args.name

	for line in open("variants.txt"):
            ad=line.split("DP4=")[1].split(";")[0]
            allele_ref = int(ad.split(",")[0]) + int(ad.split(",")[1])
            allele_alt = int(ad.split(",")[2]) + int(ad.split(",")[3])
            if (allele_ref + allele_alt) > 0:
                af = allele_alt / (allele_ref + allele_alt) * 100

                if(af >= 1):
                    line_parts = line.split("\t")
                    fixed_aa_change = line_parts[2].split(":p.")[1].split(",")[0]
                    fixed_protein = line_parts[2].split(":")[1] 
                    fixed_depth = int(allele_ref + allele_alt)
                    fixed_nuc_change = line_parts[2].split(":c.")[1].split(":")[0]
                    # if(line_parts[12].rstrip()=="-"):
                    #     mat_peptide = ""
                    #     mat_peptide_nuc_change = ""
                    #     mat_peptide_aa_change = ""
                    # else:
                    #     mat_peptide = line_parts[12].split(":")[0]
                    #     mat_peptide_nuc_change = line_parts[12].split(":")[1].rstrip().split(";")[0].strip()
                    #     mat_peptide_aa_change = line_parts[12].split(":")[1].rstrip().split(";")[1].strip()
                    #visualization.write("Sample,Position,Protein,AAChange,NucleotideChange,AlleleFreq,Depth,Type,MatPeptide,MatPeptideAAChange,MatPeptideNucChange\n")
                    fixed_file.write(sample_name + "," + line_parts[4] + "," + str(fixed_protein) + "," + str(fixed_aa_change) + "," + str(fixed_nuc_change) + "," + str(af) + "," + str(fixed_depth) + "," + line_parts[1] + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 