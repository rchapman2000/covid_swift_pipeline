import subprocess 
import argparse
import sys

def find_values(line):
	dp = line.split("DP=")[1].split(";")[0]
	dp4 = line.split("DP4=")[1].split(";")[0]
	ad = line.split("AD=")[1].split(";")[0]
	return dp,dp4,ad


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-file', help="Provide pre2.vcf.")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	file = args.file

	with open(file) as f:
		content = f.readlines()
	content = [x.strip() for x in content]

	sample_name = file.split("_pre2")[0]

	new_vcf = open(sample_name+".vcf","w+") 

	list_of_positions=[]
	list_of_variants = []

	for index, line in enumerate(content):
		# Let's correct our variants
		if "#" not in line:
			position=int(line.split("\t")[1])
			dp,dp4,ad = find_values(line)
			variant_depth = ad.split(",")[1]
			
			if(position not in list_of_positions):
				list_of_positions.append(position)
				list_of_variants.append(line)

			else:
				index = list_of_positions.index(position)
				comparison_variant = list_of_variants[index]
				comparison_ad = comparison_variant.split("AD=")[1].split(";")[0].split(",")[1]

				if(int(variant_depth) >= int(comparison_ad)):
					list_of_variants[index] = line
		else:
			new_vcf.write(line+"\n")

	for variant in list_of_variants:
		new_vcf.write(variant+"\n")