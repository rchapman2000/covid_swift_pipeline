import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('unwrap_fasta', help='')
	args = parser.parse_args()

	# Grabs name of file.
	fasta_name = str(args.unwrap_fasta) + ".fasta"
	new_fasta = open(fasta_name, "w")

	# Grabs only the fasta we want (not the reference fasta),
	# with the bases 201-29740 (per sarscov2_masterfile.txt).
	line_num = 0
	for line in open("repositioned_unwrap.fasta"):
		# Writes fasta header with real name
		if(line_num == 0):
			fasta_header = ">" + args.unwrap_fasta + "\n"
			new_fasta.write(fasta_header)
		elif(line_num==1):
			line_trimmed = line[201:29741]
			new_fasta.write(line_trimmed)
		line_num = line_num + 1