# gff_feature_filter.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:f:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = None
feature_want = "exon"


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** gff_feature_filter.py | Written by DJP, 19/11/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Adds an exon record matching the coords of a partic feature")
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-f'):
		feature_want = arg
	else:
		print("i dont know")
		sys.exit(2)


#out_file = open(out_prefix + "_" + feature_want + ".gff", "w")

in_gff = open(in_file_name)
out_file = open(in_file_name.replace(".gff", "") + "_wExon.gff", "w")

for line in in_gff:	
	line_o = line.strip()
	if line.startswith("#"):
		out_file.write(line_o + "\n")
	else:
		line = line.strip().split("\t")
		feature = line[2]
		ID = line[8].split("ID=")[1].split(";")[0]
		if feature == feature_want:
			out_file.write(line_o + "\n")
			old_desc = line[8].split(";")
			new_desc = ""
			for i in old_desc:
				if i.startswith("ID="):
					new_desc = new_desc + ";" + i + ":exon"
				elif i.startswith("Parent="):
					new_desc = new_desc + ";Parent=" + ID 
				else:
					new_desc = new_desc + ";" + i
			new_desc = new_desc.strip(";")
			out_file.write(line[0] + "\t" + line[1] + "\t" + "exon" + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_desc + "\n")
		else:
			out_file.write(line_o + "\n")
			
			


print("\n\nFinished, Elise.\n\n")










