import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:r:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_gff_name = None
read_count_filename = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** default_py.py | Written by DJP, 09/02/23 in Python 3.5 in Bangor, UK ****\n")
		print("\n**** Usage****\n")
		print("see /Users/drp22jhz/Documents/University/Lausanne/Timema_transcriptomes/Asex_males/2_read_counts.sh")
		sys.exit(2)
	elif opt in ('-i'):
		in_gff_name  = arg
	elif opt in ('-r'):
		read_count_filename = arg
	else:
		print("i dont know")
		sys.exit(2)


### get gene pos
gene_mid_pos_dict  = {}
in_gff = open(in_gff_name)
for line in in_gff:
	if not line.startswith("#"):
		line = line.rstrip("\n")
		feature = line.split("\t")[2]
		if feature == "gene":
			geneID = line.split("\t")[8].split("ID=")[1].split(";")[0]
			mid_point =  int((int(line.split("\t")[4]) + int(line.split("\t")[3])) / 2)
			gene_mid_pos_dict[geneID] = mid_point
					
in_gff.close()		

outfile = open(read_count_filename.replace(".csv","L.csv"), "w")
read_counts = open(read_count_filename)

line_N = 0
for line in read_counts:
	line_N = line_N + 1
	line = line.strip()
	if line_N == 1:
		outfile.write(line + "," + "mid_point" + "\n")	
	else:
		gene_name = line.split(",")[0]
		mid_point = gene_mid_pos_dict.get(gene_name)
	
		outfile.write(line + "," + str(mid_point) + "\n")
		
