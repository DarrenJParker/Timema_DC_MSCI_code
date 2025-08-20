import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_file_name = None
outprefix    = "testout"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** default_py.py | Written by DJP, 24/08/23 in Python 3.5 in Porto ****\n")
		print("add LG info to Tce count files")	
		print("\n**** Usage****\n")
		print(" python3 accessory_scripts/add_HiC_LG_to_Tce_counts.py -i data/readcounts/to_Tce_H2E.counts_wGLSL.csv  \n\n")
		sys.exit(2)
	elif opt in ('-i'):
		in_file_name  = arg
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)

if in_file_name == None:
	print("\n\nERROR. No input file\n\n")
	

LG_dict = {
"LG1":[1,8],
"LG2":[9,26],
"LG3":[27,41],
"LG4":[42,49],
"LG5":[50,60],
"LG6":[61,68], # and 111
"LG7":[69,79],
"LG8":[80,101],
"LG9":[102,110],
"LG10":[112,124],
"LG11":[125,134],
"LG12":[135,145],
"LG13":[146,157]
}

LG_dict_scafs = {}

scaf_in_LG = set()
for el in LG_dict:
	rec = LG_dict.get(el)
	
	for i in range(rec[0], rec[1] + 1):
		LG_dict_scafs["Tce_LRv5a_scf" + str(i)] = el
		scaf_in_LG.add("Tce_LRv5a_scf" + str(i))
		


LG_dict_scafs["Tce_LRv5a_scf" + str(111)] = "LG6"
scaf_in_LG.add("Tce_LRv5a_scf" + str(111))
	
	
### get scaf col N
scaf_index = None

in_file = open(in_file_name)
line_N = 0
for line in in_file:
	line_N = line_N + 1
	if line_N == 1:
		line = line.rstrip("\n").split(",")
		scaf_index = line.index("scaf")

in_file.close()


oldLG_to_newLG_dict = {
"LG3" : "HiCLG1",
"LG2" : "HiCLG2",
"LG8" : "HiCLG3",
"LG13" : "HiCLG4",
"LG7" : "HiCLG5",
"LG1" : "HiCLG6",
"LG5" : "HiCLG7",
"LG9" : "HiCLG8",
"LG10" : "HiCLG9",
"LG6" : "HiCLG10",
"LG4" : "HiCLG11",
"LG11" : "HiCLG12",
"LG12" : "HiCLG13"
}
	

### read in file
in_file = open(in_file_name)
out_temp = open(in_file_name.replace(".counts_wGLSL.csv", ".counts_wGLSLLG.csv"), "w")
line_N = 0
for line in in_file:
	line_N = line_N + 1
	if line_N == 1:
		out_temp.write(line.strip() + ",LG\n") 
	
	else:
		line = line.rstrip("\n").split(",")
		LG = "NA"
		if line[scaf_index] in scaf_in_LG:
			LG = LG_dict_scafs.get(line[scaf_index])
		
		out_line = ""
		for i in line:
			out_line = out_line + "," + i
		
		new_LG = oldLG_to_newLG_dict.get(LG)
		if new_LG == None:
			new_LG = "NA"
		
		out_line = out_line.strip(",") + "," + new_LG
		out_temp.write(out_line + "\n")

out_temp.close()



# out_temp = open(in_file_name.replace(".counts_wGLSL.csv", ".counts_wGLSLLG.csv.temp"))
# out_file = open(in_file_name.replace(".counts_wGLSL.csv", ".counts_wGLSLLG.csv"), "w")
# 
# LG_index = None
# line_N = 0
# 
# for line in out_temp:
# 	line = line.strip().split(",")
# 	line_N = line_N + 1
# 	if line_N == 1:
# 		LG_index = line
# 
# 
# 

