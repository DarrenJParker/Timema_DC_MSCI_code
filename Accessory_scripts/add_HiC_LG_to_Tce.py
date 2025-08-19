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
		print("\n**** default_py.py | Written by DJP, 15/08/23 in Python 3.5 in Porto ****\n")
		print("add LG info to Tce sliding window coverage files")	
		print("\n**** Usage****\n")
		print(" python3 accessory_scripts/add_HiC_LG_to_Tce.py -i /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/angsD_LR_sliding_window/Tms_M_25_15055_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt \n\n")
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
LG_to_scafs_dict = {}
seen_LG = set()
scaf_in_LG = set()
for el in LG_dict:
	rec = LG_dict.get(el)
	
	for i in range(rec[0], rec[1] + 1):
		LG_dict_scafs["Tce_LRv5a_scf" + str(i)] = el
		scaf_in_LG.add("Tce_LRv5a_scf" + str(i))
		if el not in seen_LG:
			seen_LG.add(el)
			LG_to_scafs_dict[el] = ["Tce_LRv5a_scf" + str(i)]
		else:
			rec1 = LG_to_scafs_dict.get(el)
			rec1.append("Tce_LRv5a_scf" + str(i))
			LG_to_scafs_dict[el] = rec1


LG_dict_scafs["Tce_LRv5a_scf" + str(111)] = "LG6"
scaf_in_LG.add("Tce_LRv5a_scf" + str(111))


scaf_lens_dict = {}

### read in file
in_file = open(in_file_name)
out_temp = open("temp_HiC_add.temp", "w")
line_N = 0
for line in in_file:
	line_N = line_N + 1
	if line_N == 1:
		out_temp.write(line.strip() + "\tLG\n") 
	
	else:
		line = line.rstrip("\n").split("\t")
		LG = "NA"
		#print(line)
		scaf_lens_dict[line[0]] = int(line[2])
		if line[0] in scaf_in_LG:
			LG = LG_dict_scafs.get(line[0])
		
		out_line = ""
		for i in line:
			out_line = out_line + "\t" + i
		
		out_line = out_line.strip() + "\t" + LG
		out_temp.write(out_line + "\n")

out_temp.close()



#################
LGs = ["LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13"]


seen_scafLG = set()
seen_LG = set()

for L in LGs:
	out_temp = open("temp_HiC_add.temp")
	out_temp_2 = open("temp_HiC_add" + L + ".temp2", "w")
	for line in out_temp:
		line_o = line.rstrip("\n")
		line = line.rstrip("\n").split("\t")
		LG = line[5]
		scaf = line[0]

		if LG == L:
			if scaf + LG not in seen_scafLG:
				seen_scafLG.add(scaf + LG)
				if LG not in seen_LG:
					seen_LG.add(LG)
					out_temp_2.write(line_o + "\n")
				else:
					out_temp_2.write(last_scaf + "\t" + str(e_coord) + "\t" + str(e_coord + 200000) +  "\tNA\tNA\t" + LG + "\n" +  line_o  + "\n")
			else:
				out_temp_2.write(line_o  + "\n")
			
			last_scaf = line[0]
			s_coord = int(line[1])
			e_coord = int(line[2])
		
	out_temp.close()
	out_temp_2.close()
	

outfile = open(in_file_name.replace(".txt", "_wLG.txt.temp"), "w")	
outfile.write("scaf_name\tstart\tend\ttot_feat\tcov\tLG\tLG_start\tLG_end\n")


path = "."
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(".temp2"):
			#print (os.path.join(path, name))
		
			out_temp_2 = open(os.path.join(path, name))
			
			n_end = 0
			for line in out_temp_2:
				line = line.rstrip("\n").split("\t")
				s_coord = int(line[1])
				e_coord = int(line[2])
				LG = line[5]
				scaf = line[0]
				
			
				if  e_coord < n_end:
					if line[3] == "NA":
						n_end = ((e_coord - s_coord)) + n_end
					else:
						n_end = ((e_coord - s_coord) / 2) + n_end
				else:
					n_end = e_coord
			
				n_start = n_end - (e_coord - s_coord)
				
				out_line = ""
				for i in line:
					out_line = out_line + "\t" + i
		
				out_line = out_line.strip() + "\t" + str(n_start) + "\t" + str(n_end)
				outfile.write(out_line + "\n")
				# print(line)
				# print(str(n_start) + "\t" + str(n_end))
				
			out_temp_2.close()	
			
			
out_temp = open("temp_HiC_add.temp")			
for line in out_temp:
	line_o = line.rstrip("\n")
	line = line.rstrip("\n").split("\t")
	LG = line[5]
	scaf = line[0]
	if LG == "NA":
		out_line = ""
		for i in line:
			out_line = out_line + "\t" + i
		
		out_line = out_line.strip() + "\tNA\tNA"
		outfile.write(out_line + "\n")
				# print(line)

outfile.close()

### rename LGs by size


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





outfile = open(in_file_name.replace(".txt", "_wLG.txt.temp"))	
outfile_r = open(in_file_name.replace(".txt", "_wLG.txt"), "w")	

line_N = 0
for line in outfile:
	
	line_N = line_N + 1
	line = line.strip().split("\t")
	if line_N == 1:
		outfile_r.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" +line[4] + "\t" + line[5] + "\t" +line[6] + "\t" +line[7] + "\n")
	else:
		old_LG = line[5]
		new_LG = oldLG_to_newLG_dict.get(old_LG)
		if new_LG == None:
			new_LG = "NA"
		# print(line)
		# print(old_LG)
		# print(new_LG)
		outfile_r.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" +line[4] + "\t" + new_LG + "\t" +line[6] + "\t" +line[7] + "\n")
		


#### tidy up


path = "."
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith(".temp2"):
			os.remove(os.path.join(path, name))
			
os.remove("temp_HiC_add.temp")		
os.remove(in_file_name.replace(".txt", "_wLG.txt.temp"))


# for el in LG_to_scafs_dict :
# 	rec = LG_to_scafs_dict.get(el)
# 	total_len = 0
# 	for s in rec:
# 		total_len = total_len + scaf_lens_dict.get(s)
# 		#print(s)		
# 	print(el)
# 	#print(rec)
# 	print(total_len)
# 
# 
# 
# 

