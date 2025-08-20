import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 's:a:o:hL')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_sex_file_name  = None
in_asex_file_name  = None
outprefix    = "testout"
add_linkage_info = False

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** default_py.py | Written by DJP, 02/04/20 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Template script")	
		print("\n**** Usage****\n")
		print("python3 default_py.py -i [in file] -o [outprefix] \n\n")
		sys.exit(2)
	elif opt in ('-s'):
		in_sex_file_name  = arg
	elif opt in ('-a'):
		in_asex_file_name = arg
	elif opt in ('-o'):
		outprefix = arg
	elif opt in ('-L'):
		add_linkage_info = True
	else:
		print("i dont know")
		sys.exit(2)


seen_genes = set()
seen_tiss  = set()
seen_sp    = set()

gene_info_dict = {}
gene_tiss_sp_MF_dict = {}


in_sex_file = open(in_sex_file_name)
line_N = 0
for line in in_sex_file:
	line = line.strip().split(",")
	line_N = line_N + 1
	if line_N > 1:
		gene_id   = line[0]
		scaf      = line[1]
		XA        = line[2]
		mid_point = line[3]
		meanFPKM  = line[4]
		sp        = line[5]
		tiss      = line[6]
		sex       = line[7]
		rep_m     = line[8]
		sex_XA    = line[9]
										
		seen_genes.add(gene_id)
		gene_info_dict[gene_id] = [scaf, XA, mid_point]
		gene_tiss_sp_MF_dict[(gene_id, tiss, sp, sex)] = meanFPKM 
		
		seen_tiss.add(tiss)
		seen_sp.add(sp)

in_sex_file.close()


	

in_asex_file = open(in_asex_file_name)
line_N = 0
for line in in_asex_file:
	line = line.strip().split(",")
	line_N = line_N + 1
	if line_N > 1:
		gene_id   = line[0]
		scaf      = line[1]
		XA        = line[2]
		mid_point = line[3]
		meanFPKM  = line[4]
		sp        = line[5]
		tiss      = line[6]
		sex       = line[7]
		rep_m     = line[8]
		sex_XA    = line[9]
										
		seen_genes.add(gene_id)
		gene_info_dict[gene_id] = [scaf, XA, mid_point]
		gene_tiss_sp_MF_dict[(gene_id, tiss, sp, sex)] = meanFPKM 
		
		seen_tiss.add(tiss)
		seen_sp.add(sp)
		
in_asex_file.close()
	

print(len(seen_genes))
print(len(seen_tiss))
print(len(seen_sp))

# print(gene_tiss_sp_MF_dict)

######################################################################
### order genes by position within chr

gene_l = []

for gene in seen_genes:
	rec = gene_info_dict.get(gene)
	#print(rec)
	gene_l.append((gene, int(rec[0].split("_scf")[1]), int(rec[2])))


gene_ls_1 =  sorted(gene_l,    key=lambda x: x[2])
gene_ls_2 =  sorted(gene_ls_1, key=lambda x: x[1])

####################################################################$#$
## get sp/tiss order

species_order = ["Tbi", "Tce", "Tcm", "Tpa", "Tps", "Tte", "Tms", "Tsi", "Tge", "Tdi"]

tissue_sex_order_want  = set(['A_F', 'A_M', 'B_F', 'B_M', 'DG_F', 'DG_M', 'FB_F', 'FB_M', 'Fe_F', 'Fe_M', 'Ta_F', 'Ta_M',  'Gu_F', 'Gu_M', 'Go_F', 'AG_M', 'Te_M'])
tissue_order  = ['A', 'B', 'DG', 'FB', 'Fe', 'Ta', 'Gu', 'Go', 'AG', 'Te']
sex_order  = ["F", "M"]


my_sp_order = []
for el in species_order:
	if el in seen_sp:
		my_sp_order.append(el)
print(my_sp_order)


#############################################################################
### output

out_file = open(outprefix + "_fheatmap.csv", "w")

line_N = 0

for gene in gene_ls_2:
	line_N = line_N + 1
	gene_id = gene[0]
	gene_info = gene_info_dict.get(gene_id)

	FPKM_list = ""
	FPKM_samps = ""
	for t in tissue_order:
		for s in my_sp_order:
			for x in sex_order:
				if t + "_" + x in tissue_sex_order_want:
					rec = gene_tiss_sp_MF_dict.get((gene_id, t, s, x))
					#print(rec)
					if rec == None:
						rec = "NA"

					FPKM_list =  FPKM_list + "," + rec
					FPKM_samps = FPKM_samps + "," + s + "_" + t + "_"  + x
					
	if line_N == 1:
		out_file.write("gene_id,scaf,XA,mid_point" + FPKM_samps + "\n")
	
	out_file.write(gene_id + "," + gene_info[0]  + "," +  gene_info[1]  + "," +  gene_info[2]  + FPKM_list + "\n")
	
	


#########################################################################################################
### add linkage info for Tce

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
	

if add_linkage_info == True:
	out_file = open(outprefix + "_fheatmap.csv")
	out_file_2 = open(outprefix + "_fheatmap.csv.temp", "w")
	
	LINE_N = 0
	for line in out_file:
		line = line.strip().split(",")
		LINE_N = LINE_N + 1
		outline = ""
		if LINE_N == 1:
			for i in line:
				outline = outline + "," + i
			
			outline = outline.strip(",") + "," + "LG"
		else:
			scaf = line[1]
			LG = LG_dict_scafs.get(scaf)
			if LG == None:
				LG = "NA"
			else:
				LG = oldLG_to_newLG_dict.get(LG)
			
			for i in line:
				outline = outline + "," + i
			
			outline = outline.strip(",") + "," + LG
				
			#print(line)
			#print(LG)		
		
		out_file_2.write(outline + "\n")
	out_file_2.close()
	out_file.close()
	os.rename(outprefix + "_fheatmap.csv.temp", outprefix + "_fheatmap.csv") 


