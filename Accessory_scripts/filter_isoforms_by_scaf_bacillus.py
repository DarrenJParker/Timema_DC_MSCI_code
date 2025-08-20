import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'c:i:s:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

count_file_name   = None
isoform_file_name = None
scaf_col_name = None


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** filter_isoforms_by_scaf.py| Written by DJP, 27/12/23 in Python 3.5 in Porto ****\n")
		print("filters isoforms to have only those on good scafs")	
		print("\n**** Usage****\n")
		print("python3 accessory_scripts/filter_isoforms_by_scaf_bacillus.py -c output/orthologs/Brsri_droso_annot_gene.gtf -i temp/Brsri_v3_aa_longest.fa \n see 4_orthologs.sh\n")
		sys.exit(2)
	elif opt in ('-i'):
		isoform_file_name  = arg
	elif opt in ('-c'):
		count_file_name  = arg
	elif opt in ('-s'):
		scaf_col_name  = arg
	else:
		print("i dont know")
		sys.exit(2)



scaf_col_name_index = None
all_want_scafs = set(
[
"Brsri_v3_scf1",  
"Brsri_v3_scf2",  
"Brsri_v3_scf3",  
"Brsri_v3_scf4_1", 
"Brsri_v3_scf4_2",
"Brsri_v3_scf5",   
"Brsri_v3_scf6", 
"Brsri_v3_scf7",  
"Brsri_v3_scf8",  
"Brsri_v3_scf9_1", 
"Brsri_v3_scf9_2", 
"Brsri_v3_scf10", 
"Brsri_v3_scf11",  
"Brsri_v3_scf12", 
"Brsri_v3_scf13",
"Brsri_v3_scf14",
"Brsri_v3_scf15", 
"Brsri_v3_scf16",
"Brsri_v3_scf17",
"Brsri_v3_scf18"]
)

all_kept_scafs = set()
all_kept_genes = set()
all_kept_genes_list = []

line_N = 0
count_file = open(count_file_name)
for line in count_file:
	line_N = line_N + 1	
	line = line.rstrip("\n").split("\t")
	scaf_name_c = line[0]
	gene_name = line[8].split("gene_id")[1].split(";")[0].strip().strip('"')

	if scaf_name_c in all_want_scafs:
		all_kept_genes.add(gene_name)
		all_kept_scafs.add(scaf_name_c)
		all_kept_genes_list.append(gene_name)
print(all_kept_scafs)
print(len(all_kept_scafs))
print(len(all_kept_genes))
print(line_N)

count_file.close()


### takes fasta file, unwraps it, and adds seqs to a dict
def fasta_to_dict(in_fasta_file_name):
	output_fasta_name = in_fasta_file_name + ".TEMP_extract_fasta_file" 
	
	output_file = open(output_fasta_name, "w")
	print("\nUnwrapping fasta file")
	count = 0
	in_file = open(in_fasta_file_name)
	for line in in_file:
		count = count + 1
		line = line.rstrip("\n")
		if line.startswith(">") and count == 1:
			output_file.write(line + "\n")
		elif line.startswith(">") and count > 1:
			output_file.write("\n" + line + "\n")
		else: 
			output_file.write(line)	
	
	output_file.close()
	
	
	### add seqs to dictionary
	name_list = []
	seq_list = []
	seq_dict = {}
	
	done = 0
	seq_file_1 = open(output_fasta_name)
	for line in seq_file_1:
		lineA = line.rstrip("\n")
		if lineA.startswith(">"):
			lineB = lineA.lstrip(">")
			name_list.append(lineB)
		else:
			seq_list.append(lineA)
			done = done + 1
			seq_len = len(lineA)
	
	for element in range(0,len(name_list)):
		name1 = name_list[element]
		seq1 = seq_list[element].replace(" ", "") ## remove gaps if seq comes from gblocks 
		seq_dict[name1] = seq1

	## tidyup
	seq_file_1.close()
	os.remove(output_fasta_name)
	
	print("Read " + str(done) + " sequences from " + in_fasta_file_name)
	
	return(seq_dict)


iso_seq_dict = fasta_to_dict(isoform_file_name)

seq_list = []
for s in iso_seq_dict:
	seq_list.append(s)

seq_list = sorted(seq_list)

outfile = open(isoform_file_name.replace(".fa", "") + "_onscafs.fa", "w")

N_output_seq = 0
all_output_genes = set()
for s in seq_list:
	seq = iso_seq_dict.get(s)
	gene_name_c = s.split(".")[0]
	#print(gene_name_c)
	if gene_name_c in all_kept_genes:
		outfile.write(">" + gene_name_c + "\n" + seq + "\n")
		all_output_genes.add(gene_name_c)
		N_output_seq = N_output_seq + 1

### check excluded genes		
excluded_genes = all_kept_genes.difference(all_output_genes) ## in a1 not in a2
n_ncRNA = 0
n_tRNA = 0
n_rRNA = 0
for g in excluded_genes:
	if g.startswith("ncRNA_"):
		n_ncRNA = n_ncRNA + 1
	elif g.startswith("tRNA_"):
		n_tRNA = n_tRNA + 1
	elif g.startswith("rRNA_"):
		n_rRNA = n_rRNA + 1
	else:
		print("these genes should be in the output! ERROR!")
		print(g)


print("number of seqs kept")		








# 
# ### read in file
# isoform_file = open(isoform_file_name)
# for line in in_file:
# 	line = line.rstrip("\n")
# 	print(line)
# 
# in_file.close()

