import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'f:s:x:hSR')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_file_name = None
replace_internal_stop  = False
remove_name_after_space = False
stop_codon = "."
rename_all_seqs_prefix = False

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** fasta_prot_trim_stop.py | Written by DJP, 06/01/22 in Python 3.5 in Porto, Portugal ****\n")
		print("Removes stop codons from fasta files. Default to only remove terminal 3' stops but can be extended to anywhere in the sequence with the -R option")	
		print("\n**** Usage ****\n")
		print("python3 fasta_prot_trim_stop.py -f [in fasta file] [options] \n\n")
		print("**** Options ****\n")
		print("-s\tstop symbol. Default = '.'")		
		print("-S\tSimple header in output - remove everything after a space in fasta header. Default: OFF")
		print("-X\trename all seqs with specified prefix. Default: OFF")		
		print("-R\tRemove internal stops. Default: OFF\n\n\n")
		
		sys.exit(2)
	elif opt in ('-f'):
		in_file_name = arg
	elif opt in ('-R'):
		 replace_internal_stop  = True
	elif opt in ('-S'):
		 remove_name_after_space = True
	elif opt in ('-x'):
		 rename_all_seqs_prefix = arg
	elif opt in ('-s'):
		 stop_codon = arg
	else:
		print("i dont know")
		sys.exit(2)

if in_file_name == None:
	print("\n\nERROR. No input file\n\n")

print("\nstop codon symbol is " + stop_codon + " use -s to alter\n")
stop_codon

if replace_internal_stop == True:
	print("\n-R specified. Internal stop codons will be replaced\n")

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


prot_seq_dict = fasta_to_dict(in_file_name)
N_terminal_stops = 0
seqs_with_internal_stops = 0
outfile = open(in_file_name.replace(".fasta", "").replace(".fa", "") + "_aans.fa", "w")

seen_names = set()

S_N = 0
for s in prot_seq_dict:
	seq = prot_seq_dict.get(s)
	S_N = S_N + 1
	if seq.endswith(stop_codon) :
		seq = seq.rstrip(stop_codon)
		N_terminal_stops = N_terminal_stops + 1
	if stop_codon in seq:
		seqs_with_internal_stops = seqs_with_internal_stops + 1
		if replace_internal_stop == True:
			seq = seq.replace(stop_codon, "")		
		else:
			print("internal stop codon in seq " + s + ". Not removed. To remove use -R")



	if remove_name_after_space == True:
		s = s.split(" ")[0]
	
	if s not in seen_names:
		seen_names.add(s)
	else:
		print("\n\nError - duplicate sequence names found\n\n")
		outfile.close()
		os.remove(in_file_name.replace(".fasta", "").replace(".fa", "") + "_aans.fa")
		sys.exit()
	
	if rename_all_seqs_prefix == False:
		outfile.write(">" + s + "\n" + seq + "\n")
	else:
		outfile.write(">" + rename_all_seqs_prefix + str(S_N) + "\n" + seq + "\n")
	
	
print("Number of seqs with stop codon aa at the end: " + str(N_terminal_stops))
print("Number of seqs with internal stop codon aa: " + str(seqs_with_internal_stops))


print("\n\nFinished, John Redlantern\n\n\n")









