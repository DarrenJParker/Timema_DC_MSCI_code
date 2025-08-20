import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:c:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_file_name = None
chr_to_invert = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** default_py.py | Written by DJP, 02/04/20 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Template script")	
		print("\n**** Usage****\n")
		print("python3 default_py.py -i [in file] -o [outprefix] \n\n")
		sys.exit(2)
	elif opt in ('-i'):
		in_file_name  = arg
	elif opt in ('-c'):
		chr_to_invert = arg
	else:
		print("i dont know")
		sys.exit(2)

	

### read in file
max_want_scaf = 0
in_file = open(in_file_name)
for line in in_file:
	line = line.rstrip("\n").split("\t")
	if line[0] == chr_to_invert:
		if int(line[1]) > max_want_scaf:
			max_want_scaf = int(line[1])
		if int(line[2]) > max_want_scaf:
			max_want_scaf = int(line[2]) 

in_file.close()

max_want_scaf = max_want_scaf + 1

print(max_want_scaf)

out_f = open(in_file_name.replace(".bed", "_" + chr_to_invert + ".bed"), "w")
in_file = open(in_file_name)
for line in in_file:
	line = line.rstrip("\n").split("\t")
	if line[0] == chr_to_invert:
		new_start = max_want_scaf - int(line[2])
		new_end   = max_want_scaf - int(line[1])
		# print(line)
		# print(line[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + line[3])
		out_f.write(line[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + line[3] + "\n")
	else:
		out_f.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\n")
in_file.close()


