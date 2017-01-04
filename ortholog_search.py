## This is the script I am writing to search for lncRNA orthologs based on synteny (location in the genome

# how to run:
#(command)
# dependencies

import sys
import collections
import itertools
import re
import gc

rna_blast_file = sys.argv[2]
protein_blast_file = sys.argv[3]

#sys.argv[1] = modified gff file (modified to only have gene and ncRNA)
## where/name
## made by what script?
#sys.argv[2] = ncRNA blast file of species of interest
## where/name
#sys.argv[3] = protein blast file
## where/name
#sys.argv[4] = protein file
## where/name

window_length = 4 ### this can change, theoretically. This is telling the script to find 4 upstream genes and 4 downstream genes

