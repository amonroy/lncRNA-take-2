## This script is involved, but basically read in ncRNA start and stop location, gets the chromosome (DNA) that has that location, and prints it out

import sys
from location_find import transcript_location
from chrom_functions import location_seq
from chrom_functions import read_chrom
import datetime

ncRNA_file = sys.argv[1]
chrom_file = sys.argv[2]

ncRNA_dict = transcript_location(ncRNA_file)

#print ncRNA_dict
#datetime.date.today()
today = datetime.date.today()

# '1_dmel_protein_ncRNA_%s_out.txt' %today, 'w')
chrom_set = set()
for k,v in ncRNA_dict.iteritems():
    chrom_set.add(v[0])
#print "chrom_set", chrom_set

f = open("ncRNA_DNA_dmel_seqs_%s_out.fasta" %today, "w")
f.close()

with open("ncRNA_DNA_dmel_seqs_%s_out.fasta" %today, "a") as g:
    #print "I'm in this loop"
    for i in chrom_set:
        #print "this is i", i
        full_chrom_seq = read_chrom(chrom_file, i)
        #print "full_chrom_seq", full_chrom_seq
        for key, value in ncRNA_dict.iteritems():
            if i == value[0]:
                #print "%s == %s" %(key, value[0])
                #print key, value
                ncRNA_seq = location_seq(full_chrom_seq, value[1], value[2])
                #print "ncRNA_seq", ncRNA_seq
                g.write(">"+key+";"+str(value[0])+";"+str(value[1])+";"+str(value[2])+"\n"+ncRNA_seq+"\n")
