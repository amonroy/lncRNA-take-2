## This is the script I am writing to search for lncRNA orthologs based on synteny (location in the genome

# how to run:
#(command)
# dependencies

import sys
import collections
import itertools
import re
import gc

#rna_blast_file = sys.argv[2]
#protein_blast_file = sys.argv[3]

#sys.argv[1] = modified gff file (modified to only have gene and ncRNA)
## where/name out-take-2/1_dmel_protein_ncRNA_2017-01-04_out.txt
## made by what script? parse_gff3.py (in same directory as this script) 
#sys.argv[2] = ncRNA blast file of species of interest
## where/name
#sys.argv[3] = protein blast file
## where/name
#sys.argv[4] = protein file
## where/name

window_length = 4 ### this can change, theoretically. This is telling the script to find 4 upstream genes and 4 downstream genes

def mel_gff_list():
    """This function takes the modified gff3 file (sys.argv[1]) and creates a list"""
    mod_gff3 = sys.argv[1]
    with open(mod_gff3, 'r') as f:
        gff = [line.strip().split('\t') for line in f]
    f.close() # this close in unnecessary, I think, with the way I opened the file
    return gff

def mel_ncRNA_strand(list): #2 function
    """This function takes the gff_list and makes an ncRNA_list"""
    strand_dic = {} #initiates dictionary
    for i in list:
        if i[2] == 'ncRNA':
            idRNA = i[8].split(';')[0].split('=')[1]
            strand = i[6]
            start = i[3]
            chrom = i[0]
            strand_dic[idRNA] = [strand, start,chrom]
    return strand_dic
#'FBtr0335437': ['-', '11977640', '3R'], 'FBtr0332375': ['-', '15295636', '2L']

def make_start_chrom_dictionary(list):
    """This function makes a dictionary (named for the chromosome arm)  using the start location (the base pair number) as a key and the name of the protein as the entry."""
    #chrom_dict = collections.defaultdict(dict)
    chrom_list = ['2L','2R','3L','3R', '4', 'X','Y']
    chrom_2L = {}
    chrom_2L["chrom_name"]="2L" #make sure has chrom_name key in each dictionary
    #2L_chrom = {}
    chrom_2R = {}
    chrom_2R["chrom_name"]="2R"
    chrom_3R = {}
    chrom_3R["chrom_name"]="3R"
    chrom_3L = {}
    chrom_3L["chrom_name"]="3L"
    chrom_4 = {}
    chrom_4["chrom_name"]="4"
    chrom_X = {}
    chrom_X["chrom_name"]="X"
    chrom_Y = {}
    chrom_Y["chrom_name"] = "Y"

    #print 2L_chrom
    
    for i in list:
        start = i[3]
        name = i[8].split(';')[0].split('=')[1]
        if i[0] == '2L':
    #print "2L=", i[0]
            try:
                chrom_2L[start].append(name)
            except KeyError:
                chrom_2L[start] = [name,]
        if i[0] == '2R':
         #print "2R=", i[0]
            try:
                chrom_2R[start].append(name)
            except KeyError:
                chrom_2R[start]= [name,]
        if i[0] == '3L':
            #print "3L=:", i[0]
            try:
                chrom_3L[start].append(name)
            except KeyError:
                chrom_3L[start]= [name,]
        if i[0] == '3R':
        #print '3R=', i[0]
            try:
                chrom_3R[start].append(name)
            except KeyError:
                chrom_3R[start]= [name,]
        if i[0] == '4':
        #print '4', i[0]
            try:
                chrom_4[start].append(name)
            except KeyError:
                chrom_4[start]= [name,]
        if i[0] == 'X':
        #print "X=", i[0]
            try:
                chrom_X[start].append(name)
            except KeyError:
                chrom_X[start]= [name,]
        if i[0] == "Y":
        #print "Y=", i[0]
            try:
                chrom_Y[start].append(name)
            except KeyError:
                chrom_Y[start]= [name,]
        elif i[0] not in chrom_list: #this is taking out poorly mapped regions.
            print "Don't care =", i[0]
    #print chrom_2L
    return chrom_2L, chrom_2R, chrom_3L, chrom_3R, chrom_4, chrom_X, chrom_Y

mel_gff_obj = mel_gff_list()
mel_rna_strand_dic = mel_ncRNA_strand(mel_gff_obj)
dic_2L, dic_2R, dic_3L, dic_3R, dic_4, dic_X, dic_Y = make_start_chrom_dictionary(mel_gff_obj)

def testing_functions_1to3(o1, o2, o3_2L, o3_2R, o3_3L, o3_3R, o3_4, o3_X, o3_Y):
    """This is me attempting to write a simple test to show that everything is ok... hopefully I'll figure out how to do unit tests properly, b/c I think this is going to be pretty lame."""
    import random
    rand_smpl = [ o1[i] for i in sorted(random.sample(xrange(len(o1)), 20)) ]
    for j in rand_smpl:
        #print "This is output from mel_gff_list", j
        id = j[8].split(';')[0].split('=')[1]
        if j[2]== "ncRNA":
            print id, "This means function 1 looks alright"
            print o2[id], "This means function 2 looks alright"
            if o2[id][2]=="2L":
                print o3_2L[o2[id][1]], "This means function 3 looks alright"
            if o2[id][2]== "2R":
                print o3_2R[o2[id][1]], "This means function 3 looks alright"
            if o2[id][2]== "3L":
                print o3_3L[o2[id][1]], "This means function 3 looks alright"
            if o2[id][2]=="3R":
                print o3_3R[o2[id][1]], "This means function 3 looks alright"
            if o2[id][2]== "4":
                print o3_4[o2[id][1]], "This means function 3 looks alright"
            if o2[id][2]== "X":
                print o3_X[o2[id][1]], "This means function 3 looks alright"
            if o2[id][2]== "Y":
                print o3_Y[o2[id][1]], "This means function 3 looks alright"
            elif o2[id][2] not in ["2L", "2R", "3L", "3R","4", "X", "Y"] :
                print "something might have gone screwy", j
           #'FBtr0332375': ['-', '15295636', '2L']   
           # print o2
        
testing_functions_1to3(mel_gff_obj, mel_rna_strand_dic, dic_2L, dic_2R, dic_3L, dic_3R, dic_4, dic_X, dic_Y)
