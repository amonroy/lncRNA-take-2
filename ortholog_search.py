## This is the script I am writing to search for lncRNA orthologs based on synteny (location in the genome
### Jan 31 : so things to think about ... the end ortholog bit seems to be grabbing some weird stuff, see what you do at the point where there may be more than one protein at a start, and see if there's a way to manage that. Also, not all proteins will be in the ortholog part, so you have to do a try/except for that part.



# how to run:
#(command)
# dependencies

import sys
import collections
import itertools
import re
import gc



#sys.argv[1] = modified gff file (modified to only have gene and ncRNA)
## where/name out-take-2/1_dmel_protein_ncRNA_2017-01-04_out.txt
## made by what script? parse_gff3.py (in same directory as this script)#sys.argv[2]= protein fasta file
#sys.argv[3]= ortholog map file
##
##sys.argv[2] = ncRNA blast file of species of interest # this is wrong
## where/name
#sys.argv[3] = protein blast file #this is wrong
## where/name
#sys.argv[4] = protein file # this has to change
#rna_blast_file = sys.argv[2] # this has to change
#protein_blast_file = sys.argv[3]
## where/name

window_length = 4 ### this can change, theoretically. This is telling the script to find 4 upstream genes and 4 downstream genes

def mel_gff_list():
    """This function takes the modified gff3 file (sys.argv[1]) and creates a list. This should have all dmel proteins and all dmel lncRNAs"""
    mod_gff3 = sys.argv[1]
    with open(mod_gff3, 'r') as f:
        gff = [line.strip().split('\t') for line in f]
    f.close() # this close in unnecessary, I think, with the way I opened the file
    return gff

def mel_ncRNA_strand(list): #2 function
    """This function takes the gff_list and makes an ncRNA_list """
    strand_dic = {} #initiates dictionary
    for i in list:
    	#print i
    	#quit()
        if i[2] == 'ncRNA':
            idRNA = i[8].split(';')[0].split('=')[1]
            strand = i[6]
            #start = i[3]
            if int(i[3]) > int(i[4]):
                print "THIS DID HAPPEN!!"
            	start = i[4]
            if int(i[4]) > int(i[3]):
            	start = i[3]
            	
            chrom = i[0]
            strand_dic[idRNA] = [strand, start,chrom]
    return strand_dic
#'FBtr0335437': ['-', '11977640', '3R'], 'FBtr0332375': ['-', '15295636', '2L']

def make_start_chrom_dictionary(list):
    """This function makes a dictionary (named for the chromosome arm)  using the start location (the base pair number) as a key and the name of the protein as the entry. Uses mel gff list. As of right now drops out weird chromosomes"""
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
        	continue
            #print "Don't care =", i[0]
    #print chrom_2L
    return chrom_2L, chrom_2R, chrom_3L, chrom_3R, chrom_4, chrom_X, chrom_Y


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

def indexing_location(rna_strand_dict, dict_chrom, p2g_dict): # I did something different in another script with the protein... I used blast hits... not sure why I did that, but something to think about.
    """This function makes a dictionary with upstream genes and downstream genes for each melanogaster lncRNA, I sometimes refer to it as the front-back dictionary. This function has to be run for each chromosome and/or chromosome arm."""
    fbgn_id_dict = {}
    chrom = dict_chrom["chrom_name"]# there's a "chrom_name" key in each dictionary, to make life easier
    del dict_chrom["chrom_name"] #made a copy of this and then deleted it. This is so I can repeat this for each chromosome dictionary and grab the name of the chromosome with out too much trouble.
    list_chrom = dict_chrom.keys() #making an ordered list of the key in location_dic. The keys are the start of the gene
    numbers_chrom = [int(x) for x in list_chrom]
    sorted_chrom = sorted(numbers_chrom)
    #print sorted_chrom
    #print "sorted_chrom length:", len(sorted_chrom)
    for k,v in rna_strand_dict.iteritems():
        #print "This is k, v:", k, v 
        rna_chrom = v[2]
        start = v[1]
        #print "this is chrom:", chrom
        if rna_chrom == chrom:
            #print "%s == %s" %(rna_chrom, chrom)
            rna_index = sorted_chrom.index(int(start))
            up_counter = 1
            down_counter = 1
            uplist, downlist = [], []
            while len(uplist) < window_length:
                #print "rna_index:", rna_index
                try:
                    #print sorted_chrom[rna_index- up_counter]
                    #print "What is in loc dic", dict_chrom.get(str(sorted_chrom[rna_index - up_counter]))
                    #print "What is in v", v
                    for i in dict_chrom.get(str(sorted_chrom[rna_index- up_counter])):
                        #print "This is i", i
                        try:
                            p2g_dict[i]
                        except KeyError:
                            #print "This is probably and RNA", i
                            continue
                    #print "this shit is done"
                except IndexError:
                    #print "this is up_counter", up_counter
                    pass

                try:
                    for i in dict_chrom.get(str(sorted_chrom[rna_index - up_counter])):
                        if re.search("FBpp*", i):
                            try:
                                uplist.append(p2g_dict[i]) #maybe change this to i[0] later
                                break
                            
                            except KeyError:
                                print "this is causing trouble up here", i
                    up_counter += 1
                            
                except IndexError:
                    note = 'no_more_genes'
                    uplist.append(note)
                    #print "reached end of list"
                    break
            #print "This is uplist"
            #print uplist

            while len(downlist) < window_length:
                try:
                    for i in dict_chrom.get(str(sorted_chrom[rna_index + down_counter])):
                        if re.search("FBpp*", i):
                            try:
                                downlist.append(p2g_dict[i])
                                break
                            except KeyError:
                            	print "This causing trouble", i
                    down_counter += 1

                except IndexError:
                    note1 = 'no_more_genes'
                    downlist.append(note1)
                    #print "This number gave me a problem", rna_index + down_counter
                    break
            fbgn_id_dict[k]= uplist, downlist

   # print "THIS IS FB ME:", fbgn_id_dict
    #print len(fbgn_id_dict)
    return fbgn_id_dict

def match_pp_to_gn():
    """This reads in the protein file, and finds the gn name associated with the protein"""
    pp_to_gn_dict = {}
    protein = sys.argv[2]
    with open(protein, 'r') as p:
        for line in p:
            if line.startswith('>'):
                data = line.strip().split(';')
                #print data
                chrom = re.findall('loc=(\S+):', data[1])
                #print chrom
                #print data[2]
                pp = data[2].split('=')[1]
                gn = data[4].split(',')[0].split('=')[1]
                pp_to_gn_dict[pp]= [gn, chrom[0]]

            else:
                continue
   # print "P2G", pp_to_gn_dict
    return pp_to_gn_dict
    
mel_gff_obj = mel_gff_list()
mel_rna_strand_dic = mel_ncRNA_strand(mel_gff_obj)
dic_2L, dic_2R, dic_3L, dic_3R, dic_4, dic_X, dic_Y = make_start_chrom_dictionary(mel_gff_obj)       
testing_functions_1to3(mel_gff_obj, mel_rna_strand_dic, dic_2L, dic_2R, dic_3L, dic_3R, dic_4, dic_X, dic_Y)

protein_to_gene_dict = match_pp_to_gn()

up_down_protein_dic_2L = indexing_location(mel_rna_strand_dic, dic_2L, protein_to_gene_dict) # Have to do this for each chromosome/chromosome arm. a little tedious.


print "*" * 8
print "finished up_down_protein_dic_2L"
up_down_protein_dic_2R= indexing_location(mel_rna_strand_dic, dic_2R, protein_to_gene_dict)
up_down_protein_dic_3L= indexing_location(mel_rna_strand_dic, dic_3L, protein_to_gene_dict)
up_down_protein_dic_3R= indexing_location(mel_rna_strand_dic, dic_3R, protein_to_gene_dict)
up_down_protein_dic_4= indexing_location(mel_rna_strand_dic, dic_4, protein_to_gene_dict)
up_down_protein_dic_X= indexing_location(mel_rna_strand_dic, dic_X, protein_to_gene_dict)
up_down_protein_dic_Y= indexing_location(mel_rna_strand_dic, dic_Y, protein_to_gene_dict)

set_2L = set(up_down_protein_dic_2L.iterkeys())
set_2R = set(up_down_protein_dic_2R.iterkeys())
set_3L =set(up_down_protein_dic_3L.iterkeys())
set_3R = set(up_down_protein_dic_3R.iterkeys())
set_4 = set(up_down_protein_dic_4.iterkeys())
set_X= set(up_down_protein_dic_X.iterkeys())
set_Y= set(up_down_protein_dic_Y.iterkeys())
u = set.intersection(set_2L, set_2R, set_3L, set_3R, set_4, set_X, set_Y) #just want to verify that there are no overlapping keys
print "This is intersections:",  u


dicts = up_down_protein_dic_2L, up_down_protein_dic_2R, up_down_protein_dic_3L, up_down_protein_dic_3R, up_down_protein_dic_4, up_down_protein_dic_X, up_down_protein_dic_Y
super_dict = {}
for d in dicts:
    for k, v in d.iteritems():  # d.items() in Python 3+
        super_dict[k] = v

#print "This is super_dict", super_dict
#'FBtr0343761': ([['FBgn0000479', 'X'], ['FBgn0000479', 'X'], ['FBgn0003374', 'X'], ['FBgn0003086', 'X']], [['FBgn0000479', 'X'], ['FBgn0000479', 'X'], ['FBgn0029656', 'X'], ['FBgn0000479', 'X']]), 'FBtr0343762': ([['FBgn0003867', '3R'], ['FBgn0038903', '3R'], ['FBgn0038902', '3R'], ['FBgn0027575', '3R']], [['FBgn0013759', '3R'], ['FBgn0013759', '3R'], ['FBgn0041229', '3R'], ['FBgn0045470', '3R']]
#quit()
#print "This is length", len(super_dict) # 2899


def find_best_blast_hit(blast_file):
    """This reads in your protein blast output (in format 6) and gives you the data for the best hit per protein. Makes a dictionary. Important: I am filtering out proteins and lncRNAs orthologs that have more than one 'best hit'"""
    #print "In find_best_hit()"
    import collections
    best_hit_dict = {}
    with open(blast_file, 'rU') as f:
        for line in f:
            data = line.strip().split('\t')
            key = data[0]
            details = [data[1], data[8], data[9], data[10]]
    #print "these are details", details
            if key in best_hit_dict:
                if float(best_hit_dict[key][0][3]) < float(data[10]):
                    continue
                if float(best_hit_dict[key][0][3])> float(data[10]):
                    del best_hit_dict[key]
                    best_hit_dict.setdefault(key, []).append(details)
            #protein_dict[key]= [data[1], data[8], data[9], data[10]]
                if float(best_hit_dict[key][0][3]) == float(data[10]):
                    best_hit_dict.setdefault(key, []).append(details)
            else:
                best_hit_dict.setdefault(key, []).append(details)
        #protein_dict[key]= [data[1], data[8], data[9], data[10]]
        #print protein_dict
        a = []
        for k,v in best_hit_dict.iteritems():
            if len(v) == 1:
                a.append(k)
            bestest_hit_dict = {}
            for i in a:
                bestest_hit_dict[i] = best_hit_dict[i][0]

                    #print bestest_hit_dict
                    #print len(best_hit_dict)
                    #print len(bestest_hit_dict) (rna 2837: 2575 ; protein 30337: 17263)
            print bestest_hit_dict
                    #'FBtr0334812': ['Scf_2L', '2081514', '2080524', '0.0'], 'FBtr0334815': ['Scf_3R', '17902553', '17901770', '0.0'], 'FBtr0334816': ['Scf_2L', '7123428', '7123185', '1.15e-117'], 'FBtr0332776': ['Scf_3L', '19790127', '19789715', '3.48e-160'], 'FBtr0344464': ['Scf_X', '20090299', '20089726', '0.0'], 'FBtr0347118': ['Scf_3L', '3148141', '3148551', '0.0'], 
                    #'FBpp0298339': ['Scf_2R', '13592458', '13592060', '1.66e-70'], 'FBpp0112465': ['Scf_NODE_3978', '2785', '2522', '1.11e-52'], 'FBpp0297425': ['Scf_3L', '19713597', '19711975', '1.88e-157'], 'FBpp0086930': ['Scf_2R', '9424819', '9424103', '1.87e-145']
            return bestest_hit_dict
                #'FBtr0342626': ['Scf_2L', '12192374', '12192839', '0.0'], 
                #'FBtr0342625': ['Scf_X', '12617057', '12616737', '1.11e-121'], 
                #'FBtr0339184': ['Scf_3R', '25355617', '25355121', '0.0'], 
                #'FBtr0339187': ['Scf_3R', '26860087', '26860485', '7.04e-146'], 
                #'FBtr0339186': ['Scf_3R', '26860019', '26860485', '0.0'],
                
                #'FBpp0078599': ['Scf_3R', '104719', '103040', '0.0'], 
                #'FBpp0074019': ['Scf_X', '15309593', '15309057', '5.91e-114'], 
                #'FBpp0071254': ['Scf_X', '8703331', '8704488', '0.0'], 
                #'FBpp0071255': ['Scf_X', '8731616', '8732224', '1.09e-67'], 
                #'FBpp0300231': ['Scf_X', '9387607', '9387074', '3.19e-77'], 

#need to write function that identifies id (FBtr0350598 and links it to "loc = Scf_X:join etc". I'll probably simplify that by grabbing the smallest and largest locations
#>FBtr0350598 type=mRNA; loc=Scf_X:join(2802014..2802920,2814693..2814783,2822502..2822764,2822821..2823025,2825516..2825646,2829718..2832641,2832711..2835932,2835999..2836308,2836384..2836818,2836898..2837513); ID=FBtr0350598; name=Dsim\N-RB; dbxref=FlyBase_Annotation_IDs:Dsim\N-RB,FlyBase_Annotation_IDs:GD16631-RB,FlyBase:FBtr0350598,GNOMON:Dsim_gnomon_101_rna.16769621; MD5=d91a94abf402b28a62307616efffae81; length=9104; parent=FBgn0012846; release=r2.02; species=Dsim; 
#GAAACAGATCGCTTTTTTCCAGTGGACGAAACGGTTGTGAAAGCGGACGAGCGTAAGGCAGACGAACCTGGAAAGCGCAG
#AGCACAGTTCTCAACGTTTTTTTTTTTGGAAGTGAGTGCAACAACGCACGCAAACCGCGCAGCCAACAGGATATACAAAC
#AAATCAATCACAGCAAGCAAATGCCATGAAATGAAAAGGATGGCCCCAGCGGGAAAGCCGTTCAGCAAGAGCAAGGAGTG
#CCTGTCGAAGGGATAGCAACGAGAGAGAGAGAGAGAGGAAGAGAGAAACAAGGATTTTCGAAAAGTGTATCTACCTCGAG
#TCGCGCGTGTGTGAGAGTGAGACGCAAGCCGAGTGCAAGAAGCGCAATACGCAAGCGTGCGGCGTCGGTTTGAATTTGAA
#TTTGTGCTCGATCCTCGCGAAGAGAAAAGCAAGCAAAAGATACACGAAAAGCGTTCTTTTTTTGCCACTTTTTTTTTATG

                
def transcript_location(some_file):
    """This function reads in transcript file of interest and grabs name of transcript and the lowest and highest location and the scaffold"""
    #need to sort this lists to get the lowest and highest
    import re
    with open(some_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                data = line.strip().split(';')
                #print "This is data", data
                id = re.findall('>(\S+)\w', data[0])
                #print "This is id", id
                chrom = re.findall('\W(loc=\S+):', data[1])
                #might need to do some fancy stuff here. Like, if "join" is in location, do such and such", else" this other thing"                                                                         
                location = re.findall('\Wloc=\S+:(\S+)', data[1])
                #print "This is location", location
                match = re.search('join', location[0])
                if match:
                    
                    ##MATCH join(8781076..8781272,8783212..8784440,8784501..8784970)
                    join_loc = re.findall('join\((\S+)\)', location[0])
                    #print join_loc
                    #print type(join_loc) # list
                    #loc_list =[]
                    
                    for i in join_loc:
                    #print "i before split", i
                        no_comma = i.split(',')
                        #print "This is no comma", no_comma
                    final_list = []
                    for j in no_comma:
                        loc_list = j.split('..')
                        #print "This is loc_list", loc_list
                        for k in loc_list:
                            final_list.append(k)
                    print "this is final_list", final_list
            
                else:
                    #print "NOOOOOOO", location[0]
                    loc = location[0].split('..')
                    print loc
            
def mel_gene_set(dict): # this uses the flanking genes, specifically
    """This function finds unique mel genes, and puts them in a set (what is returned), so we don't get the same coords twice. It takes fbgn_id_dict. This is so we have the mel genes that we need coordinates for in the non-mel species', ie we're using this to find the orthologs that we care about"""
    mel_gene_set = set()
    for k, v in dict.iteritems():
        #v[0] is up, v[1] is down
        print "this is v:", v
        print "this is v[0]", v[0]
        for mg in v[0]:
            mel_gene_set.add(mg[0])
        for mj in v[1]:
            mel_gene_set.add(mj[0])
    return mel_gene_set

mel_set =  mel_gene_set(super_dict)
#print "This is mel_set", mel_set

def map_mel_gene_to_nonmel_ortho(set):
    ### good news: this appears to work. Something to consider, however.  Perhaps what you should do is have an additional point for "strandedness" (+ or -) and make sure that coord[0] is in fact the 'earliest base'
    """This function maps genes from mel_set to orthologs of nonmel species found in ortholog file, from flybase"""
   ### FBgn_ID[0]      GeneSymbol[1]      Arm/Scaffold[2]    Location[3]        Strand[4]  Ortholog_FBgn_ID[5]        Ortholog_GeneSymbol[6]     Ortholog_Arm/Scaffold[7]   Ortholog_Location[8]      Ortholog_Strand[9] OrthoDB_Group_ID[10]
#so far it looks like this is already written as coord[0] is always smaller than coord[1], so you just did that and don't need to use it.
    mapping = dict()
    fly = "Dsim" #this is hard coded, change that once you've finished debugging
    with open(sys.argv[3], 'r') as orthos:
        not_count = 0
        yes_count = 0
        for line in orthos:
            if not line.startswith('#') and not line.startswith('\n'):
                data = line.strip().split('\t')
                if fly in data[6]:
#                    print data
 #                   quit()
                    if data[0] in set:
                        coord = data[8].split("..")
                        direction =data[9]
                        ortho_group = data[10]
                        try:
                            if 'nonp' in mapping[data[0]]:
                                mapping[data[0]][data[5]] = [data[5], data[7], coord[0], coord[1], direction, ortho_group]
                        except KeyError:
                            mapping[data[0]] = {}
                            mapping[data[0]]['nonp'] = [data[5],data[7], coord[0], coord[1], direction, ortho_group]
                        yes_count = yes_count + 1
                    if data[0] not in set:
                        not_count = not_count +1 
    return mapping

orthos_gene_map = map_mel_gene_to_nonmel_ortho(mel_set)
#print "this is orthos_gene_map", orthos_gene_map
for sk, sv in super_dict.iteritems():
    print sk
    for si in sv[0]:
        print "upstream!"
        upstr = sv[0]
        for ui in upstr:
            print ui[0], ui[1]
            print orthos_gene_map[ui[0]]
           # print "that should have been the otrholog"
        print upstr
    for sj in sv[1]:
        print "downstream"
        dwstr = sv[1]
        print dwstr
        for di in dwstr:
            print di[0], di[1]
            print orthos_gene_map[ui[0]]
           # print "that should have been the downstream ortholog"

quit()
