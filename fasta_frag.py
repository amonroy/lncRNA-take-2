###practice on HSR-omega, which is located on 3R
### this script is designed to pull out specific starts and stop on chromosome and the 
### sequence in between

import sys

#"~/Desktop/Projects/RNA_annotation/FlyBase/Dyak/dyak_3R_prac.fa"

chrom_file = sys.argv[1]
# chromosome_seqs


# ex/ 7.1_sort_prac.txt


###Here we are going to read in the sorted 7.1 file. Right now I am starting with '3R' for simplicity --> got from copy paste of sorted
out = open ('8_3R_fafrag_dyak_prac.fa', 'w')

with open('7.1_sort_prac.txt', 'r') as s: 
	# this list will be a list of tuple coordinates
	list = []
	for line in s:
		pre_coord = line.strip().rsplit()
		#print "This: %s should be greater than this %s" %(pre_coord[2], pre_coord[3])
		name = pre_coord[0]
		x = pre_coord[2]
		y = pre_coord[3]
		if x < y:
			list.append((name,x,y))
		elif y < x:
			list.append((name,y,x))
		else:
			print "This is screwy", x, y
		
	
for tuple in list:
	dmel_name = tuple[0]
	end_start = tuple[1]
	end_stop = tuple[2]		
	start = int(end_start)-1  
	stop = int(end_stop)-1

	with open (chrom_file, 'r') as i:
	#Change this to make sure name changes with fly species
		header = i.readline().rstrip()
		for line in i:
			#header = line
			#print header
			sequence = line.strip('\n')
			out.write(">%s\tstart:%s\tstop: %s\n" %(dmel_name, int(end_start), int(end_stop)))
			#string = sequence[start:stop]
			out.write(sequence[start:stop]+'\n')
			#out.write('%s\n' %(string))
out.close
		#out.write(">%s\tstart:%s\tstop: %s\n" %(dmel_name, int(end_start), int(end_stop)))
		#string = sequence[start:stop]
		#out.write("%s\n" %(string))
			
	
	