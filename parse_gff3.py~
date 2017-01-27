import sys
import datetime

#argv[1] should be the gff3 file, latest release "dmel-all-r6.11.gff.gz" 
#on ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/
#input is gff file
#output is gff-like file, but with only the lines concerning 'ncRNA' or 'protein'

def parse_gff_nc_genes():
	"""This function takes the Dmel gff file from argv[1] and parses it. 
	It returns lines with 'ncRNA' and 'genes' in them. And removes the sequences.
	It uses sys and datetime"""
	
	gff3 = sys.argv[1]

	today = datetime.date.today()

	ncRNA_gene_out = open('1_dmel_protein_ncRNA_%s_out.txt' %today, 'w')

	f= open(gff3, 'r')

	while True:
		line = f.readline()
		if not line.startswith('#'):
			data = line.strip().split("\t")
	
			if True == True:
	
				try:
		#originally had this as 'gene', will now make it 'protein'
					if data[2] == 'protein':
						ncRNA_gene_out.write('%s' %line)
			
					elif data[2] == 'ncRNA':
						ncRNA_gene_out.write('%s' %line)
				
					else:
						continue
				
				except IndexError:
					if line.startswith('>'):
						break
					elif line.startswith('A'):
						break
					elif line.startswith('T'):
						break
					elif line.startswith('G'):
						break
					elif line.startswith('C'):
						break
				
			if not line:
				break

	ncRNA_gene_out.close()			
	f.close()			

parse_gff_nc_genes()