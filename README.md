# lncRNA-take-2
#being watched:
#ortholog_search.py
#parse_gff3.py


#parse_gff3.py is explained in "ortholog_search.py"; its needs to be run on the gff3 file before ortholog_search is run

# Blasts need to be run before running "ortholog_search.py"
# shell script used to do this (run on longleaf) is now being watched
# query ncRNA dmel against subject non-mel transcripts, blastn
# a possibility is also to have CDS against CDS, but that should be encompassed in "gene_orthologs_..." file from flybase