def transcript_location(some_file):
    """This function reads in transcript file of interest and grabs name of transcript and the lowest and highest location and the scaffold"""
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
					#print "MATCH", location[0]
##					MATCH join(8781076..8781272,8783212..8784440,8784501..8784970)
					join_loc = re.findall('join\((\S+)\)', location[0])
					#print join_loc
					#print type(join_loc) # list
					#loc_list =[]
					###theres a problem here .. not getting all of the items in the list
					for i in join_loc:
					#	print "i before split", i
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
					#8834315..8834682