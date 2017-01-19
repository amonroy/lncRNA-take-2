def transcript_location(some_file):
    """This function reads in transcript file of interest and grabs name of transcript and the lowest and highest location and the scaffold"""
    import re
    with open(some_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                pre_data = line.strip().split(';')
#                print "pre_data", pre_data
                data = []
                for l in pre_data:
                    m= l.strip(' ')
                    data.append(m)
#                print "This is data", data
                id = re.findall('>(\S+)\s', data[0]) #\S is any non-whitespace, and \s is any whitespace
                print "This is id", id
                chrom = re.findall('(loc=\S+):', data[1])
                print "this is chrom", chrom
#                print "this is data[1]", data[1]
                location = re.findall('loc=\S+:(\S+)', data[1])
#                print "This is location", location
                match = re.search('\(', location[0])

                if match:
                    print "MATCH", location[0]
                    inside = re.findall('\D+\((\S+)\)', location[0])
                    print "inside", inside
                    for i in inside:
                        no_comma = i.split(',')
                        final_list = []
                        for j in no_comma:
                            loc_list = j.split('..')
                            for k in loc_list:
                                final_list.append(k)
                    print "final list", final_list
                    final_list.sort()
                    print "sorted", final_list
                else:
                    loc = location[0].split('..')                                                                                                                                                     
                    print "no match", loc
                    loc.sort()
                    print "sorted", loc
                                        #8834315..8834682                                                                                                                                                   
#                                        loc.sort()                                                                                                                                                         
#                                        print "sorted", loc                                                                                                                                                





##					MATCH join(8781076..8781272,8783212..8784440,8784501..8784970)
#					join_loc = re.findall('join\((\S+)\)', location[0])
					#print join_loc
					#print type(join_loc) # list
					#loc_list =[]
#					###theres a problem here .. not getting all of the items in the list
#					for i in join_loc:
#					#	print "i before split", i
#						no_comma = i.split(',')
#						#print "This is no comma", no_comma
#					final_list = []
#					for j in no_comma:
#						loc_list = j.split('..')
						#print "This is loc_list", loc_list
#						for k in loc_list:
#							final_list.append(k)
#					print "this is final_list", final_list					
#                                        final_list.sort()
#                                       print "sorted", final_list
#               elif cmatch:
#                   print "cmatch", location[0]
                    #>FBtr0347178 type=ncRNA', ' loc=3L:complement(20207866..20208801)
#               else:
					#print "NOOOOOOO", location[0]
#            loc = location[0].split('..')
#					print loc
					#8834315..8834682
#                                        loc.sort()
#                                        print "sorted", loc
