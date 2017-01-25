def transcript_location(some_file):
    """This function reads in transcript file of interest and grabs name of transcript and the lowest and highest location and the scaffold"""
#there's something weird going on
    import re
    tscrpt_locs= {}
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
               # print "This is id", id
                chrom = re.findall('(loc=\S+):', data[1])
               # print "this is chrom", chrom
                final_chrom =chrom[0].split('=')[1]
               # print "This is final chrom", final_chrom
#                print "this is data[1]", data[1]
                location = re.findall('loc=\S+:(\S+)', data[1])
#                print "This is location", location
                match = re.search('\(', location[0])
                if match:
                   # print "MATCH", location[0]
                    inside = re.findall('\D+\((\S+)?\)', location[0])
                   # print "inside", inside
                    for i in inside:
                        no_comma = i.split(',')
                        final_list = []
                        for j in no_comma:
                            loc_list = j.split('..')
                            for k in loc_list:
                                final_list.append(k)
                   # print "final list", final_list
                    final_list.sort()
                    #print "sorted", final_list
                    #print "mfirst and last", final_list[0], final_list[-1]
                    tscrpt_locs[id[0]]=[final_chrom, final_list[0], final_list[-1]]
                else:
                    loc = location[0].split('..')                                                                                                                                                     
                    #print "no match", loc
                    loc.sort()
                    #print "sorted", loc
                    tscrpt_locs[id[0]]= [final_chrom, final_list[0], final_list[-1]]
                    #print "first and last", final_list[0], final_list[-1]
#    print tscrpt_locs
    return tscrpt_locs
