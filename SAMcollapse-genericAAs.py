#!/usr/bin/env python2

# Collapses SAM files

import sys

SAMList = list(line.strip().split("\t") for line in open(sys.argv[1]))
Output = open(sys.argv[2], "w")

myDict = dict()

for samLine in SAMList:
    line = "\t".join(samLine)
    if samLine[0].startswith("@"): #if a header
        Output.write(line + "\n")
    else:
        read = str(samLine[0])
        if read in myDict.keys():
            myDict[read].append(line)
        else:
            myList = [line]
            myDict[read] = (myList) 

tRNAgroupDict = dict()
counter = 0
counterGroup = 0
for k,v in myDict.iteritems():
    #print("\n", k)
    if v[0].split("\t")[2].startswith("chr"):
        feature = v[0].split("\t")[2].split("-")[1]
    else:
        feature = v[0].split("\t")[2]
    database_matches = list()
    for i in v:
        match = i.split("\t")[2]
        database_matches.append(match)
    if all(feature in x for x in database_matches) == True:  # if every value in the dict contains the feature species
        newAlignment = v[0].split(":")
        newAlignment[-1:] = ["1"] # Replace last element (NH flag DB match #) with 1
        newAlignment = ":".join(newAlignment)
        Output.write(newAlignment + "\n")     
        counter = counter + 1
    else:
        if v[0].split("\t")[2].startswith("chr"):
            tRNAgroup = (v[0].split("\t")[2].split("-")[1])[:3] # e.g. 'Pro' from a Proline tRNA
        else:
            tRNAgroup = v[0].split("\t")[2]
        if all(tRNAgroup in x for x in database_matches) == True:  # if every value in the dict contains the tRNA group
            counterGroup = counterGroup + 1
            for alignment in v:
                aligned = alignment.strip().split("\t")
                ### Determine read alignment orientation.
                ### This script should only have to deal with either '0' or '16' flags in the 2nd column of SAM
                ### as these indicate +/- orientation for primary alignment. I have added the 256/272 catchalls just in
                ### case something weird happened to the SAM file.
                ### This step is necessary to make sure the correct plot is generated
                if int(aligned[1]) == 0: # read aligned in + orientation
                    newLine = aligned
                    break
                elif int(aligned[1]) == 16: # read aligned in - orientation
                    newLine = aligned
                elif int(aligned[1]) == 256: # read aligned in + orientation (secondary alignment) 
                    newLine = aligned
                    break
                elif int(aligned[1]) == 272: # read aligned in - orientation (secondary alignment) 
                    newLine = aligned
                else:
                    print ("ERROR in SAMcollapse.py: Alignment below does not have correct flag:\n%s" % (aligned))
            #newLine = tRNAgroup + "\t" + newLine[3] + "\t" + str(len(newLine[9]) + int(newLine[3])) + "\t" + newLine[0]
            # Create new fake name for the DB match for this read mapping instance
            tRNAmatch = v[0].split("\t")[2].split("-")[1]
            new_tRNA_name = "chr100.trna1000-" + tRNAmatch[:3]
            newAlignment = v[0].split("\t")
            #print(newAlignment)
            newAlignment[2] = new_tRNA_name # Replace 3rd element in list with new tRNA name
            newAlignment[14] = "NH:i:1" # Replace no. of alignments with 1
            #print(newAlignment)
            #newAlignment = newAlignment.split(":")
            #print(newAlignment)
            #newAlignment[-1:] = ["1"] # Replace last element (NH flag DB match #) with 1
            newAlignment = "\t".join(newAlignment)
            #print(newAlignment)
            Output.write(newAlignment + "\n")     

            #print(newLine)
            #if tRNAgroup in tRNAgroupDict.keys():
            #    tRNAgroupDict[tRNAgroup].append(newLine)
            #else:
            #    myList = [newLine]
            #    tRNAgroupDict[tRNAgroup] = (myList)

Output.close()

#print counter, counterGroup

#leftovers = open(sys.argv[2] + "_tRNAs-almost-mapped.txt", "w")
#for k,v in tRNAgroupDict.iteritems():
#    for values in v:
#        leftovers.write(values + "\n")
#leftovers.close()

