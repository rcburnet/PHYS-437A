import numpy as np
from sklearn.neighbors import BallTree

#This script will read my list of 315 primaries and compare it to Ryan's list
# to make sure all of Ryan's primaries are in my list (which they are).
# Identical to script in Assignment_3.


#read files
file_284_primary = open('CasJobs_315_primaries_in_SDSS.txt')
list_284_primary = file_284_primary.readlines()
file_284_primary.close()

file_274_primary = open('Ryans_list_of_274_primaries.txt')
list_274_primary = file_274_primary.readlines()
file_274_primary.close()

#organize lists
for i in range(len(list_284_primary)):
    list_284_primary[i] = list_284_primary[i].split(',')

for i in range(len(list_274_primary)):
    list_274_primary[i] = list_274_primary[i].split('\t')

#find extras (find every entry in my list that isn't in Ryan's list)
length = 0
for i in range(len(list_284_primary)):
    count = 0
    for j in range(len(list_274_primary)):
        if list_284_primary[i][0] == list_274_primary[j][0]:
            count += 1
    if count == 0:
        print list_284_primary[i]
        length += 1

print length

#I get back 23 primaries that were in my list but not Ryan's. 23 extra, not 10.

#Compare Ryan's list to mine. Which primaries are in Ryan's list, but not my list.
for i in range(len(list_274_primary)):
    count = 0
    for j in range(len(list_284_primary)):
        if list_274_primary[i][0] == list_284_primary[j][0]:
            count += 1
    #if count == 0:
    #    print list_274_primary[i]

#I get back 13 primaries that were in Ryan's list but not mine. 13 primaries
# I cut that Ryan did not.

#This means I cut 13 that I shouldn't have cut, leaving me with 297 primaries,
# not 284 primaries after applying the second cut, and 23 primaries extra which
# are in badly masked regions or regions of incomplete coverage (297 - 23 = 274)
