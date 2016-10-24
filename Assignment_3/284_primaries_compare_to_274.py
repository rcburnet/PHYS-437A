import numpy as np
from sklearn.neighbors import BallTree

#This script will read my list of 284 primaries that I got from applying Ryan's
# cuts (except for cuts to galaxies in badly masked regions or regions of
# incomplete coverage) to Ryan's list of 274 primaries, to find the extras in
# my list and see why Ryan cut them. They are most likely the primaries
# in badly masked regions or regions of incomplete coverage in the SDSS. The
# script also find primaries that were in Ryan's list but not my final list
# (ie. primaries I cut that Ryan did not). It is found that I have 23 extra
# primaries in my list, not the expected 10, meaning 23 primaries are in badly
# masked regions or regions of incomplete coverage. I know for a fact that these
# 23 extra are a result of being in badly masked regions or regions of incomplete
# coverage because I checked last week's script to make sure the first and
# third cuts cut primaries that weren't in Ryan's list (ie. the first and
# third cuts were being applied properly). It is also found that 13 primaries
# in Ryan's list were not in my list. This is a result of the CrossID tool not
# being sufficient at providing the primaries that are in SDSS. I cut 13 that
# I should not have cut as they were in SDSS (I checked that they were using
# the navigate tool:
# http://skyserver.sdss.org/dr8/en/tools/chart/navi.asp )


#read files
file_284_primary = open('284_primaries_after_third_cut.txt')
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
for i in range(len(list_284_primary)):
    count = 0
    for j in range(len(list_274_primary)):
        if list_284_primary[i][0] == list_274_primary[j][0]:
            count += 1
    #if count == 0:
    #    print list_284_primary[i]

#I get back 23 primaries that were in my list but not Ryan's. 23 extra, not 10.

#Compare Ryan's list to mine. Which primaries are in Ryan's list, but not my list.
for i in range(len(list_274_primary)):
    count = 0
    for j in range(len(list_284_primary)):
        if list_274_primary[i][0] == list_284_primary[j][0]:
            count += 1
    if count == 0:
        print list_274_primary[i]

#I get back 13 primaries that were in Ryan's list but not mine. 13 primaries
# I cut that Ryan did not.

#This means I cut 13 that I shouldn't have cut, leaving me with 297 primaries,
# not 284 primaries after applying the second cut, and 23 primaries extra which
# are in badly masked regions or regions of incomplete coverage (297 - 23 = 274)
