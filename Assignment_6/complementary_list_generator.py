#This script generates complementary list to Ryan's list (galaxies that have
# neighbours within 1.5 mpc, but are not near Virgo, Leo, or Coma and aren't
# Andromeda).

import numpy as np
from sklearn.neighbors import BallTree

#new_file is comma delimited parent sample list
new_file = open('new_parent_sample_with_DA.txt')
primary_list = new_file.readlines()
new_file.close()
for i in range(len(primary_list)):
    primary_list[i] = primary_list[i].split(',')

#list_356 is ryan's list after first cut
file_356 = open('356_primaries_after_first_cut.txt')
list_356 = file_356.readlines()
new_file.close()
for i in range(len(list_356)):
    list_356[i] = list_356[i].split(',')

#complementary list of Ryan's list (871 - 356 - 1 = 514 primaries)
complementary_list = []
for i in range(len(primary_list)):
    if primary_list[i] not in list_356:
        if primary_list[i][0] != 'NGC0224':
            complementary_list.append(primary_list[i])

def ang_sep(array1, array2):
    '''
    2 lists of [R.A., Dec., optional: distance] -> angle (deg)
    Inputs two lists which correspond to two galaxies coordinates and
    calculates the angular separation of the two galaxies in degrees.
    '''
    RA1 = array1[0]
    DEC1= array1[1]
    RA2 = array2[0]
    DEC2 = array2[1]
    RA1 = float(RA1) * 2*np.pi/360.0
    DEC1 = float(DEC1) * 2*np.pi/360.0
    RA2 = float(RA2) * 2*np.pi/360.0
    DEC2 = float(DEC2) * 2*np.pi/360.0
    phi = np.arccos(np.sin(DEC1)*np.sin(DEC2) + np.cos(DEC1)*np.cos(DEC2)*np.cos(RA1-RA2))
    return phi*180.0/np.pi

def D_dist(array1, array2):
    '''
    2 list of [R.A., Dec., distance] -> distance (Mpc)
    Inputs two lists which correspond to two galaxies coordinates and
    calculates the 3D distance between the two galaxies in Mpc.
    '''
    d1 = array1[2]
    d2 = array2[2]
    phi = ang_sep(array1, array2)*np.pi/180.0
    d = np.sqrt(d1**2.0 + d2**2.0 - 2*d1*d2*np.cos(phi))
    return d

#formula to calculate angular separation of 1.0 Mpc at distance of specific
# primary in arcminutes.
def DA(d):
    d = float(d)
    theta = np.arctan(1.0/d)*10800/np.pi
    return theta

#Second cut is applied in CasJobs. Let's move on to third cut.

#Now apply third cut

#Cut all galaxies that are within 5 degrees of the center of the Virgo cluster
# (R.A., Dec.) = (186.75, 12.7167) or 3 degrees of the center of the Coma
# (R.A., Dec.) = (194.9542, 27.9806) or Leo (R.A., Dec.) = (176.1542, 19.7589)
# clusters

third_cut_list = [] #list of indices of galaxies that are to be cut
virgo_num = 0 #number of galaxies within 5 degrees of Virgo cluster
coma_num = 0 #number of galaxies within 3 degrees of Coma cluster
leo_num = 0 #number of galacies within 3 degrees of Leo cluster
for i in range(len(complementary_list)):
    phi_virgo = ang_sep([complementary_list[i][1],complementary_list[i][2]],[186.75, 12.7167])
    phi_coma = ang_sep([complementary_list[i][1],complementary_list[i][2]],[194.9542, 27.9806])
    phi_leo = ang_sep([complementary_list[i][1],complementary_list[i][2]],[176.1542, 19.7589])
    if phi_virgo < 5.0:
        virgo_num += 1
        third_cut_list.append(i)
    if phi_coma < 3.0:
        coma_num += 1
        third_cut_list.append(i)
    if phi_leo < 3.0:
        leo_num += 1
        third_cut_list.append(i)

#Need to actually apply the cut (remove 93 primaries that are near clusters).
ind_list4 = sorted(third_cut_list)[::-1]
for i in ind_list4:
    del complementary_list[i]

#This leaves me with complementary list of Ryan's list. Note: I still need to
# apply second cut. I send resulting text file to CasJobs to query it to
# determine if a primary is in SDSS. Complementary list is now 421 primaries long

#Write primary_list to new text file.
complementary_list_third_cut = []
for i in range(len(complementary_list)):
    complementary_list_third_cut.append(','.join(complementary_list[i]))

new_primary_list_file = open('421_complementary_list.txt', 'w')
for i in range(len(complementary_list_third_cut)):
    new_primary_list_file.write(complementary_list_third_cut[i])

new_primary_list_file.close()

