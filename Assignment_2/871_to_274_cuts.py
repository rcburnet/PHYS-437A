import numpy as np
from sklearn.neighbors import BallTree

#Note: scikit-learn must be v0.16 or less. v0.17 or higher results in
# deprecation error for passing array to ball tree.

#The below was used to read the parent sample from Atlas3D, which came
# in an ugly space delimited table, and write a new text file that
# is easily readable in python. Commented out because it has already been
# done.

######################################################################
#primary_file = open('list_of_871_primaries.txt')
#new_file = open('new_parent_sample.txt','w')

#primary_list = primary_file.readlines()

##for loop to remove whitespaces, comma deliminate columns, and write
## new primary list that is easily readible in python
#for i in range(len(primary_list)):
#    clean_string = ' '.join(primary_list[i].split())
#    primary_list[i] = clean_string.replace(' ',',')
#    new_file.write(primary_list[i]+'\n')
    
#new_file.close()
######################################################################

#new_file is comma delimited parent sample list
new_file = open('new_parent_sample.txt')
primary_list = new_file.readlines()
for i in range(len(primary_list)):
    primary_list[i] = primary_list[i].split(',')

#Note: Columns are: Galaxy name, RA (deg), DEC(deg), SBF, NED-D, Virgo,
# VHel (km/s), D (Mpc), M_K (mag), A_B (mag), T-type, and log(Re) (")
# in that order. For the cuts that Ryan applied, we care about RA, DEC,
# and D.

#Functions to calculate 3D distance between two points using Haversine
# formula to calculate angular separation, then cosine law to find
# distance separation.
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

#coord is list of (RA,DEC,D) coordinates of primaries, 871 points in 3D.
coord = []
for i in range(len(primary_list)):
    coord.append([float(primary_list[i][1]),float(primary_list[i][2]),float(primary_list[i][7])])

#Use Ball Tree with the above Haversine function as metric to find 
# nearest neighbors quickly. See sklearn documentation for details.
tree = BallTree(coord, leaf_size=2, metric=D_dist)

ind_list = [] #list of indices of the nearest neighbors
dist_list = [] #list of distances (rad) of nearest neighbors
for i in range(len(coord)):
    dist, ind = tree.query(coord[i], k = 2) #k=2 since first "nearest neighbor" is sometimes itself.
    #sometimes, D_dist will return NaN when calculating D_dist for same
    # point. In this case, ind[0][0] will be the nearest neighbor index.
    # Otherwise ind[0][0] is i, and so ind[0][1] would be the nearest neighbor
    # index
    if ind[0][0] == i:
        ind_list.append(ind[0][1])
        dist_list.append(dist[0][1])
    else:
        ind_list.append(ind[0][0])
        dist_list.append(dist[0][0])

#Note: ind_list is list of indices of i's nearest neighbour. Eg. ind_list[0]=382
# means that primary_list[0]'s nearest neighbor is primary_list[382]. dist_list
# is list of distances to i's nearest neighbor. Eg. dist_list[0] is 0.7381...
# meaning that primary_list[0] is 0.7381... Mpc from its nearest neighbor,
# primary_list[382]

#Apply first cut: Remove members that are within 1.5 Mpc of another member
first_cut_list = [] #list that will hold indices of primaries to be cut
for i in range(len(ind_list)):
    if dist_list[i] <= 1.5:
        first_cut_list.append(i)    #append index of primary
        first_cut_list.append(ind_list[i])  #append index of primary's nearest neighbor

#Remove duplicates from list        
first_cut_list = list(set(first_cut_list))

#Now first_cut_list is list of indices which correspond to primary_list indices
# that are primaries who are within 1.5 Mpc of another primary. These galaxies
# are to be cut.

#Note: first_cut_list is 514 elements long. This leaves us with 357 primaries
# left after the first cut, NOT 356 primaries as detailed in Ryan's paper. This
# is because we still need to remove Andromeda (NGC0224) from the list:
for i in range(len(primary_list)):
    if primary_list[i][0] == 'NGC0224':
        first_cut_list.append(i)

first_cut_list = sorted(first_cut_list)[::-1] #list is reversed so I can delete
                                              # elements in reverse order

#Now first_cut_list is 515 elements long with Andromeda, leaving us with 356
# primaries after the first cut like in Ryan's paper.

#Apply the cut
for i in first_cut_list:
    del primary_list[i]
    del coord[i]

#Now primary_list contains the 356 primaries that weren't cut

#create new text file for new primary list after first cut:
primary_list_first_cut = []
for i in range(len(primary_list)):
    primary_list_first_cut.append(','.join(primary_list[i]))

new_primary_list_file = open('356_primaries_after_first_cut.txt', 'w')
for i in range(len(primary_list_first_cut)):
    new_primary_list_file.write(primary_list_first_cut[i])

new_primary_list_file.close()



#Now apply second cut.

#First, remove objects that are not within SDSS survey area.

#Convert celestial coordinates of primaries to Survey coordinate system, then
# compare primaries Survey coordinates to stripe regions to see if galaxies are
# within survey coverage.
ind_list2 = [] #second index list corresponding to indices of primaries to be cut
survey_coord = [] #list of coordinates in survey coordinate system (lambda, eta)
for i in range(len(primary_list)):
    lambda_coord = np.arcsin(-1*np.cos((coord[i][0]-95)*np.pi/180.0)*np.cos(coord[i][1]*np.pi/180.0))*180.0/np.pi
    eta_coord = np.arccos(np.sin((coord[i][0]-95)*np.pi/180.0)*np.cos(coord[i][1]*np.pi/180.0)/np.cos(lambda_coord*np.pi/180.0))*180.0/np.pi-32.5
    survey_coord.append([lambda_coord, eta_coord])
    #Northern Galactic Cap:
    if -53.75 > eta_coord > -56.25:
        if not 7.1 > lambda_coord > -35.5:
            ind_list2.append(i)
    elif -51.25 > eta_coord > -53.25:
        if not 19.8 > lambda_coord > -42.8:
            ind_list2.append(i)
    elif -48.75 > eta_coord > -51.25:
        if not 28.3 > lambda_coord > -47.2:
            ind_list2.append(i)
    elif -46.25 > eta_coord > -48.75:
        if not 34.7 > lambda_coord > -50.4:
            ind_list2.append(i)
    elif -43.75 > eta_coord > -46.25:
        if not 39.6 > lambda_coord > -52.8:
            ind_list2.append(i)
    elif -41.25 > eta_coord > -43.75:
        if not 43.6 > lambda_coord > -54.6:
            ind_list2.append(i)
    elif -38.75 > eta_coord > -41.25:
        if not 46.8 > lambda_coord > -56.1:
            ind_list2.append(i)
    elif -36.25 > eta_coord > -38.75:
        if not 49.4 > lambda_coord > -57.6:
            ind_list2.append(i)
    elif -33.75 > eta_coord > -36.25:
        if not 51.7 > lambda_coord > -58.8:
            ind_list2.append(i)
    elif -31.25 > eta_coord > -33.75:
        if not 53.6 > lambda_coord > -59.6:
            ind_list2.append(i)
    elif -28.75 > eta_coord > -31.25:
        if not 55.2 > lambda_coord > -60.4:
            ind_list2.append(i)
    elif -26.25 > eta_coord > -28.75:
        if not 56.6 > lambda_coord > -61.2:
            ind_list2.append(i)
    elif -23.75 > eta_coord > -26.25:
        if not 57.8 > lambda_coord > -61.9:
            ind_list2.append(i)
    elif -21.25 > eta_coord > -23.75:
        if not 58.9 > lambda_coord > -62.4:
            ind_list2.append(i)
    elif -18.75 > eta_coord > -21.25:
        if not 59.8 > lambda_coord > -62.8:
            ind_list2.append(i)
    elif -16.25 > eta_coord > -18.75:
        if not 60.6 > lambda_coord > -63.1:
            ind_list2.append(i)
    elif -13.75 > eta_coord > -16.25:
        if not 61.2 > lambda_coord > -63.4:
            ind_list2.append(i)
    elif -11.25 > eta_coord > -13.75:
        if not 61.8 > lambda_coord > -63.6:
            ind_list2.append(i)
    elif -8.75 > eta_coord > -11.25:
        if not 62.3 > lambda_coord > -63.7:
            ind_list2.append(i)
    elif -6.25 > eta_coord > -8.75:
        if not 62.7 > lambda_coord > -63.8:
            ind_list2.append(i)
    elif -3.75 > eta_coord > -6.25:
        if not 63.1 > lambda_coord > -63.7:
            ind_list2.append(i)
    elif -1.25 > eta_coord > -3.75:
        if not 63.3 > lambda_coord > -63.7:
            ind_list2.append(i)
    elif 1.25 > eta_coord > -1.25:
        if not 63.5 > lambda_coord > -63.5:
            ind_list2.append(i)
    elif 3.75 > eta_coord > 1.25:
        if not 63.7 > lambda_coord > -63.3:
            ind_list2.append(i)
    elif 6.25 > eta_coord > 3.75:
        if not 63.7 > lambda_coord > -63.1:
            ind_list2.append(i)
    elif 8.75 > eta_coord > 6.25:
        if not 63.8 > lambda_coord > -62.7:
            ind_list2.append(i)
    elif 11.25 > eta_coord > 8.75:
        if not 63.7 > lambda_coord > -62.3:
            ind_list2.append(i)
    elif 13.75 > eta_coord > 11.25:
        if not 63.6 > lambda_coord > -61.8:
            ind_list2.append(i)
    elif 16.25 > eta_coord > 13.75:
        if not 63.4 > lambda_coord > -61.2:
            ind_list2.append(i)
    elif 18.75 > eta_coord > 16.25:
        if not 63.1 > lambda_coord > -60.6:
            ind_list2.append(i)
    elif 21.25 > eta_coord > 18.75:
        if not 62.8 > lambda_coord > -59.8:
            ind_list2.append(i)
    elif 23.75 > eta_coord > 21.25:
        if not 62.4 > lambda_coord > -58.9:
            ind_list2.append(i)
    elif 26.25 > eta_coord > 23.75:
        if not 61.9 > lambda_coord > -57.8:
            ind_list2.append(i)
    elif 28.75 > eta_coord > 26.25:
        if not 61.2 > lambda_coord > -56.6:
            ind_list2.append(i)
    elif 31.25 > eta_coord > 28.75:
        if not 60.4 > lambda_coord > -55.2:
            ind_list2.append(i)
    elif 33.75 > eta_coord > 31.25:
        if not 59.6 > lambda_coord > -53.6:
            ind_list2.append(i)
    elif 36.25 > eta_coord > 33.75:
        if not 58.8 > lambda_coord > -51.7:
            ind_list2.append(i)
    elif 38.75 > eta_coord > 36.25:
        if not 57.6 > lambda_coord > -49.4:
            ind_list2.append(i)
    elif 41.25 > eta_coord > 38.75:
        if not 56.1 > lambda_coord > -46.8:
            ind_list2.append(i)
    elif 43.75 > eta_coord > 41.25:
        if not 54.6 > lambda_coord > -43.6:
            ind_list2.append(i)
    elif 46.25 > eta_coord > 43.75:
        if not 52.8 > lambda_coord > -39.6:
            ind_list2.append(i)
    elif 48.75 > eta_coord > 46.25:
        if not 50.4 > lambda_coord > -34.7:
            ind_list2.append(i)
    elif 51.25 > eta_coord > 48.75:
        if not 47.2 > lambda_coord > -28.3:
            ind_list2.append(i)
    elif 53.75 > eta_coord > 51.25:
        if not 42.8 > lambda_coord > -19.8:
            ind_list2.append(i)
    elif 56.25 > eta_coord > 53.75:
        if not 35.5 > lambda_coord > -7.1:
            ind_list2.append(i)
            
    #Southern Galactic Cap (don't have data for):
    #elif ...

    #Otherwise, galaxies are not within survey area if their eta_coord don't
    # correspond to a stripe's eta_coord range:
    else:
        ind_list2.append(i)

#Need to actually apply the cut. Note: I did not apply the cut yet as I could
# not successfully find the galaxies to cut. After applying cut, need to write
# left over primaries (after cut) to new text file.



#Now apply third cut

#Cut all galaxies that are within 5 degrees of the center of the Virgo cluster
# (R.A., Dec.) = (186.75, 12.7167) or 3 degrees of the center of the Coma
# (R.A., Dec.) = (194.9542, 27.9806) or Leo (R.A., Dec.) = (176.1542, 19.7589)
# clusters

third_cut_list = [] #list of indices of galaxies that are to be cut
virgo_num = 0 #number of galaxies within 5 degrees of Virgo cluster
coma_num = 0 #number of galaxies within 3 degrees of Coma cluster
leo_num = 0 #number of galacies within 3 degrees of Leo cluster
for i in range(len(coord)):
    phi_virgo = ang_sep(coord[i],[186.75, 12.7167])
    phi_coma = ang_sep(coord[i],[194.9542, 27.9806])
    phi_leo = ang_sep(coord[i],[176.1542, 19.7589])
    if phi_virgo < 5.0:
        virgo_num += 1
        third_cut_list.append(i)
    if phi_coma < 3.0:
        coma_num += 1
        third_cut_list.append(i)
    if phi_leo < 3.0:
        leo_num += 1
        third_cut_list.append(i)

#Need to actually apply the cut. Note: I did not apply the cut yet as I could
# not successfully apply the second cut above. After applying cut, need to write
# left over primaries (after cut) to new text file.
