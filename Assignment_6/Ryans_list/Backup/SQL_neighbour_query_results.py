##DEPRECATED. Kept in case I need to go back to individual SQL queries. Currently
# using single SQL query to find nearby neighbours for all primaries.

#Script will first create directories for each primary of Ryan's 356 primaries
# that survived the first cut and are also in SDSS (315 in total). The list of
# primaries are in the text file called CasJobs_315_primaries_in_SDSS.txt as
# well as my CasJobs DB.

#Afterwards, the script will print the corresponding SQL query that I will need
# to send to CasJobs for each of the 315 primaries to find the neigbours within
# 2 degrees of the primary and save them as a text file in their corresponding
# directory.

#First: make directories

import os, sys
import numpy as np

#open text file and generate list of 315 primaries
primaries_file = open('CasJobs_315_primaries_in_SDSS.txt')
primaries_list = primaries_file.readlines()

primaries_list = primaries_list[1:] #remove first row (column name row)


#create directories if they don't exist and generate individual SQL CasJobs
# query to find objects within 2 degrees of primary for each primary.
combined_query = '' #combined query of all primaries
for i in range(len(primaries_list)):
    primaries_list[i] = primaries_list[i].split(',')
    primary_name = primaries_list[i][0]
    RA = primaries_list[i][1]
    DEC = primaries_list[i][2]
    #create directories if they don't exist
    if not os.path.exists('./'+primary_name):
        os.makedirs('./'+primary_name)
    #generate individual SQL queries to find nearby objects
    write_to_file = ''
    write_to_file +='SELECT\n'
    write_to_file += ' p1.objid, p1.z, p1.zErr, p2.ra, p2.dec, p2.clean, p2.flags, p2.type, p2.expRad_u, p2.expRad_g, p2.expRad_r, p2.expRad_i, p2.expRad_z, p2.expRadErr_u, p2.expRadErr_g, p2.expRadErr_r, p2.expRadErr_i, p2.expRadErr_z, p2.cmodelMag_u, p2.cmodelMag_g, p2.cmodelMag_r, p2.cmodelMag_i, p2.cmodelMag_z, p2.cmodelMagErr_u, p2.cmodelMagErr_g, p2.cmodelMagErr_r, p2.cmodelMagErr_i, p2.cmodelMagErr_z, p2.extinction_u, p2.extinction_g, p2.extinction_r, p2.extinction_i, p2.extinction_z\n'
    write_to_file += 'FROM PhotoObj AS p2\n'
    write_to_file += ' JOIN Photoz AS p1 ON p1.objid = p2.objid\n'
    write_to_file += 'JOIN dbo.fGetNearbyObjEq('+str(RA)+','+str(DEC)+',120) AS p0 ON p2.objid = p0.objid WHERE p1.z < 0.15 AND p1.z > 0'
    file_to_write_to = open('./'+primary_name+'/'+primary_name+'_SQL_query.txt', 'w')
    file_to_write_to.write(write_to_file)
    file_to_write_to.close()
    #add write_to_file string to combined_query
    if i != 0:
        combined_query += '\n'
        combined_query += '\n'
    combined_query += write_to_file

#write combined_query to file
file_to_write_to = open('combined_SQL_query.txt','w')
file_to_write_to.write(combined_query)
file_to_write_to.close()
#SQL query will look like:

'''
SELECT 
  p1.objid, p1.z, p1.zErr, p2.ra, p2.dec, p2.clean, p2.?ags, p2.type, p2.expRad_u, p2.expRad_g, p2.expRad_r, 
  p2.expRad_i, p2.expRad_z, p2.expRadErr_u, p2.expRadErr_g, p2.expRadErr_r, p2.expRadErr_i, p2.expRadErr_z, p2.cmodelMag_u, 
  p2.cmodelMag_g, p2.cmodelMag_r, p2.cmodelMag_i, p2.cmodelMag_z, p2.cmodelMagErr_u, p2.cmodelMagErr_g, p2.cmodelMagErr_r, 
  p2.cmodelMagErr_i, p2.cmodelMagErr_z, p2.extinction_u, p2.extinction_g, p2.extinction_r, p2.extinction_i, p2.extinction_z
FROM PhotoObj AS p2 
  JOIN Photoz AS p1 ON p1.objid = p2.objid 
JOIN dbo.fGetNearbyObjEq(RA,DEC,120) AS p0 ON p2.objid = p0.objid WHERE p1.z < 0.15 AND p1.z > 0
'''

#You can also make a query to count number of objects (ie. number of rows from
# output of above query) :

'''
SELECT count(*) as 'total',
  sum( case when (p2.type=3) then 1 else 0 end) as 'Galaxies',
  sum( case when (p2.type=6) then 1 else 0 end) as 'Stars',
  sum( case when (p2.type not in (3,6)) then 1 else 0 end) as 'Other'
FROM PhotoPrimary AS p2
 JOIN Photoz AS p1 ON p1.objid = p2.objid
JOIN dbo.fGetNearbyObjEq(175.077,9.009861,120) AS p0 ON p2.objid = p0.objid WHERE p1.z < 0.15 AND p1.z > 0
'''



