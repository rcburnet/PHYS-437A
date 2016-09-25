#Query NED. Astroquery can also query other databases (such as SIMBAD, NASA ADS, etc) in case you need to query other databases for data. Look up astroquery for more details.
'''
from astroquery.ned import Ned
result_table = Ned.query_refcode('2011MNRAS.413..813C')
print result_table

from astropy.cosmology import WMAP9
print WMAP9.luminosity_distance(0.006181)
'''

