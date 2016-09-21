from astroquery.ned import Ned
result_table = Ned.query_refcode('2011MNRAS.413..813C')
print result_table
