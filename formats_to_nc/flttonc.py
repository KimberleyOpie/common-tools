# REFORMAT_BOM_ASCII.PY
#
# Master script to split a BoM ArcView ASCII grid file into its component
# parts (header, data, tail metadata).
#
# The header defines elements of the data array and can be optionally written
# to an ASCII file.
#
# The data is converted to a numpy array, resampled on to a standard (SILO)
# grid, precision of values reduced to 1 decimal place (for consistency), and
# the missing value is redefined to -999.0.
# The numpy array is optionally written to both a flat binary f32 file,
# and/or to a CF and ACDD compliant netCDF file.
#
# The tail metadata is written to the netCDF file as the history attribute.
# Optionally, it may also be written to an ASCII file.


import re
import sys
import os.path
import numpy as np
import netCDF3
import nc3_handler as nh
import numpy_routines as nr
from datetime import datetime


def flttonc(fltstem, ncfile, varname, yyyymmdd, units=None):
    """Main function to process a binary flt and hdr with name FLTSTEM
       to netCDF file NCFILE with variable name VARNAME.  Time is set to
       the date corresponding to YYYYMMDD.  Units can be provided optionally.
    """

    d = datetime.strptime(yyyymmdd, "%Y%M%d")

    # Read the metadata from the header
    meta = {}
    r = re.compile("^(\S+)\s+(\S+)$")
    f = open(fltstem+'.hdr','r')
    for line in f:
        line=line.strip()
        m = re.match(r, line)
        if m:
            if m.group(1).lower() == 'byteorder':
                meta[m.group(1).lower()] = m.group(2)
            else:
                meta[m.group(1).lower()] = float(m.group(2))
    f.close()

    # Get the data and shape it appropriately
    a=np.fromfile(fltstem+'.flt',dtype=np.float32)
    a.reshape(int(meta['ncols']), int(meta['nrows']), 1)

    # Make the lon and lat coordinate variables
    if 'xllcenter' in meta:
        lon = nr.create_vector(meta['xllcenter'], \
                               meta['ncols'], \
                               meta['cellsize'])
    else:
        lon = nr.create_vector(meta['xllcorner']+0.5*meta['cellsize'], \
                               meta['ncols'], \
                               meta['cellsize'])
    if 'yllcenter' in meta:
        lat = nr.create_vector(meta['yllcenter'], \
                               meta['nrows'], \
                               meta['cellsize'])
    else:
        lat = nr.create_vector(meta['yllcorner']+0.5*meta['cellsize'], \
                               meta['nrows'], \
                               meta['cellsize'])
    # Reverse the latitude elements so they run from North to South
    lat = lat[::-1] 

    # Make a :history global attribute
    history = []
    now = datetime.now()
    history.append("flttonc at %04d-%02d-%02dT%02d:%02d:%02d" % (now.year, now.month, now.day, now.hour, now.minute, now.second))
    history.append("Input file: "+fltstem+".flt")

    attr = {}
    if 'nodata_value' in meta:
        attr[varname+':_FillValue'] = meta['nodata_value']
    if units:
        attr[varname+':units'] = units
    attr['history'] = '\n'.join(history)

    # Write the netCDF file
    ncobj = nh.nc3_open(ncfile,'w')
    nh.nc3_set_timelatlon(ncobj,1,len(lat),len(lon))
    nh.nc3_set_var(ncobj,varname)
    nh.nc3_set_var(ncobj,'wgs84',dims=())  # Grid mapping container
    nh.nc3_add_time(ncobj,[d])
    nh.nc3_add_data(ncobj,'latitude',lat)
    nh.nc3_add_data(ncobj,'longitude',lon)
    nh.nc3_add_data(ncobj,varname,a)
    nh.nc3_set_attributes(ncobj,attr)
    nh.nc3_close(ncobj)
    return ncfile


if __name__ == '__main__':

    if len(sys.argv) != 5 and len(sys.argv) != 6:
        print "Usage: python flttonc.py FLTSTEM NCFILE VARNAME YYYMMDD"
        sys.exit(1)

    if len(sys.argv) == 5:
        flttonc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    if len(sys.argv) == 6:
        flttonc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
