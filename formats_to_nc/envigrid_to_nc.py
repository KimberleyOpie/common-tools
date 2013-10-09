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

import sys, os, re
import numpy as np
import hashlib
from datetime import timedelta
from datetime import datetime
try: from collections import OrderedDict
except ImportError:
    try: from ordereddict import OrderedDict
    except ImportError:
        print "Require OrderedDict, https://pypi.python.org/pypi/ordereddict"
        raise
import netcdf_builder as nb
import json_handler as jh
import numpy_routines as nr

def envihdr_to_dict(header):
    """
    Read an ENVI header file into an OrderedDict. All keys are returned as
    lowercase. All values are returned as strings or a list of strings. Values
    in braces, {}, are multi-line joined if required and the braces are removed.
    The description value is split on whitespace and reloined with a single
    whitespace between words to form a string. Other values in braces are split
    on comma to form a list.
    """
    Other 
    meta = OrderedDict()
    r = re.compile("^\s*(\S+)\s*=\s*(\S+)$")
    f = open(header,'r')
    for line in f:
        line = line.strip()
        m = re.match(r, line)
        if m:
            meta[m.group(1).lower()] = m.group(2)
        elif len(meta) > 0:
            k = meta.keys()[-1]
            meta[k] += ' '+line
    f.close()
    # Remove braces, extra whitespace and split on comma
    r = re.compile("^\s*{\s*(.+)\s*}\s*$")
    for k,v in meta.iteritems():
        m = re.match(r,v)
        if m:
            if k == 'description':
                # No splitting required, just remove extra whitespace
                meta[k] = ' '.join(v.split())
            else:
                meta[k] = re.split("\s*,\s*",v)
    return meta

def archdr_to_dict(header):
    """
    Read an Arc float header file into an OrderedDict. All keys are returned as
    lowercase. All values are returned as strings or a list of strings.
    """
    # Read the metadata from the header
    meta = OrderedDict()
    r = re.compile("^\s*(\S+)\s+(\S+)\s*$")
    f = open(header,'r')
    for line in f:
        line = line.strip()
        m = re.match(r, line)
        if m:
            meta[m.group(1).lower()] = m.group(2)
            # Old method. New method retains strings by default
            #if m.group(1).lower() == 'byteorder':
            #    meta[m.group(1).lower()] = m.group(2)
            #else:
            #    meta[m.group(1).lower()] = float(m.group(2))
    f.close()
    return meta

def make_latlon_vectors_geo(meta):
    """
    Create latitude and longitude vectors from the values contained in the
    ENVI header dictionary, meta. Returned as a tuple of (latitude,longitude).
    """
    #samples = 9580
    #lines   = 7451
    #map info = {Projection name,
    #            reference pixel x location (in file coordinates),
    #            reference pixel y location (in file coordinates),
    #            pixel x location (in map coordinates),
    #            pixel y location (in map coordinates),
    #            x pixel size,
    #            y pixel size,
    #            projection zone,
    #            units=Units (map coordinate unites),
    #            North or South (UTM only) }
    #map info = {Geographic Lat/Lon, 1.0000, 1.0000, 110.00000000, -10.00000000, 4.6974247390e-003, 4.6974247390e-003, WGS-84, units=Degrees}
    #map info = {UTM, 1, 1, 295380.000, 4763640.000, 30.000000, 30.000000, 13, North}
    

def make_latlon_vectors(meta):
    """
    Create latitude and longitude vectors from the values contained in the
    Arc header dictionary, meta. Returned as a tuple of (latitude,longitude).
    """
    # Convert corners to centers if required
    if 'xllcorner' in meta:
        meta['xllcenter'] = float(meta['xllcorner'])+float(meta['cellsize'])/2
    if 'yllcorner' in meta:
        meta['yllcenter'] = float(meta['yllcorner'])+float(meta['cellsize'])/2

    # Create vectors
    input_lon = nr.create_vector(float(meta['xllcenter']), \
                                 float(meta['ncols']), \
                                 float(meta['cellsize']))
    input_lat = nr.create_vector(float(meta['yllcenter']), \
                                 float(meta['nrows']), \
                                 float(meta['cellsize']))
    input_lat = input_lat[::-1]    # reverse elements
    return input_lat,input_lon

def dtype_envi_to_numpy(dt):
    """
    Return the corresponding numpy data type for a given ENVI data type.
    """
    envi2np = {}
    # 1=8-bit byte;
    # 2=16-bit signed integer;
    # 3=32-bit signed long integer;
    # 4=32-bit floating point;
    # 5=64-bit double-precision floating point;
    # 6=2x32-bit complex, real-imaginary pair of double precision;
    # 9=2x64-bit double-precision complex, real-imaginary pair of double precision;
    # 12=16-bit unsigned integer;
    # 13=32-bit unsigned long integer;
    # 14=64-bit signed long integer;
    # 15=64-bit unsigned long integer.


def flt_to_numpy(flt, meta):
    dtype = np.float32
    if 


    # Copy data into a numpy array and shape it appropriately
    input_data = np.fromfile(flt,dtype=np.float32)
    input_data.reshape(int(meta['ncols']), int(meta['nrows']), 1)

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

def flttonc(fltstem, ncfile, varname, yyyymmdd, units=None):
    """Main function to process a binary flt and hdr with name FLTSTEM
       to netCDF file NCFILE with variable name VARNAME.  Time is set to
       the date corresponding to YYYYMMDD.  Units can be provided optionally.
    """
    pass
    d = datetime.strptime(yyyymmdd, "%Y%M%d")


if __name__ == '__main__':

    if len(sys.argv) != 5 and len(sys.argv) != 6:
        print "Usage: python flttonc.py FLTSTEM NCFILE VARNAME YYYMMDD"
        sys.exit(1)

    if len(sys.argv) == 5:
        flttonc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    if len(sys.argv) == 6:
        flttonc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
