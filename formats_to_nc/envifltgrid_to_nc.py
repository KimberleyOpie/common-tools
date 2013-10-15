# ENVIFLTGRID_TO_NC.PY
#
# Author:
# Matt Paget, CSIRO Marine and Atmospheric Research, Canberra, Australia.
#
# Available from:
# https://github.com/cmar-rs/common-tools/formats_to_nc/envifltgrid_to_nc.py
#
# License:
# Apache License, 2.0 http://www.apache.org/licenses/LICENSE-2.0
#
# Version:
# 1.0   15 October 2013
#       Adapted from arcfltgrid_to_nc.py
#       Not working code yet.
#
# Basically the header file is parsed to a dictionary and  the float file is
# parsed to a numpy array. A netcdf file is built from the coordinate variables
# and the data array (using netcdf_builder.py).
#
# The time dimension is generally parsed from the filename as its not readily
# provided in the header information and is rarely added as a parsable element
# in the file. As such the parsing regular expression and time period value
# may need to be adjusted for your application.
#
# Extra hooks are provided to capture per-file metadata in JSON format files.
# The JSON metadata can be edited, collated and added to the netcdf file as
# attributes as a seperate process (see json_handler.py).
#
# Requires:
#  NumPy
#  OrderedDict
#  https://github.com/cmar-rs/common-lib/python/netcdf_builder.py
# Optional:
#  https://github.com/cmar-rs/common-lib/python/numpy_routines.py
#  https://github.com/cmar-rs/common-lib/python/json_handler.py

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
