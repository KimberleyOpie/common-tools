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
import datetime as dt
try: from collections import OrderedDict
except ImportError:
    try: from ordereddict import OrderedDict
    except ImportError:
        print "Require OrderedDict, https://pypi.python.org/pypi/ordereddict"
        raise
import netcdf_builder as nb
import json_handler as jh
import numpy_routines as nr

def arcflthdr_to_dict(header):
    """
    Read an Arc float header file into an OrderedDict. All keys are returned as
    lowercase. All values are returned as strings or a list of strings.
    """
    # Read the metadata from the header
    meta = OrderedDict()
    r = re.compile("^(\S+)\s+(.+)$")
    f = open(header,'r')
    for line in f:
        #line = line.strip()
        m = re.match(r, line.strip())
        if m:
            meta[m.group(1).lower()] = m.group(2)
    f.close()
    return meta

def set_latlon(meta, datadict=None):
    """
    Create latitude and longitude vectors from the values contained in the
    Arc header dictionary, meta. Returned as a tuple of
    (latitude,longitude,datadict).
    """
    # Convert corners to centers if required
    if 'xllcorner' in meta:
        meta['xllcenter'] = float(meta['xllcorner'])+float(meta['cellsize'])/2
    if 'yllcorner' in meta:
        meta['yllcenter'] = float(meta['yllcorner'])+float(meta['cellsize'])/2

    # Create vectors
    lonvec = nr.create_vector(float(meta['xllcenter']), \
                              float(meta['ncols']), \
                              float(meta['cellsize']))
    latvec = nr.create_vector(float(meta['yllcenter']), \
                              float(meta['nrows']), \
                              float(meta['cellsize']))
    latvec = latvec[::-1]    # reverse elements

    # Create/update a datadict, which includes standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({
        'xmin':min(lonvec),
        'xmax':max(lonvec),
        'xstep':meta['cellsize'],
        'xnum':meta['ncols'],
        'xunits':'degrees_east',
        'ymin':min(latvec),
        'ymax':max(latvec),
        'ystep':meta['cellsize'],
        'ynum':meta['nrows'],
        'yunits':'degrees_north'
    })

    # Return a tuple
    return latvec,lonvec,datadict

def fltgrid_to_numpy(fltfile, meta, datadict=None):
    """
    Convert a binary file to a 3D numpy array.
    Arg datadict is not used but is retained for consistency.
    """
    # Default datatype for Arc flt grids
    dtype = '<f4'  # little-endian float32

    # Consider byteorder
    #   lsbfirst -> little-endian
    #   msbfirst -> big-endian
    if 'byteorder' in meta:
        if re.match('[mM]',meta['byteorder']): dtype = '>f4'

    # Copy data into a 3D array
    # Want shape to align with time,latitude,longitude
    data = np.fromfile(fltfile,dtype=dtype)
    data.reshape(1, int(meta['nrows']), int(meta['ncols']))

    # Return a tuple
    return data,datadict

def resample_array(input_data, input_lat, input_lon, input_dict=None):
    """
    Code snippets for resampling the input array on to a new grid, changing
    the missing value or reducing the data type or floating point precision
    of the array.

    Edit the components of this routine for your particular needs. Just check
    that input_data and output_data arrays are being used and referenced as
    required.
    """
    pass

def set_datetime(fname, datadict=None):
    """
    Code example for creating the date/time information. Generally, the
    date/time information will be in the filename as its not included in the
    standard header information.

    Edit the components of this routine for your particular needs.
    """
    m = re.search('^(\d{4})(\d\d)',fname)
    d1year = int(m.group(1))
    d1month = int(m.group(2))
    d1day = 1

    # Monthly = (d1month + 1) - (1 day)
    d2year = d1year
    d2month = d1month + 1
    if d2month == 13:
        d2month = 1
        d2year = d2year + 1

    d1 = dt.datetime(d1year,d1month,d1day)
    d2 = dt.datetime(d2year,d2month,d1day) - dt.timedelta(days=1)

    # Create/update a datadict, which adds some standardised labels for later
    if datadict is None: datadict = dict()

    # Metadata elements
    # See http://en.wikipedia.org/wiki/ISO_8601
    datadict.update({
        'dcreate':"unknown",
        'dmodify':dt.datetime.utcnow().strftime("%Y%m%dT%H%M%S"),
        'tmin':d1.strftime("%Y-%m-%d"),
        'tmax':d2.strftime("%Y-%m-%d"),
        'tduration':"P1M",
        'tresolution':"P1M"
    })

    # Retrun as a tuple
    return d1,d2,datadict

def set_varname(fname, datadict=None):
    """
    Code example for creating the variable name. Generally, the variable name
    will be in the filename as its not included in the standard header
    information.

    Edit the components of this routine for your particular needs.
    """
    varname = "varname1"
    varlong = "long variable name"
    varunit = "unit"

    # Create/update a datadict, which adds some standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({'varname':varname,
                     'varlong':varlong,
                     'varunits':varunit})
    return datadict

def set_attributes(fname, meta, datadict=None):
    # Add a default no_data value if required. If you want the default value
    # set outside this routine then pass it in as meta['nodata_value'].
    miss = -999.0
    if 'nodata_value' in meta: miss = meta['nodata_value']

    # Create/update a datadict, which includes standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({
        'missing':miss
        'sha1_arcfloat':hashlib.sha1(open(fltfile,'r').read()).hexdigest()
    }) 

    # Make a :history global attribute
    history = []
    now = dt.datetime.now()
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
    d = dt.datetime.strptime(yyyymmdd, "%Y%M%d")


if __name__ == '__main__':

    if len(sys.argv) != 5 and len(sys.argv) != 6:
        print "Usage: python flttonc.py FLTSTEM NCFILE VARNAME YYYMMDD"
        sys.exit(1)

    if len(sys.argv) == 5:
        flttonc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    if len(sys.argv) == 6:
        flttonc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
