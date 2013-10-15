# BOMASCII_TO_CSIROFLT.PY
#
# Author:
# Matt Paget, CSIRO Marine and Atmospheric Research, Canberra, Australia.
#
# Available from:
# https://github.com/cmar-rs/common-tools/raw_data_tools/bomascii_to_csiroflt.py
#
# License:
# Apache License, 2.0 http://www.apache.org/licenses/LICENSE-2.0
#
# Version:
# 1.0   15 October 2013
#       Adapted from arcasciigrid_to_nc.py
#
# Basically the ascii text is parsed into a list of strings for the header and
# for the data. The header defines the coordinate variables and the missing
# data value. The data is converted to a numpy array. Float32 and header files
# are built from the coordinate variables and the data array.
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
import datetime as dt
import copy
try: from collections import OrderedDict
except ImportError:
    try: from ordereddict import OrderedDict
    except ImportError:
        print "Require OrderedDict, https://pypi.python.org/pypi/ordereddict"
        raise
import netcdf_builder as nb
import json_handler as jh
import numpy_routines as nr

def split_asciigrid(fname):
    """
    Split an Arc ASCII grid file into its header, data and (possibly) tail
    components.  End-of-line characters and whitespace at either end of a
    line are striped off.  The head and tail components are returned as a
    list of lines.  The data component is returned as a list of lists, with
    each inner list being a list of (string) values in a line.

    The number of rows and columns in the grid are determined from the
    head component as it is read ('ncols' and 'nrows', respectively).
    A ValueError is raised if either the number of data elements in a line
    is not equal to the expected number of columns, or if the number of
    data lines is not equal to the expected number of rows.
    """
    headmeta = []
    datarows = []
    tailmeta = []
    ncols = -1
    nrows = -1
    isheader = 1
    f = open(fname,'r')
    for line in f:
        line=line.strip()
        if not line: continue
        a = line.split()
        # Assume header lines contain only 2 whitespace-split elements
        if (len(a)==2) and isheader:
            headmeta.append(line)
            if   a[0]=='ncols': ncols = float(a[1])
            elif a[0]=='nrows': nrows = float(a[1])
        # Else if the number of elements==ncols the line is probably data
        elif len(a)==ncols:
            isheader = 0
            datarows.append(a)
        # Else if we don't have the expected number of data rows there is
        # an error
        elif len(datarows)!=nrows:
            raise ValueError("Line contains data but not of length "+str(ncols))
        # Anything else is tail metadata
        else:
            tailmeta.append(line)
    f.close()
    return headmeta,datarows,tailmeta

def arcasciihdr_to_dict(header):
    """
    Read an Arc ASCII header list of strings (from split_asciigrid) into an
    OrderedDict. All keys are returned as lowercase. All values are returned as
    strings or a list of strings.
    """
    meta = OrderedDict()
    r = re.compile("^(\S+)\s+(.+)$")
    for line in header:
        m = re.match(r, line.strip())
        if m:
            meta[m.group(1).lower()] = m.group(2)
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
        'xstep':float(meta['cellsize']),
        'xnum':int(meta['ncols']),
        'xunits':'degrees_east',
        'ymin':min(latvec),
        'ymax':max(latvec),
        'ystep':float(meta['cellsize']),
        'ynum':int(meta['nrows']),
        'yunits':'degrees_north'
    })

    # Return a tuple
    return latvec,lonvec,datadict

def asciigrid_to_numpy(datarows, meta, datadict=None):
    """
    Convert a list of lists of (string) values to a 3D NumPy array.
    Arg meta is not used but is retained for consistency.
    """
    # Copy datarows into a 3D array
    # Want shape to align with time,latitude,longitude
    data = np.array(datarows,dtype=np.float32,ndmin=3)

    # Create/update a datadict, which includes standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({'datatype':'f4'})  # Float32

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
    # Check missing value
    # Some files gave an incorrect missing value, so we check for the
    # actual missing value before using/replacing it.
    checkfor = 99999.9
    if np.sum(np.where(input_data==checkfor,1,0)) > input_data.size*0.7:
        input_dict.update({'missing':checkfor})

    # Copy dict
    if input_dict is None:
        output_dict = dict()
    else:
        output_dict = copy.deepcopy(input_dict)

    # Output dimensions
    (xs, xn, xc) = (112, 841, 0.05)  # start, number, cellsize
    (ys, yn, yc) = (-44, 681, 0.05)  # start, number, cellsize
    output_lon = nr.create_vector(xs, xn, xc)
    output_lat = nr.create_vector(ys, yn, yc)
    output_lat = output_lat[::-1]    # reverse elements

    # Create output array
    output_data = np.zeros((output_lat.size,output_lon.size)) \
                                      + input_dict['missing']

    # Copy data onto output grid
    output_data = nr.copy_grids(input_data[0,:,:],input_lon,input_lat,
                                output_data,output_lon,output_lat)
    output_data.shape = (1,yn,xn)
    output_dict.update({'xmin':min(output_lon),
                        'xmax':max(output_lon),
                        'xstep':xc,
                        'xnum':xn,
                        'ymin':min(output_lat),
                        'ymax':max(output_lat),
                        'ystep':yc,
                        'ynum':yn,})

    # Reduce precision of values to 1 decimal place and convert to f32
    output_data = output_data.round(decimals=1)
    output_data = np.float32(output_data)
    output_dict.update({'datatype':'f4'})

    # Change missing value
    miss = -999.0
    if input_dict['missing'] != miss:
        output_data = nr.replace_values(output_data,input_dict['missing'],miss)
        output_dict.update({'missing':miss})

    return (output_data, output_lat, output_lon, output_dict)

def set_attributes(datadict):
    """
    See http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=raster_to_float_%28conversion%29
    """
    fltmeta = OrderedDict()
    fltmeta['ncols'] = "{0:d}".format(datadict['xnum'])
    fltmeta['nrows'] = "{0:d}".format(datadict['ynum'])
    fltmeta['xllcorner'] = "{0:f}".format(datadict['xmin']-datadict['xstep']/2)
    fltmeta['yllcorner'] = "{0:f}".format(datadict['ymin']-datadict['ystep']/2)
    fltmeta['cellsize'] = "{0:f}".format(datadict['ystep'])
    fltmeta['nodata_value'] = "{0:f}".format(datadict['missing'])
    if sys.byteorder == "little":
        fltmeta['byteorder'] = 'LSBFIRST'
    else:
        fltmeta['byteorder'] = 'MSBFIRST'
    return fltmeta
    
def asciigrid_to_flt(arcfilename,fileroot):
    """
    The main routine that calls the calls other routines to prepare the data
    and metadata and create the Arc Float file.
    """
    # Read ascii grid file
    asciihead,asciidata,asciitail = split_asciigrid(arcfilename)
    meta = arcasciihdr_to_dict(asciihead)

    # Create numpy components
    latvec,lonvec,datadict = set_latlon(meta)
    data,datadict = asciigrid_to_numpy(asciidata,meta,datadict)

    # Resample or edit array
    # Add a default no_data value if required.
    miss = -999.0
    if 'nodata_value' in meta: miss = float(meta['nodata_value'])
    datadict['missing'] = miss
    data,latvec,lonvec,datadict = resample_array(data,latvec,lonvec,datadict)

    # Prepare metadata
    attributes = set_attributes(datadict)

    # Write Arc Float files
    fltfile = fileroot+'.flt'
    hdrfile = fileroot+'.hdr'
    if os.path.exists(fltfile): os.remove(fltfile)
    if os.path.exists(hdrfile): os.remove(hdrfile)
    data.tofile(fltfile)
    f = open(hdrfile,'w')
    for k in attributes:
        f.write(k+"  "+attributes[k]+'\n')
    f.close
    print 'Wrote:',fltfile
    print 'Wrote:',hdrfile

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:"
        print "  ", sys.argv[0], "path/to/arc_ascii_grid_file"
        print "Notes:"
        print "  You may want to you add or amend some of the array manipulation"
        print "  routines given in the subroutine resample_array. This subroutine is not"
        print "  invoked by default."
        exit()
    else:
        fileroot = re.sub('.\w+$','',sys.argv[1])
        asciigrid_to_flt(sys.argv[1],fileroot)

