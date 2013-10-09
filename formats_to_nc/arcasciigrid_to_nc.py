# ARCASCIIGRID_TO_NC.PY
#
# Author:
# Matt Paget, CSIRO Marine and Atmospheric Research, Canberra, Australia.
#
# The software is open source and is released under Creative Commons
# Attribution (CC-BY).  Comments, questions, revisions and recommendations
# can be directed to <data at auscover.org.au>.
#
# Version:
# 1.0   12 April 2013
# 1.1   6 May 2013
#       Changed '+ncols' to '+str(ncols)' in split_asciigrid raise statement.
#       Added more helpful usage text.
#
# This code was adapted from code written to convert Bureau of Meteorology
# data in Arc ascii grid format to netcdf by CSIRO.
#
# Basically the ascii text is parsed into a list of strings for the header and
# for the data. The header defines the coordinate variables and the missing
# data value. The data is converted to a numpy array. A netcdf file is built
# from the coordinate variables and the data array (using netcdf_builder.py).
#
# The time dimension is generally parsed from the filename as its not readily
# provided in the header information and is rarely added as a parsable element
# in the file. As such the parsing regular expression and time period value
# may need to be adjusted for your application.
#
# Extra hooks are provided to capture per-file metadata in JSON format files.
# The JSON metadata can be edited, collated and added to the netcdf file as
# attributes as a seperate process (see json_handler.py).

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
            if   a[0]=='ncols': ncols=float(a[1])
            elif a[0]=='nrows': nrows=float(a[1])
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


def asciigrid_to_numpy(datarows,metalist):
    """
    Convert a list of lists of (string) values to a 2D NumPy array.  Create
    dimension vectors based on extent and cell-size information in 'metalist'.
    """
    # Parse metalist into a dict
    # python 2.6 safe
    meta = dict( (y[0].lower(),float(y[1]))
                 for y in [x.split() for x in metalist])
    # python 2.7 preferred
    # meta = {y[0].lower():float(y[1]) for y in [x.split() for x in metalist]}
    miss = -999.0
    if 'nodata_value' in meta: miss = meta['nodata_value']

    # Convert corners to centers if required
    if 'xllcorner' in meta:
        meta['xllcenter'] = meta['xllcorner'] + meta['cellsize']/2.0
    if 'yllcorner' in meta:
        meta['yllcenter'] = meta['yllcorner'] + meta['cellsize']/2.0

    # ASCII dimensions
    input_lon = nr.create_vector(meta['xllcenter'], \
                                 meta['ncols'], \
                                 meta['cellsize'])
    input_lat = nr.create_vector(meta['yllcenter'], \
                                 meta['nrows'], \
                                 meta['cellsize'])
    input_lat = input_lat[::-1]    # reverse elements

    # Copy datarows into a 2D array
    input_data = np.array(datarows,dtype=np.float64)

    # Create a dict for the output dimensions
    input_dict = {'xmin':min(input_lon),
                  'xmax':max(input_lon),
                  'xstep':meta['cellsize'],
                  'xnum':meta['ncols'],
                  'xunits':'degrees_east',
                  'ymin':min(input_lat),
                  'ymax':max(input_lat),
                  'ystep':meta['cellsize'],
                  'ynum':meta['nrows'],
                  'yunits':'degrees_north',
                  'missing':miss}

    return input_data,input_lat,input_lon,input_dict


def resample_array(input_data, input_lat, input_lon, input_dict):
    """
    Code snippets for resampling the input array on to a new grid, changing
    the missing value or reducing the data type or floating point precision
    of the array.

    Edit the components of this routine for your particular needs. Just check
    that input_data and output_data arrays are being used and referenced as
    required.
    """
    # Copy dict
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
    output_data = nr.copy_grids(input_data,input_lon,input_lat,
                             output_data,output_lon,output_lat)

    # Update dict
    output_dict.update({'xmin':min(output_lon),
                        'xmax':max(output_lon),
                        'xstep':xc,
                        'xnum':xn,
                        'ymin':min(output_lat),
                        'ymax':max(output_lat),
                        'ystep':yc,
                        'ynum':yn,})

    # Reduce precision of values to 1 decimal place and convert to f32
    output_data = input_data.round(decimals=1)
    output_data = np.float32(input_data)

    # Change missing value
    miss = -999
    output_data = nr.replace_values(input_data,input_dict['missing'],miss)
    output_dict.update({'missing':miss})

    return (output_data, output_lat, output_lon, output_dict)


def get_datetime(fname,datadict):
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

    d1 = datetime(d1year,d1month,d1day)
    d2 = datetime(d2year,d2month,d1day) - timedelta(days=1)

    # Metadata elements
    # See http://en.wikipedia.org/wiki/ISO_8601
    datadict.update({'dcreate':"unknown",
                     'dmodify':datetime.utcnow().strftime("%Y%m%dT%H%M%S"),
                     'tmin':d1.strftime("%Y-%m-%d"),
                     'tmax':d2.strftime("%Y-%m-%d"),
                     'tduration':"P1M",
                     'tresolution':"P1M"})

    # d1 will be used by make_netcdf as a list, so its possible d1 should be
    # returned as a list here
    return d1,d2,datadict


def get_varname(fname,datadict):
    """
    Code example for creating the variable name. Generally, the variable name
    will be in the filename as its not included in the standard header
    information.

    Edit the components of this routine for your particular needs.
    """
    varname = "varname1"
    varlong = "long variable name"
    varunit = "unit"
    datadict.update({'varname':varname,
                     'varlong':varlong,
                     'varunits':varunit})
    return datadict


def get_attributes(datadict):
    """
    """
    # Define a new metadata dict to control the order of elements
    meta = OrderedDict()

    # Create some history text. Added in a list so the join character can be
    # changed easily.
    history = [datadict['dmodify']+": Reformatted to NetCDF."]
    # Add some details about any modifications to the data
    # history.extend("Reduced precision of values to 1 decimal place.")
    # history.extend("Created a no-data value of {0}.".format(datadict['missing']))
    history.extend(["If present, the \"sha1_arcasciigrid\" attribute of a variable is the SHA1 hex digest of the input Arc ASCII Grid. This allows the Arc ASCII Grid and netCDF files to be uniquely linked irrespective of filename changes."])

    # The date*, geospatial* and time* attributes come from the Attribute
    # Convention for Dataset Discovery (ACDD). See,
    # http://www.unidata.ucar.edu/software/netcdf/conventions.html
    # These are optional for a NetCDF file but no harm in having them.
    vname = datadict['varname']
    meta['history'] = ' '.join(history)
    meta['date_created'] = datadict['dcreate']
    meta['date_modified'] = datadict['dmodify']
    meta['geospatial_lat_min'] = "{0:.2f}".format(datadict['ymin'])
    meta['geospatial_lat_max'] = "{0:.2f}".format(datadict['ymax'])
    meta['geospatial_lat_step'] = "{0:.2f}".format(datadict['ystep'])
    meta['geospatial_lat_units'] = datadict['yunits']
    meta['geospatial_lon_min'] = "{0:.2f}".format(datadict['xmin'])
    meta['geospatial_lon_max'] = "{0:.2f}".format(datadict['xmax'])
    meta['geospatial_lon_step'] = "{0:.2f}".format(datadict['xstep'])
    meta['geospatial_lon_units'] = datadict['xunits']
    meta['time_coverage_start'] = datadict['tmin']
    meta['time_coverage_end'] = datadict['tmax']
    meta['time_coverage_duration'] = datadict['tduration']
    meta['time_coverage_resolution'] = datadict['tresolution']
    meta[vname+':long_name'] = datadict['varlong']
    meta[vname+':units'] = datadict['varunits']
    meta[vname+':_FillValue'] = datadict['missing']
    meta[vname+':grid_mapping'] = 'crs'
    meta[vname+':sha1_arcasciigrid'] = datadict['sha1_ascii']
    meta['latitude:units'] = datadict['yunits']
    meta['longitude:units'] = datadict['xunits']
    meta['crs:grid_mapping_name'] = "latitude_longitude"
    meta['crs:longitude_of_prime_meridian'] = 0.0
    meta['crs:semi_major_axis'] = 6378137.0
    meta['crs:inverse_flattening'] = 298.257223563
    return meta


def asciigrid_to_nc(arcfilename,fileroot):
    """
    The main routine that calls the calls other routines to prepare the data
    and metadata and create the netcdf file.
    """
    # Read ascii grid file and prepare data
    asciihead,asciidata,asciitail = split_asciigrid(arcfilename)
    data,latvec,lonvec,datadict = asciigrid_to_numpy(asciidata,asciihead)
    #data,latvec,lonvec,datadict = resample_array(data,latvec,lonvec,datadict)

    # Prepare metadata
    datadict['sha1_ascii'] = hashlib.sha1(
                               open(arcfilename,'r').read()).hexdigest()
    d1,d2,datadict = get_datetime(arcfilename,datadict)
    timevec = [d1]
    datadict = get_varname(arcfilename,datadict)
    varname = datadict['varname']
    vartype = 'f4'
    attributes = get_attributes(datadict)

    # Netcdf options
    # http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    nc_format = 'NETCDF4_CLASSIC'
    nc_compress = True
    debug = False

    # Write netcdf file
    if os.path.exists(fileroot+'.nc'): os.remove(fileroot+'.nc')
    ncobj = nb.nc_open(fileroot+'.nc','w',format=nc_format)
    nb.nc_set_timelatlon(ncobj,None,len(latvec),len(lonvec)) # unlimited time
    nb.nc_set_var(ncobj,varname,dtype=vartype,
                  fill=datadict['missing'],zlib=nc_compress)
    nb.nc_set_var(ncobj,'crs',dims=(),dtype="i4")  # Grid mapping container
    nb.nc_add_time(ncobj,timevec)
    nb.nc_add_data(ncobj,'latitude',latvec)
    nb.nc_add_data(ncobj,'longitude',lonvec)
    if debug:
        print varname,data.shape
        nb.nc_show_dims(ncobj)
    # nc_add_data should work but is presently broken. Use direct method
    #nb.nc_add_data(ncobj,varname,data)
    ncobj.variables[varname][0,:,:] = data
    nb.nc_set_attributes(ncobj,attributes)
    nb.nc_close(ncobj)
    print 'Wrote:',fileroot+'.nc'

    # Write metadata to json
    if os.path.exists(fileroot+'.json'): os.remove(fileroot+'.json')
    jh.json_dump(attributes,fileroot+'.json')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:"
        print "  ", sys.argv[0], "path/to/arc_ascii_grid_file"
        print "Notes:"
        print "  At least two subroutines need to be edited for each type of data, namely:"
        print "    get_datetime - Create a time coordinate value, possibly as a regular"
        print "                   expression match on the input filename or hardcoded."
        print "    get_varname  - Create the netcdf variable name, possibly as a regular"
        print "                   expression match on the input filename or harcoded."
        print "  Additionally, you may want to you and amend some of the array manipulation"
        print "  routines given in the subroutine resample_array. This subroutine is not"
        print "  invoked by default."
        exit()
    else:
        fileroot = re.sub('.\w+$','',sys.argv[1])
        asciigrid_to_nc(sys.argv[1],fileroot)
