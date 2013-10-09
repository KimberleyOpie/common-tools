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

def asciigrid_to_numpy(datarows, meta, datadict=None):
    """
    Convert a list of lists of (string) values to a 3D NumPy array.
    Args meta and datadict are not used but are retained for consistency.
    """
    # Copy datarows into a 3D array
    # Want shape to align with time,latitude,longitude
    data = np.array(datarows,dtype=np.float64,ndmin=3)

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
    # Copy dict
    output_dict = copy.deepcopy(input_dict)

    # Create/update a datadict, which adds some standardised labels for later
    if datadict is None: datadict = dict()

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
    output_dict.update({'datatype':'f4'})

    # Change missing value
    miss = -999
    output_data = nr.replace_values(input_data,input_dict['missing'],miss)
    output_dict.update({'missing':miss})

    return (output_data, output_lat, output_lon, output_dict)

def set_datetime(fname, datadict=None):
    """
    Code example for creating the date/time information. Generally, the
    date/time information will be in the filename as its not included in the
    standard header information.

    Edit the components of this routine for your particular needs.
    """
    print fname
    m = re.search('(\d{4})(\d\d)(\d\d)(\d{4})(\d\d)(\d\d)',fname)
    d1year = int(m.group(1))
    d1month = int(m.group(2))
    d1day = int(m.group(3))
    d2year = int(m.group(4))
    d2month = int(m.group(5))
    d2day = int(m.group(6))

    d1 = dt.datetime(d1year,d1month,d1day)
    d2 = dt.datetime(d2year,d2month,d2day)

    dur = "P1M"
    if d1day == d2day: dur = "P1D"

    # Create/update a datadict, which adds some standardised labels for later
    if datadict is None: datadict = dict()

    # Metadata elements
    # See http://en.wikipedia.org/wiki/ISO_8601
    datadict.update({
        'dcreate':"unknown",
        'dmodify':dt.datetime.utcnow().strftime("%Y%m%dT%H%M%S"),
        'tmin':d1.strftime("%Y-%m-%d"),
        'tmax':d2.strftime("%Y-%m-%d"),
        'tduration':dur,
        'tresolution':dur
    })

    # Return as a tuple
    return d1,d2,datadict

def set_varname(fname, datadict=None):
    """
    Code example for creating the variable name. Generally, the variable name
    will be in the filename as its not included in the standard header
    information.

    Edit the components of this routine for your particular needs.
    """
    part = re.search('\w+',fname)
    varname = "varname1"
    varlong = "long variable name"
    varunit = "unit"

    # Create/update a datadict, which adds some standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({'varname':varname,
                     'varlong':varlong,
                     'varunits':varunit})

    # Capture the SHA1 digest of the input file
    datadict['sha1'] = hashlib.sha1(open(arcfilename,'r').read()).hexdigest()

    return datadict

def set_attributes(fname, meta, datadict):
    """
    """
    # Define a new metadata dict to control the order of elements
    ncmeta = OrderedDict()

    # Create some history text. Added in a list so the join character can be
    # changed easily.
    history = [datadict['dmodify']+": Reformatted to NetCDF."]
    history.extend(["Input file: "+fname])
    # Add some details about any modifications to the data
    # history.extend("Reduced precision of values to 1 decimal place.")
    # history.extend("Created a no-data value of {0}.".format(datadict['missing']))
    history.extend(["If present, the \"sha1_arcasciigrid\" attribute of a variable is the SHA1 hex digest of the input Arc ASCII Grid. This allows the Arc ASCII Grid and netCDF files to be uniquely linked irrespective of filename changes."])

    # The date*, geospatial* and time* attributes come from the Attribute
    # Convention for Dataset Discovery (ACDD). See,
    # http://www.unidata.ucar.edu/software/netcdf/conventions.html
    # These are optional for a NetCDF file but no harm in having them.
    vname = datadict['varname']
    ncmeta['history'] = ' '.join(history)
    ncmeta['date_created'] = datadict['dcreate']
    ncmeta['date_modified'] = datadict['dmodify']
    ncmeta['geospatial_lat_min'] = "{0:.2f}".format(datadict['ymin'])
    ncmeta['geospatial_lat_max'] = "{0:.2f}".format(datadict['ymax'])
    ncmeta['geospatial_lat_step'] = "{0:.2f}".format(datadict['ystep'])
    ncmeta['geospatial_lat_units'] = datadict['yunits']
    ncmeta['geospatial_lon_min'] = "{0:.2f}".format(datadict['xmin'])
    ncmeta['geospatial_lon_max'] = "{0:.2f}".format(datadict['xmax'])
    ncmeta['geospatial_lon_step'] = "{0:.2f}".format(datadict['xstep'])
    ncmeta['geospatial_lon_units'] = datadict['xunits']
    ncmeta['time_coverage_start'] = datadict['tmin']
    ncmeta['time_coverage_end'] = datadict['tmax']
    ncmeta['time_coverage_duration'] = datadict['tduration']
    ncmeta['time_coverage_resolution'] = datadict['tresolution']
    ncmeta[vname+':long_name'] = datadict['varlong']
    ncmeta[vname+':units'] = datadict['varunits']
    ncmeta[vname+':_FillValue'] = datadict['missing']
    ncmeta[vname+':grid_mapping'] = 'crs'
    ncmeta[vname+':sha1_arcasciigrid'] = datadict['sha1']
    ncmeta['latitude:units'] = datadict['yunits']
    ncmeta['longitude:units'] = datadict['xunits']
    ncmeta['crs:grid_mapping_name'] = "latitude_longitude"
    ncmeta['crs:longitude_of_prime_meridian'] = 0.0
    ncmeta['crs:semi_major_axis'] = 6378137.0
    ncmeta['crs:inverse_flattening'] = 298.257223563
    return ncmeta

def asciigrid_to_nc(arcfilename,fileroot):
    """
    The main routine that calls the calls other routines to prepare the data
    and metadata and create the netcdf file.
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
    if 'nodata_value' in meta: miss = meta['nodata_value']
    datadict['missing'] = miss
    #data,latvec,lonvec,datadict = resample_array(data,latvec,lonvec,datadict)

    # Prepare time, variable name and metadata
    d1,d2,datadict = set_datetime(arcfilename,datadict)
    datadict = set_varname(arcfilename,datadict)
    attributes = set_attributes(arcfilename,meta,datadict)

    # Netcdf options
    # http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    nc_format = 'NETCDF4_CLASSIC'
    nc_compress = True
    debug = False

    # Write netcdf file
    if os.path.exists(fileroot+'.nc'): os.remove(fileroot+'.nc')
    varname = datadict['varname']
    vartype = datadict['datatype']
    fillval = datadict['missing']
    timevec = [d1]
    ncobj = nb.nc_open(fileroot+'.nc','w',format=nc_format)
    nb.nc_set_timelatlon(ncobj,None,len(latvec),len(lonvec)) # unlimited time
    nb.nc_set_var(ncobj,varname,dtype=vartype,fill=fillval,zlib=nc_compress)
    nb.nc_set_var(ncobj,'crs',dims=(),dtype="i4")  # Grid mapping container
    nb.nc_add_time(ncobj,timevec)
    nb.nc_add_data(ncobj,'latitude',latvec)
    nb.nc_add_data(ncobj,'longitude',lonvec)
    if debug:
        print varname,data.shape
        nb.nc_show_dims(ncobj)
    # nc_add_data should work but is presently broken. Use direct method
    #nb.nc_add_data(ncobj,varname,data)
    #ncobj.variables[varname][0,:,:] = data  # 2D numpy array
    ncobj.variables[varname][:] = data  # 3D numpy array
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
        print "  Additionally, you may want to you add or amend some of the array manipulation"
        print "  routines given in the subroutine resample_array. This subroutine is not"
        print "  invoked by default."
        exit()
    else:
        fileroot = re.sub('.\w+$','',sys.argv[1])
        asciigrid_to_nc(sys.argv[1],fileroot)

