# ARCFLTGRID_TO_NC.PY
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

def floatgrid_to_numpy(fltfile, meta, datadict=None):
    """
    Convert a binary file to a 3D NumPy array.
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
    data.shape = (1, int(meta['nrows']), int(meta['ncols']))

    # Create/update a datadict, which includes standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({'datatype':'f4'})  # Float32. Endian-ness is managed by
                                        # the netCDF library.

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
    # See arcasciigrid_to_nc.py for code snippet examples
    pass

def set_datetime(fname, datadict=None):
    """
    Code example for creating the date/time information. Generally, the
    date/time information will be in the filename as its not included in the
    standard header information.

    Edit the components of this routine for your particular needs.
    """
    # Example: ..../arcfltgrid_WRel2_20130113.{flt,hdr}
    m = re.search('(\d{4})(\d\d)(\d\d)',fname)
    d1year = int(m.group(1))
    d1month = int(m.group(2))
    d1day = int(m.group(3))

    # Monthly = (d1month + 1) - (1 day)
    d2year = d1year
    d2month = d1month + 1
    if d2month == 13:
        d2month = 1
        d2year = d2year + 1

    d1 = dt.datetime(d1year,d1month,d1day)
    d2 = dt.datetime(d2year,d2month,d1day) - dt.timedelta(days=1)

    tdur = "P1M"

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
    part = re.search('\w+',fname)
    varname = "varname1"
    varlong = "long variable name"
    varunit = "unit"

    # Create/update a datadict, which adds some standardised labels for later
    if datadict is None: datadict = dict()
    datadict.update({'varname':varname,
                     'varlong':varlong,
                     'varunits':varunit})

    # Capture the SHA1 digest of the input flt file
    datadict['sha1'] = hashlib.sha1(open(fname,'r').read()).hexdigest()

    return datadict

def set_attributes(fltfile, meta, datadict=None):
    """
    Generate a dictionary with keys representing the actual attribute names
    to be added to the netCDF file. Most of the values come from the datadict
    dictionary, which has sanitised and collated most of the required and
    available information. Some additional keys/values are added here, such as
    the history attribute.
    The returned dictionary can be added directly to a netCDF file (via
    netcdf_builder.set_attributes()).
    """
    # Define a new metadata dict to control the order of elements
    ncmeta = OrderedDict()

    # Create some history text. Added in a list so the join character can be
    # changed easily.
    history = [datadict['dmodify']+": Reformatted to NetCDF."]
    history.extend(["Input file: "+fltfile])
    # Add some details about any modifications to the data
    # history.extend("Reduced precision of values to 1 decimal place.")
    # history.extend("Created a no-data value of {0}.".format(datadict['missing']))
    history.extend(["If present, the \"sha1_arcfltgrid\" attribute of a variable is the SHA1 hex digest of the input Arc Float Grid. This allows the Arc Float Grid and netCDF files to be uniquely linked irrespective of filename changes."])

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
    ncmeta['geospatial_lat_units'] = datadict['yunits']
    ncmeta['geospatial_lat_resolution'] = "{0:.2f}".format(datadict['ystep'])
    ncmeta['geospatial_lon_min'] = "{0:.2f}".format(datadict['xmin'])
    ncmeta['geospatial_lon_max'] = "{0:.2f}".format(datadict['xmax'])
    ncmeta['geospatial_lon_units'] = datadict['xunits']
    ncmeta['geospatial_lon_resolution'] = "{0:.2f}".format(datadict['xstep'])
    ncmeta['time_coverage_start'] = datadict['tmin']
    ncmeta['time_coverage_end'] = datadict['tmax']
    ncmeta['time_coverage_duration'] = datadict['tduration']
    ncmeta['time_coverage_resolution'] = datadict['tresolution']
    ncmeta[vname+':long_name'] = datadict['varlong']
    ncmeta[vname+':units'] = datadict['varunits']
    ncmeta[vname+':grid_mapping'] = 'crs'
    ncmeta[vname+':sha1_arcfloatgrid'] = datadict['sha1']
    ncmeta['latitude:units'] = datadict['yunits']
    ncmeta['longitude:units'] = datadict['xunits']
    ncmeta['crs:grid_mapping_name'] = "latitude_longitude"
    ncmeta['crs:longitude_of_prime_meridian'] = 0.0
    ncmeta['crs:semi_major_axis'] = 6378137.0
    ncmeta['crs:inverse_flattening'] = 298.257223563
    return ncmeta

def floatgrid_to_nc(arcfilename,fileroot):
    """
    The main routine that calls the calls other routines to prepare the data
    and metadata and create the netCDF file.
    """
    # Set input filenames
    fltfile = re.sub('\.\w+$','.flt',arcfilename)
    hdrfile = re.sub('\.\w+$','.hdr',arcfilename)
    
    # Read metadata data
    meta = arcflthdr_to_dict(hdrfile)

    # Create numpy components
    latvec,lonvec,datadict = set_latlon(meta)
    data,datadict = floatgrid_to_numpy(fltfile,meta,datadict)

    # Resample or edit array
    # Add a default no_data value if required.
    miss = -999.0
    if 'nodata_value' in meta: miss = meta['nodata_value']
    datadict['missing'] = miss
    #data,latvec,lonvec,datadict = resample_array(data,latvec,lonvec,datadict)

    # Prepare time, variable name and metadata
    d1,d2,datadict = set_datetime(fltfile,datadict)
    datadict = set_varname(fltfile,datadict)
    attributes = set_attributes(fltfile,meta,datadict)

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
    ncobj = nb.ncopen(fileroot+'.nc','w',format=nc_format)
    nb.set_timelatlon(ncobj,None,len(latvec),len(lonvec)) # unlimited time
    nb.set_variable(ncobj,varname,dtype=vartype,fill=fillval,zlib=nc_compress)
    nb.set_variable(ncobj,'crs',dims=(),dtype="i4")  # Grid mapping container
    nb.add_time(ncobj,timevec)
    nb.add_data(ncobj,'latitude',latvec)
    nb.add_data(ncobj,'longitude',lonvec)
    if debug:
        print varname,data.shape
        nb.show_dimensions(ncobj)
    # nb.add_data should work but is presently broken. Use direct method
    #nb.add_data(ncobj,varname,data)
    #ncobj.variables[varname][0,:,:] = data  # 2D numpy array
    ncobj.variables[varname][:] = data  # 3D numpy array
    nb.set_attributes(ncobj,attributes)
    nb.ncclose(ncobj)
    print 'Wrote:',fileroot+'.nc'

    # Write metadata to json
    if os.path.exists(fileroot+'.json'): os.remove(fileroot+'.json')
    jh.json_dump(attributes,fileroot+'.json')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:"
        print "  ", sys.argv[0], "path/to/arc_float_grid_file"
        print "Notes:"
        print "  At least two subroutines need to be edited for each type of data, namely:"
        print "    set_datetime - Create a time coordinate value, possibly as a regular"
        print "                   expression match on the input filename or hardcoded."
        print "    set_varname  - Create the netcdf variable name, possibly as a regular"
        print "                   expression match on the input filename or harcoded."
        print "  Additionally, you may want to you add or amend some of the array manipulation"
        print "  routines given in the subroutine resample_array. This subroutine is not"
        print "  invoked by default."
        exit()
    else:
        fileroot = re.sub('.\w+$','',sys.argv[1])
        floatgrid_to_nc(sys.argv[1],fileroot)

