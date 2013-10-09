# BOM_SOLAR_ASCII_TO_NETCDF.PY
#
# The input proprietary format supported by this example is a Bureau of
# Meteorology (BoM) ArcView ASCII grid file.  The format contains a header
# (defines the shape and size of the grid), space-separated ASCII data, and
# tail metadata that contains internal BoM processing information.
#
# The ASCII data is converted to a NumPy array.  The precision of values within
# the grid are reduced to 1 decimal place (for consistency between variables),
# and the missing value is redefined to -999.0.
#
# File-specific metadata (from the header, the tail, and from the processing in
# this script) is written to a JSON file.  In a second-step, the JSON file is
# added to the NetCDF file by netcdf_json_wrapper.py.  This two-step
# process allows manual editing or checking of the metadata if required.
#
# By default the following files are created:
# - [verbose filename].json
# - [verbose filename].nc
# The selection of outputs can be customised by parameters to
# bom_ascii_to_flt_nc().



import os
import numpy as np
from scipy import interpolate
import re
import hashlib
from datetime import timedelta
from datetime import datetime
from collections import OrderedDict
import netcdf_builder as nb
import json_handler as jh
import numpy_routines as nr
from solar_obstime_lookup import get_lookup_list


def split_bom_file(fname):
    """
    Split a BoM ASCII grid file into its header, data and tail components
    End-of-line characters and whitespace at either end of a line are
    striped off.  The head and tail components are returned as a list of
    lines.  The data component is returned as a list of lists, with each
    inner list being a list of (string) values in a line.

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
            if   a[0].lower()=='ncols': ncols=float(a[1])
            elif a[0].lower()=='nrows': nrows=float(a[1])
        # Else if the number of elements==ncols the line is probably data
        elif len(a)==ncols:
            isheader = 0
            datarows.append(a)
        # Else if we don't have the expected number of data rows there is
        # an error
        elif len(datarows)!=nrows:
            raise ValueError("Line contains data but not of length "+ncols)
        # Anything else is tail metadata
        else:
            tailmeta.append(line)
    f.close()
    return (headmeta,datarows,tailmeta)

    
def resample_data(datarows,metalist):
    """
    Convert a list of lists of (string) values to a 2D NumPy array.  Create
    dimension vectors based on extent and cell-size information in 'metalist'.
    Create target dimension vectors based on pre-defined extents.
    Check the correctness of the given missing value (in metalist) against
    the data.  Replace all missing values with a standard missing value.

    Copy data from the input array to an output array defined by the
    target dimension vectors.  The precision of values in the output array are
    reduced to 1 decimal place (for consistency).  Return the output array,
    target dimension vectors and a dict that gives the extents of the target
    dimension vectors and the new missing value.
    """
    # Parse metalist into a dict
    meta = {y[0].lower():float(y[1]) for y in [x.split() for x in metalist]}
    miss = -999.0

    # ASCII dimensions
    input_lon = nr.create_vector(meta['xllcorner']+meta['cellsize']/2.0, \
    				 meta['ncols'], \
    				 meta['cellsize'])
    input_lat = nr.create_vector(meta['yllcorner']+meta['cellsize']/2.0, \
    				 meta['nrows'], \
    				 meta['cellsize'])
    input_lat = input_lat[::-1]    # reverse elements
    # Output dimensions - not remapping, so commented out
    #(xs, xn, xc) = (112, 841, 0.05)  # start, number, cellsize
    #(ys, yn, yc) = (-44, 681, 0.05)  # start, number, cellsize
    #output_lon = nr.create_vector(xs, xn, xc)
    #output_lat = nr.create_vector(ys, yn, yc)
    #output_lat = output_lat[::-1]    # reverse elements

    # Copy datarows into a 2D array
    input_data = np.array(datarows,dtype=np.float64)
    meta['nodata_value'] = \
        check_bom_missing(input_data[0,:], 99999.9, meta['nodata_value'])
    if meta['nodata_value']!=miss:
        input_data = nr.replace_values(input_data,meta['nodata_value'],miss)
        print "Replaced missing data %s with %s" % (meta['nodata_value'],miss)
    
    # Create output array
    #output_data = np.zeros((output_lat.size,output_lon.size))+miss
    
    # Copy data onto output grid
    #output_data = nr.copy_grids(input_data,input_lon,input_lat,
    #                         output_data,output_lon,output_lat)

    # Reduce precision of values to 1 decimal place and convert to f32
    #output_data = output_data.round(decimals=1)
    #output_data = np.float32(output_data)
    input_data = input_data.round(decimals=1)
    input_data = np.float32(input_data)

    input_dict = {'xmin':min(input_lon).round(decimals=3),
                  'xmax':max(input_lon).round(decimals=3),
                  'xstep':meta['cellsize'],
                  'xnum':meta['ncols'],
                  'xunits':'degrees_east',
                  'ymin':min(input_lat).round(decimals=3),
                  'ymax':max(input_lat).round(decimals=3),
                  'ystep':meta['cellsize'],
                  'ynum':meta['nrows'],
                  'yunits':'degrees_north',
                  'missing':miss}

    return (input_data, input_lat, input_lon, input_dict)


def check_bom_missing(arr,test,default):
    """
    Test for the possibility of a different missing value in 'arr'
    rather than assuming that the given missing value is correct.

    If occurence of 'test' is > 70% of size of 'arr' then return 'test'.
    >>>check_missing(numpy.zeros(5)+99999.9, 99999.9, -99.99)
    99999.9
    >>>

    If occurence of 'test' is <= 70% of size of 'arr' then return 'default'.
    >>>check_missing(numpy.zeros(5), 99999.9, -99.99)
    -99.99
    >>>
    """
    if np.sum(np.where(arr==test,1,0)) > arr.size*0.7:
        return test
    else:
        return default


def create_obs_time(latvec,lonvec,d1):
    """
    Create a data layer of the observation time for each pixel from a lookup
    table. The lookup table defines the minutes after the start hour (d1) at
    five degree increments of latitude from -10 to -44 N.
    Linearly interpolate from the five degree increment values to the latitude
    vector.
    A 2D data layer the same size as latvec,lonvec is returned.
    """
    # Choose lookup table based on date and hour
    xlook,ylook = None,None
    for obj in get_lookup_list():
        if obj.is_applicable(d1):
            xlook,ylook = obj.get_lookup()
            break
    if xlook is None:
        raise ValueError("Date-Hour not found in solar_observation_time_lookup: "+d1.strftime("%Y-Ym-%dT%HZ"))
        return None
    # The lookup table values are nearly linear but not quite (R^2>=0.9944)
    # So going to choose cubic spline interpolation
    # http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
    # Want monotonicaly increasing xlook (latitudes).
    xlook = xlook[::-1]
    ylook = ylook[::-1]
    tck = interpolate.splrep(xlook,ylook,s=0)
    obsvec = interpolate.splev(latvec,tck,der=0)
    # Replicate obsvec along width of lonvec
    X,obsarr = np.meshgrid(lonvec,obsvec)
    # Reduce precision to that given by the lookup table and convert to float
    obsarr = obsarr.round(decimals=1)
    obsarr = np.float32(obsarr)
    return obsarr


def create_meta(d1,d2,datadict):
    """
    Create an OrderedDict of global and variable metadata specific to this
    file from various sources of information including the start and end dates,
    the data extents dict and the history information.
    """
    meta = OrderedDict()

    # Make creation and modification dates
    dmodify = datetime.utcnow().strftime("%Y%m%dT%H%M%S")

    # Create multiline history text
    history = """
{:s}: Reformatted to NetCDF
Solar data are written with 1 decimal precision although all the original data
values are integers.
Observation times were interpolated from the date-satellite-latitude lookup
table to the latitude vector with a cubic spline method.  The interpolated
vector was then replicated across the grid for each longitude.  Observation
times were written with 1 decimal precision for consistency with the lookup
table.
Solar and Observation time data have a consistent no-data value of {:.0f}.
""".format(dmodify,datadict['missing'])

    # The date*, geospatial* and time* attributes come from the Attribute
    # Convention for Dataset Discovery (ACDD). See,
    # http://www.unidata.ucar.edu/software/netcdf/conventions.html
    meta['history'] = history
    meta['date_created'] = "unknown"
    meta['date_modified'] = dmodify
    meta['geospatial_lat_min'] = "{:.2f}".format(datadict['ymin'])
    meta['geospatial_lat_max'] = "{:.2f}".format(datadict['ymax'])
    meta['geospatial_lat_step'] = "{:.2f}".format(datadict['ystep'])
    meta['geospatial_lat_units'] = datadict['yunits']
    meta['geospatial_lon_min'] = "{:.2f}".format(datadict['xmin'])
    meta['geospatial_lon_max'] = "{:.2f}".format(datadict['xmax'])
    meta['geospatial_lon_step'] = "{:.2f}".format(datadict['xstep'])
    meta['geospatial_lon_units'] = datadict['xunits']
    meta['time_coverage_start'] = d1.strftime("%Y-%m-%dT%HZ")
    meta['time_coverage_end'] = d2.strftime("%Y-%m-%dT%HZ")
    meta['time_coverage_duration'] = "P1M"
    meta['time_coverage_resolution'] = "PT1H"
    meta['solar_dni:long_name'] = "Gridded hourly solar direct normal irradiance"
    meta['solar_dni:standard_name'] = "surface_downwelling_shortwave_flux_in_air_due_to_direct_fraction"
    meta['solar_dni:units'] = "W m-2"
    meta['solar_dni:grid_mapping'] = 'crs'
    meta['solar_ghi:long_name'] = "Gridded hourly solar global horizontal irradiance"
    meta['solar_ghi:standard_name'] = "surface_downwelling_shortwave_flux_in_air"
    meta['solar_ghi:units'] = "W m-2"
    meta['solar_ghi:grid_mapping'] = 'crs'
    meta['obs_time:long_name'] = 'Interpolated instantaneous pixel observation time relative to time dimension value'
    meta['obs_time:standard_name'] = 'time'
    meta['obs_time:units'] = 'minutes'
    meta['obs_time:grid_mapping'] = 'crs'
    meta['crs:grid_mapping_name'] = 'latitude_longitude'
    meta['crs:long_name'] = 'WGS 1984 datum'
    meta['crs:longitude_of_prime_meridian'] = str(0.0)
    meta['crs:semi_major_axis'] = str(6378137.0)
    meta['crs:inverse_flattening'] = str(298.257223563)
    return meta


def bom_ascii_to_nc(year,dates,froot):

    latvec = None
    lonvec = None
    ncobj = None
    latslice = slice(None,None,None)
    lonslice = slice(None,None,None)
    missed_dates = []

    for dti,dt in enumerate(dates):
        dni,ghi = get_solar_files(str(year),dt)
        index = (slice(dti,None,None),latslice,lonslice)

        if dni is not None:
            # Split the input file into metadata and data components
            dni_head,dni_rows,dni_history = split_bom_file(dni)
            # Resample the data
            dni_data,dni_lat,dni_lon,dni_dict = resample_data(dni_rows,dni_head)

        if ghi is not None:
            # Split the input file into metadata and data components
            ghi_head,ghi_rows,ghi_history = split_bom_file(ghi)
            # Resample the data
            ghi_data,ghi_lat,ghi_lon,ghi_dict = resample_data(ghi_rows,ghi_head)

        # Skip initial dates until we get a valid file because we need lat,lon
        if latvec is None and dni is not None:
            latvec = dni_lat
            lonvec = dni_lon
            miss = dni_dict['missing']
        if latvec is None and ghi is not None:
            latvec = ghi_lat
            lonvec = ghi_lon
            miss = ghi_dict['missing']
        if latvec is None:
            missed_dates.append(dt)
            continue

        # Initialise the netcdf object
        if ncobj is None:
            dni_var = 'solar_dni'
            ghi_var = 'solar_ghi'
            obs_var = 'obs_time'
            ncobj = nb.nc_open(froot+'.nc','w')
            #ncobj = nb.nc_open(froot+'.nc','w',format='NETCDF4_CLASSIC')
            nb.nc_set_timelatlon(ncobj,None,len(latvec),len(lonvec))
            nb.nc_set_var(ncobj,dni_var,fill=miss)
            nb.nc_set_var(ncobj,ghi_var,fill=miss)
            nb.nc_set_var(ncobj,obs_var,fill=miss)
            #nb.nc_set_var(ncobj,dni_var,fill=miss,zlib=True)
            #nb.nc_set_var(ncobj,ghi_var,fill=miss,zlib=True)
            #nb.nc_set_var(ncobj,obs_var,fill=miss,zlib=True)
            nb.nc_set_var(ncobj,'crs',dims=(),dtype="i4")  # Grid mapping container
            nb.nc_add_data(ncobj,'latitude',latvec)
            nb.nc_add_data(ncobj,'longitude',lonvec)

            dni_nc = ncobj.variables[dni_var]
            ghi_nc = ncobj.variables[ghi_var]
            obs_nc = ncobj.variables[obs_var]

        # Add observation time layers for any missed dates
        for mi,md in enumerate(missed_dates):
            obs_data = create_obs_time(latvec,lonvec,md)
            #nb.nc_add_data(ncobj,obs_var,obs_data,
            #               index=(slice(mi,None,None),latslice,lonslice))
            obs_nc[mi,:,:] = obs_data
            missed_dates = []

        # Calculate the observation time layer
        obs_data = create_obs_time(latvec,lonvec,dt)
        obs_nc[dti,:,:] = obs_data

        # Add data
        if dni is not None:
            #nb.nc_add_data(ncobj,dni_var,dni_data,index=index)
            dni_nc[dti,:,:] = dni_data
        if ghi is not None:
            #nb.nc_add_data(ncobj,ghi_var,ghi_data,index=index)
            ghi_nc[dti,:,:] = ghi_data
         
    # Add time values
    nb.nc_add_time(ncobj,dates)

    # Create an ordered metadata dictionary
    meta = create_meta(dates[0],dates[-1],dni_dict)

    # Output the file-specific metadata to a JSON file regardless of the
    # requested output formats
    jh.json_dump(meta, froot+'.json')
    
    # Setting attributes here is optional in this example because it is
    # known that netcdf_json_wrapper.py will be called next with
    # 'froot+".json"' as one of the input files.
    nb.nc_set_attributes(ncobj,meta)
    print "Added attributes"
    nb.nc_close(ncobj)
    print "Wrote: "+froot+".nc"
    return froot


def get_dates(year,month):
    # There are no ascii files for 12-17UT.
    # A month of files extend from 18UT the day before the first day of the
    # month to 11 UT on the last day of the month.
    # 18UT the day before
    d1 = datetime(year,month,1,0) - timedelta(hours=6)
    d2y = year
    d2m = month + 1
    if d2m == 13:
        d2y = d2y + 1
        d2m = 1
    # 12UT the day before
    d2 = datetime(d2y,d2m,1,0) - timedelta(hours=12)
    dt = []
    while d1 < d2:
        if d1.hour <= 11 or d1.hour >= 18: dt.append(d1)
        d1 = d1+timedelta(hours=1)
    return dt


def get_solar_files(year,dt):
    base = '/data/remotesensing/MTSAT-BoM/ascii-20121220'
    dni_path = os.path.join(base, 'time_series_hourly_dni')
    ghi_path = os.path.join(base, 'time_series_hourly_ghi')
    dtstr = dt.strftime('%Y%m%d_%HUT')
    dni = os.path.join(dni_path,year,'solar_dni_'+dtstr+'.txt')
    ghi = os.path.join(ghi_path,year,'solar_ghi_'+dtstr+'.txt')
    if not os.path.exists(dni): dni = None
    if not os.path.exists(ghi): ghi = None
    return dni,ghi


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print len(sys.argv)
        print "Usage:"
        print "  ", sys.argv[0], "year month [out_path]"
        exit()
    else:
        outpath = '.'
        year,month = int(sys.argv[1]),int(sys.argv[2])
        if len(sys.argv) >= 4:
            outpath = sys.argv[3]
            if not os.path.exists(outpath): os.makedirs(outpath)
            if not os.path.exists(outpath):
                exit("Could not create path", outpath)
        outpath = re.sub('/$','',outpath)
        dates = get_dates(year,month)

        outroot = 'solar_dni_ghi_{:04d}{:02d}'.format(year,month)
        outroot = os.path.join(outpath,outroot)

        bom_ascii_to_nc(year,dates,outroot)

