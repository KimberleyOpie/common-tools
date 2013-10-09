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
    meta['time_coverage_duration'] = "P1D"
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


def get_basics():
    ncfile = '/data/remotesensing/MTSAT-BoM/nc2/2000/solar_dni_ghi_20000101.nc'
    jsonfile = re.sub('.nc$','.json',ncfile)
    ncobj = nb.nc_open(ncfile,'r')
    latvec = ncobj.variables['latitude'][:]
    lonvec = ncobj.variables['longitude'][:]
    #meta = nb.nc_get_attributes(ncobj)
    meta = jh.json_load(fname=jsonfile)
    nb.nc_close(ncobj)
    d = {}
    d['ymin'] = float(meta['geospatial_lat_min'])
    d['ymax'] = float(meta['geospatial_lat_max'])
    d['ystep'] = float(meta['geospatial_lat_step'])
    d['yunits'] = meta['geospatial_lat_units']
    d['xmin'] = float(meta['geospatial_lon_min'])
    d['xmax'] = float(meta['geospatial_lon_max'])
    d['xstep'] = float(meta['geospatial_lon_step'])
    d['xunits'] = meta['geospatial_lon_units']
    #d['missing'] = ncobj.variables['solar_dni']._FillValue
    d['missing'] = -999
    return latvec,lonvec,d


def bom_ascii_to_nc(year,dates,froot):

    debug = False
    latvec = None
    lonvec = None
    ncobj = None
    latslice = slice(None,None,None)
    lonslice = slice(None,None,None)
    missed_dates = []
    adict = {}

    latvec,lonvec,adict = get_basics()
    miss = adict['missing']

    for dti,dt in enumerate(dates):
        index = (slice(dti,None,None),latslice,lonslice)

        # Initialise the netcdf object
        if ncobj is None:
            if debug: print "Define ncobj:",dt
            dni_var = 'solar_dni'
            ghi_var = 'solar_ghi'
            obs_var = 'obs_time'
            ncobj = nb.nc_open(froot+'.nc','w',format='NETCDF4_CLASSIC')
            nb.nc_set_timelatlon(ncobj,None,len(latvec),len(lonvec))
            nb.nc_set_var(ncobj,dni_var,fill=miss,zlib=True)
            nb.nc_set_var(ncobj,ghi_var,fill=miss,zlib=True)
            nb.nc_set_var(ncobj,obs_var,fill=miss,zlib=True)
            nb.nc_set_var(ncobj,'crs',dims=(),dtype="i4")  # Grid mapping container
            nb.nc_add_data(ncobj,'latitude',latvec)
            nb.nc_add_data(ncobj,'longitude',lonvec)

            dni_nc = ncobj.variables[dni_var]
            ghi_nc = ncobj.variables[ghi_var]
            obs_nc = ncobj.variables[obs_var]

        # Calculate the observation time layer
        obs_data = create_obs_time(latvec,lonvec,dt)
        obs_nc[dti,:,:] = obs_data
         
    # Add time values
    nb.nc_add_time(ncobj,dates)

    # Create an ordered metadata dictionary
    meta = create_meta(dates[0],dates[-1],adict)

    # Output the file-specific metadata to a JSON file regardless of the
    # requested output formats
    jh.json_dump(meta, froot+'.json')
    
    # Setting attributes here is optional in this example because it is
    # known that netcdf_json_wrapper.py will be called next with
    # 'froot+".json"' as one of the input files.
    nb.nc_set_attributes(ncobj,meta)
    if debug: print "Added attributes"
    nb.nc_close(ncobj)
    print "Wrote: "+froot+".nc"
    return froot


def get_dates(year,month,day):
    # There are no ascii files for 12-17UT.
    # A month of files extend from 18UT the day before the first day of the
    # month to 11 UT on the last day of the month.
    # 18UT the day before
    d1 = datetime(year,month,day,0) - timedelta(hours=6)
    # 12UT the day
    d2 = datetime(year,month,day,12)
    dt = []
    while d1 < d2:
        if d1.hour <= 11 or d1.hour >= 18: dt.append(d1)
        d1 = d1+timedelta(hours=1)
    return dt


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4:
        print len(sys.argv)
        print "Usage:"
        print "  ", sys.argv[0], "year month day [out_path]"
        exit()
    else:
        outpath = '.'
        year,month,day = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
        if len(sys.argv) >= 5:
            outpath = sys.argv[4]
            if not os.path.exists(outpath): os.makedirs(outpath)
            if not os.path.exists(outpath):
                exit("Could not create path", outpath)
        outpath = re.sub('/$','',outpath)
        dates = get_dates(year,month,day)

        outroot = 'solar_dni_ghi_{:04d}{:02d}{:02d}'.format(year,month,day)
        outroot = os.path.join(outpath,outroot)

        bom_ascii_to_nc(year,dates,outroot)

