import sys
import numpy
import netCDF3
import nc3_handler as nh
from datetime import datetime

# Given netCDF data files with two lat/lon dims and a redundant (n=1)
# time dimension, average the contents in time, with due regard to _FillValue
# and write the output as a netcdf file of the same shape.  Return the
# attributes from each file as a dict so that the caller can select
# which metadata to propagate.

def ncaverage(varname, files, outname):
    """
    Create a netCDF file OUTNAME by averaging the lon/lat grids stored in
    VARNAME in the input files FILES.
    Return a dictionary of the attributes from the input files
    so that the caller can create addiitonal application-specific metadata
    to insert into the output file.

    Checks are performed to make sure that input units are either consistent
    or non-existent, and that the input grid shapes are the same, and the
    endpoints of the coordinate variables match.
    """
    history = []
    now = datetime.now()
    history.append("ncaverage at %04d-%02d-%02dT%02d:%02d:%02d, files..." % (now.year, now.month, now.day, now.hour, now.minute, now.second))
    units = None
    attributes = dict();
    first = 1
    for file in files:

        print "File: "+file
        history.append(" + "+file)
        ncobj = nh.nc3_open(file,'r')
        attr = nh.nc3_get_attributes(ncobj)
        if varname+':units' in attr:
            newunit = attr[varname+':units']
            if units:
                if newunit != units:
                    raise ValueError("Inconsistent units ("+newunits+") on var "+varname+" in file "+file+", previous files were ("+units+")")
            units = newunit
        attributes[file] = attr

        # We need these for the first file and then on each subsequent file
        a = ncobj.variables[varname]
        latvec = ncobj.variables['latitude']
        lonvec = ncobj.variables['longitude']
        tvec = ncobj.variables['time']

        # The first file is used to 1) set the shape of the output
        # 2) capture some basic metadata used in writing the output file,
        # which we start to write here while we've got access to the 
        # dimension variables.
        if first:
            first = 0
            mean = numpy.zeros(a.shape, dtype=numpy.float32)
            n = numpy.zeros(a.shape, dtype=numpy.int32)
            if 'time:calendar' in attributes[file]:
                t1 = netCDF3.num2date(tvec[0],units=tvec.units,calendar=tvec.calendar)
            else:
                t1 = netCDF3.num2date(tvec[0],units=tvec.units)
            t2 = t1
            ncout = nh.nc3_open(outname,'w')
            nh.nc3_set_timelatlon(ncout,1,len(latvec),len(lonvec))
            nh.nc3_add_data(ncout,'latitude',latvec)
            nh.nc3_add_data(ncout,'longitude',lonvec)
            # remember the endpoints of each CV and compute the step size for each
            lat1 = latvec[0]
            lat2 = latvec[len(latvec)-1]
            lon1 = lonvec[0]
            lon2 = lonvec[len(lonvec)-1]
            dlat = abs(lat2-lat1)/len(latvec)
            dlon = abs(lon2-lon1)/len(lonvec)
       

        # Guard against incompatible files - shapes must be the 
        # same AND CV endpoint values must match within a fraction of a pixel.
        if a.shape != mean.shape:
            raise ValueError("Data in "+file+" is a different shape from first file")
        if (abs(latvec[0]-lat1) > 0.01*dlat or
            abs(latvec[len(latvec)-1]-lat2) > 0.01*dlat or
            abs(lonvec[0]-lon1) > 0.01*dlon or
            abs(lonvec[len(lonvec)-1]-lon2) > 0.01*dlon):
            raise ValueError("Coordinate endpoints in "+file+" do not match first file")

        # Keep track of the earliest and latest times
        if 'time:calendar' in attributes[file]:
            t = netCDF3.num2date(tvec[0],units=tvec.units,calendar=tvec.calendar)
        else:
            t = netCDF3.num2date(tvec[0],units=tvec.units)
        if t < t1: t1 = t
        if t > t2: t2 = t

        # b is 1 unless there is a missing value, where it is 0
        fillValue = ncobj.variables[varname]._FillValue
        b = numpy.where(a == fillValue, 0, 1)
        n = n+b
        mean = mean+b*a
        nh.nc3_close(ncobj)

    # Before computing mean, set all values with 0 count to _FillValue, and
    # the corresponding count to 1.  This avoids division by zero and leaves 
    # those pixels set to _FillValue.
    mean = numpy.where(n==0, fillValue, mean)
    n = numpy.where(n==0, 1, n)
    mean = mean/n

    nh.nc3_set_var(ncout,varname)
    nh.nc3_add_time(ncout,[t1])
    nh.nc3_add_data(ncout,varname,mean)
    attr = {}
    attr['history'] = '\n'.join(history)
    attr[varname+':_FillValue'] = fillValue
    if units:
        attr[varname+':units'] = units
    nh.nc3_set_attributes(ncout,attr)
    nh.nc3_close(ncout)
    return attributes

if __name__ == '__main__':
    files = ('test/C-bawap.D1-20110101.D2-20110101.I-P1D.V-rain_day.P-raw.DC-20110630T075313.DM-20111207T222219.nc', 'test/C-bawap.D1-20110102.D2-20110102.I-P1D.V-rain_day.P-raw.DC-20110701T075318.DM-20111207T222304.nc', 'test/C-bawap.D1-20110103.D2-20110103.I-P1D.V-rain_day.P-raw.DC-20110702T074008.DM-20111207T222346.nc', 'test/C-bawap.D1-20110104.D2-20110104.I-P1D.V-rain_day.P-raw.DC-20110703T080425.DM-20111207T222430.nc')
    var = 'rain_day'
    files = ('test/20110421_Etot.nc', 'test/20110422_Etot.nc', 'test/20110420_Etot.nc')
    var = 'Etot'
    try:
        attrs = ncaverage(var, files, "testavg.nc")
    except ValueError as e:
        print "Incompatible files: %s" % e
        sys.exit(1)

    for k in attrs.keys():
        print '----- Attributes -----'
        print k
        a = attrs[k]
        for w in a.keys():
            print " "+w

    sys.exit(0)
