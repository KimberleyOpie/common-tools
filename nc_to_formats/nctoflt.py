# Convert a single layer from a netCDF variable to a flt file.
# It is assumed that the netCDF variable dimensions are longitude, 
# latitude and, optionally, a third dimension such as time or 
# channel/waveength

import sys
import os.path
import numpy
import netCDF3
import nc3_handler as nh

def nctoflt(ncfile, fltstem, varname, iz=0):
    """Main function to process a netCDF file to binary flt

       Output files have the stem name and suffix .flt and .hdr

       If varname is 3D, then iz is the index of the first dimension used
       to extract a 2D slice.

       If the latitude runs south to north, then the grid is flipped before
       being written
    """

    ncobj = nh.nc3_open(ncfile,'r')
    a = ncobj.variables[varname]
    # Copy out into a numpy array and make sure we have only
    # 2 dimensions and type float32.
    b = numpy.float32(ncobj.variables[varname])
    if len(b.shape) < 2 or len(b.shape) > 3:
        raise ValueError("Only 2D and 3D data allowed (not "+len(b.shape)+"D)")
    if len(b.shape) == 3:
        b = numpy.float32(b[iz,::,::].reshape(b.shape[1], b.shape[2]))
    fillValue = numpy.float32(ncobj.variables[varname]._FillValue)

    latvec = ncobj.variables['latitude']
    lonvec = ncobj.variables['longitude']
    lat1 = latvec[0]
    lat2 = latvec[len(latvec)-1]
    # Reverse if latitude runs South to North
    if lat1 < lat2:
        x = lat2
        lat2 = lat1
        lat1 = x
        b = b[::-1,]
    lon1 = lonvec[0]
    lon2 = lonvec[len(lonvec)-1]
    
    dlat = abs(lat1-lat2)/(len(latvec)-1)
    dlon = abs(lon2-lon1)/(len(lonvec)-1)
    xll = lon1-dlon*0.5
    yll = lat2-dlat*0.5


    fltname = fltstem+'.flt'
    if os.path.exists(fltname): os.unlink(fltname)
    b.tofile(fltname)
    f = file(fltstem+".hdr","w")
    
    f.write("ncols %d\n" % b.shape[1])
    f.write("nrows %d\n" % b.shape[0])
    f.write("xllcorner %f\n" % xll)
    f.write("yllcorner %f\n" % yll)
    f.write("cellsize %f\n" % dlon)
    f.write("NODATA_value %f\n" % fillValue)
    if sys.byteorder == "little":
        f.write("byteorder LSBFIRST\n")
    else:
        f.write("byteorder LSBLAST\n")
    f.close()
    attr = nh.nc3_get_attributes(ncobj)
    nh.nc3_close(ncobj)
    return attr

if __name__ == '__main__':

    if 0:
        ncfile = 'test/20110421_Etot.nc'
        varname = 'Etot'
    else:
        ncfile = 'test/C-bawap.D1-20110101.D2-20110101.I-P1D.V-rain_day.P-raw.DC-20110630T075313.DM-20111207T222219.nc'
        varname = 'rain_day'
    try:
        nctoflt(ncfile, 'testflt', varname)
    except ValueError as e:
        print "Wrong number of dimensions for variable=%s: %s" % (varname, e)
