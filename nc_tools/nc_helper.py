#!/usr/bin/python

import netCDF4 as nc4
import numpy as np
import os
import re
from collections import namedtuple

def get_ncobj(path_or_files):
    """
    Return a netCDF4 object. Path_or_files can be a URL, individual
    file on the filesystem, or a path and shell regex for a set of files on
    the file system.
    """
    if re.match('http',path_or_files):
        return nc4.Dataset(path_or_files)
    if os.path.isfile(path_or_files):
        return nc4.Dataset(path_or_files)
    return nc4.MFDataset(path_or_files)

def get_timedim(ncobj,tname):
    """
    Return a datetime version of a netcdf object's time dimension.
    Requires time dimension to have units and calendar attributes.
    """
    time = ncobj.variables[tname]
    dt = nc4.num2date(time[:],units=time.units,calendar=time.calendar)
    return dt

def get_timedim_allfiles(files,verbose=False):
    """
    Concatenate the time dimension for all files in the given list.
    Return a tuple of the times list and files index list, where
    	times[i] comes from findex[i].
    """
    tname = None
    times = []
    findex = []
    for i,f in enumerate(files):
        if verbose: print i,":",f
        ncobj = get_ncobj(f)
        if tname is None: tname = get_dimname(ncobj,'t')
        temp = get_timedim(ncobj,tname)
        times.extend(temp)
        findex.extend([i]*len(temp))
        ncobj.close()
    return times,findex
    
def subset_by_mask(ncobj,varname,tmask,ymask,xmask):
    """
    Return a numpy object given mask arrays for t,y,x dimensions.
    If any of tmask,ymask,xmask are equal to True then all values for that
    dimension are extracted.
    If any of tmask,ymask,xmask are None then no data can be extracted and
    None is returned.
    """
    masks = [tmask,ymask,xmask]
    names = ['tmask','ymask','xmask']
    for i,v in enumerate(masks):
        if v is None:
            print 'nc_helper.subset_by_mask:',names[i],'is None. Returning None.'
            return None
        elif v is True:
            name = get_dimname(ncobj,names[i][0])
            masks[i] = np.array([1]*len(ncobj.dimensions[name]),dtype=np.bool)
        elif not isinstance(v,np.ndarray):
            masks[i] = np.array(v,dtype=np.bool)
    ncsub = ncobj.variables[varname][masks[0],masks[1],masks[2]]
    return ncsub

def get_dimname(ncobj,axis):
    """
    Return the variable name for the axis code, e.g. t,z,y,x.
    """
    standard_map = {'t':'time','y':'latitude','x':'longitude'}
    axis = axis.lower()
    ncdims = ncobj.dimensions
    if axis in standard_map:
        if standard_map[axis] in ncdims: return standard_map[axis]
    for k in ncdims:
        kdict = ncobj.variables[k].__dict__
        if 'axis' in kdict:
            if kdict['axis'].lower() == axis: return str(k)
    print 'Axis not found in nc object: ',axis
    return None

def get_dims_tyx(file_or_ncobj):
    """
    Returns each of the 't','y','x' dimension names and vectors in a
    namedtuple, Dims(tname,yname,xname,tvec,yvec,xvec).
    """
    ncobj = file_or_ncobj
    if os.path.exists(ncobj): ncobj = get_ncobj(ncobj)
    tname = get_dimname(ncobj,'t')
    yname = get_dimname(ncobj,'y')
    xname = get_dimname(ncobj,'x')
    tvec = np.array(get_timedim(ncobj,tname))
    yvec = np.array(ncobj.variables[yname])
    xvec = np.array(ncobj.variables[xname])
    ncobj.close()
    Dims = namedtuple('Dims','tname yname xname tvec yvec xvec')
    return Dims(tname,yname,xname,tvec,yvec,xvec)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print """
Usage: {:s} params
Where:  params          A param.
""".format(sys.argv[0])
        exit()

