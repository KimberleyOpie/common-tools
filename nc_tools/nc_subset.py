#!/usr/bin/python

import nc_helper as nch
import numpy as np
#import netCDF4 as nc4
#import datetime as dt
#import os
#import re
#import csv
#from collections import OrderedDict

def extract_ts(files,varnames,tmask,ymask,xmask, \
               maskfunc=None,verbose=False,times=None):
    """
    Extract timseries data from a set of nc files and return a corresponding
    times vector and a data structure array.

    Inputs:
    files	List of filenames, e.g., from glob.glob(path_regex).sort().
    varnames	List of variables names to extract data for.
    tmask	Boolean array corresponding to the t dimension. This is used
    		by nc_helper.subset_by_mask(). The tmask can be tricky e.g.
    		for files with different length unlimited dimensions. Use of a
    		mask is based on position, whereas for an unlimited dimension
    		you will probably want to filter based on value. If this is
    		case then set tmask to True and manage additional filtering
    		outside this function.
    		If tmask is True then all times from all files will be
                extracted as per nc_helper.subset_by_mask().
    ymask,xmask	Boolean arrays corresponding to the y,x dimensions
    		respectively. These are used by nc_helper.subset_by_mask().
    maskfunc	Optional. Name of a function to apply to the output of
    		nc_helper.subset_by_mask(), e.g. to calculate the mean of a
    		spatial mask.
    verbose	Optional. If True print some useful information.
    times	Optional. Provide a times array as previously computed by
    		extract_ts. Obviously the files and tmask parameters need to
    		be the same as for when times was generated to give
    		meaningful results.

    Outputs:
    times	List of datetime objects corresponding to the rows of data.
    data	Structure array with nrows==len(times) and
    		ncolumns==len(varnames). The dtype of the structure is built
    		from the first nc file and varnames such that
    			data[varname] = timseries of varname
    			data.dtype.names = [varnames]
    			data.dtype.fields = {varnames}
    			data.dtype.fields[varname][title] = ncvar.__dict__
    """
    # Get variable names from first file
    ncobj = nch.get_ncobj(files[0])
    ncvars = ncobj.variables
    tname = nch.get_dimname(ncobj,'t')

    # Build the dtype for the requested variable
    # Dtype is of the form ((ncvar.dict, varname), ncvar.dtype)
    dtype = [('reftime','string?')] # Can include times in data structure?
    dtype = []  # Otherwise..
    for v in varnames[:]:
        if v not in ncvars:
            print 'Requested variable,',v,'not in first nc object. Excluded.'
            varnames.remove(v)
            continue
        dtype += [((ncvars[v].__dict__,v),ncvars[v].dtype)]

    # Close the first file
    ncobj.close()

    # Create times numpy array
    do_time = True
    if isinstance(times,np.ndarray):
        do_time = False
        ntimes = len(times)
    elif isinstance(tmask,(list,tuple,np.ndarray)):
        ntimes = len(files)*sum(tmask)  # True values
        times = np.zeros(ntimes, object)  # datetime objects
    elif tmask is True:
        if verbose: print 'Retrieving the time dimension for all files'
        times,tindex = nch.get_timedim_allfiles(files,verbose)
        times = np.array(times)
        ntimes = len(times)
    else:
        print 'Unknown type or value of tmask:',tmask,type(tmask)
        return None,None

    # Create data numpy array
    data = np.zeros(ntimes, dtype)
    if verbose:
        print 'num files =',len(files)
        print 'num variables =',len(varnames)
        print 'times.shape =',times.shape
        print 'data.shape =',data.shape
        print 'data.dtype =',data.dtype

    # Work through each file
    counter = 0
    for fname in files:
        if verbose:
            print 'file:',fname
        ncobj = nch.get_ncobj(fname)
        # Add to the times array if required
        if do_time and isinstance(tmask,(list,tuple,np.ndarray)):
            temp = nch.get_timedim(ncobj,tname)
            temp = [x for i,x in enumerate(temp) if tmask[i]]
            times[counter:counter+len(temp)] = temp
        # Add to the data array
        for v in varnames:
            temp = nch.subset_by_mask(ncobj,v,tmask,ymask,xmask)
            if maskfunc is not None: temp = maskfunc(temp)
            temp = temp.flatten()
            if verbose: print 'Num zero values:',v,sum(temp==0)
            data[v][counter:counter+len(temp)] = temp
        # Cleanup before next iteration
        ncobj.close()
        counter += len(temp)

    # We're done
    return times,data



if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print """
Usage: {:s} nc_path_regex fileout_prefix sites_csv [variables]
Where:	nc_path_regex	Path and regex for the input netcdf files.
	fileout_prefix	Output file prefix including path.
	sites_csv	CSV file of sites. Columns can be in any order.
			sitename,minlat,maxlat,minlon,maxlon,mintime,maxtime
	variables	Optional. Default is all variables.
			Comma separated list of variables.
""".format(sys.argv[0])
        exit()

    
