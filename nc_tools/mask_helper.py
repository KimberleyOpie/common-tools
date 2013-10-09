#!/usr/bin/python

import numpy_routines as npr
import numpy as np
import datetime as dt
import re

def vector_mask(vec,idx0,idx1):
    mask = [False]*len(vec)
    mask[idx0:(idx1+1)] = [True]*(idx1-idx0+1)
    return np.array(mask)

def apply_dt_format(val,fmt):
    if val is not None and not \
            isinstance(val,(dt.time,dt.date,dt.datetime,dt.timedelta)):
        # No try block because we want it to die if there's and exception
        val = dt.datetime.strptime(val,fmt)
    return val

def values_to_indices_time(vec,ref=None,less=None,more=None,
                           reffmt=None,lessfmt=None,morefmt=None):
    """
    Return min,max indices into a datetime vector.
    The use of ref,less,more depends on the combination given.
    If ref only		min=max=index(ref)
    If less only	min=index(less), max=len(vec)-1
    If more only	min=0, max=index(more)
    If ref and less	min=index(ref-less), max=len(vec)-1
    If ref and more	min=0, max=index(ref+more)
    If less and more	min=index(less), max=index(more)
    If ref,less,more	min=index(ref-less), max=index(ref+more)
    Format strings should be compatible with datetime.datetime.strptime.
    """
    idx0 = 0
    idx1 = len(vec)-1
    # Apply formats where required, raise exception on failure
    ref = apply_dt_format(ref,reffmt)
    less = apply_dt_format(less,lessfmt)
    more = apply_dt_format(more,morefmt)
    # Sanity check
    if isinstance(ref,dt.timedelta):
        print 'ref is type timedelta. Require an absolute value.'
        return None
    if isinstance(less,dt.timedelta) and ref is None:
        print 'less is type timedelta. Require an absolute value when ref=None.'
        return None
    if isinstance(more,dt.timedelta) and ref is None:
        print 'more is type timedelta. Require an absolute value when ref=None.'
        return None
    # Tests
    if ref is not None and [less,more]==[None,None]:
        idx0 = npr.index_of_nearest(vec,ref)
        idx1 = idx0
    if less is not None and ref is None:
        idx0 = npr.index_of_nearest(vec,less)
    if more is not None and ref is None:
        idx1 = npr.index_of_nearest(vec,more)
    if ref is not None and less is not None:
        if not isinstance(less,dt.timedelta):
            td0 = dt.datetime.strptime('00','%S')
            less = less-td0  # timedelta
        idx0 = npr.index_of_nearest(vec,ref-less)
    if ref is not None and more is not None:
        if not isinstance(more,dt.timedelta):
            td0 = dt.datetime.strptime('00','%S')
            more = more-td0  # timedelta
        idx1 = npr.index_of_nearest(vec,ref+more)
    if idx0 > idx1:
        tmp = idx0
        idx0 = idx1
        idx1 = tmp
    return idx0,idx1

def apply_spatial_unit(val,unit):
    if val is not None:
        if unit is None: unit = 'deg'
        if re.search('deg',unit): pass
        elif unit.lower() == 'km': val = float(val)*0.01
        elif unit.lower() == 'm':  val = float(val)*0.01/1000
        else:
            raise ValueError('unknown spatial unit: '+unit+' ['+val+']')
    return val

def values_to_indices_spatial(vec,ref=None,less=None,more=None,
                              lessunit=None,moreunit=None):
    """
    Return min,max indices into a datetime vector.
    The use of ref,less,more depends on the combination given.
    If ref only		min=max=index(ref)
    If less only	min=index(less), max=len(vec)-1
    If more only	min=0, max=index(more)
    If ref and less	min=index(ref-less), max=len(vec)-1
    If ref and more	min=0, max=index(ref+more)
    If less and more	min=index(less), max=index(more)
    If ref,less,more	min=index(ref-less), max=index(ref+more)
    Units should be one of 'km', 'm'. Default is decimal degrees (same as vec).
    """
    idx0 = 0
    idx1 = len(vec)-1
    # Apply formats where required, raise exception on failure
    # If ref is None then assume less and more are in vec units
    if ref is not None:
        less = apply_spatial_unit(less,lessunit)
        more = apply_spatial_unit(more,moreunit)
    # Sanity check
    if ref is not None: ref = float(ref)
    if less is not None: less = float(less)
    if more is not None: more = float(more)
    # Tests
    if ref is not None and [less,more]==[None,None]:
        idx0 = npr.index_of_nearest(vec,ref)
        idx1 = idx0
    if less is not None and ref is None:
        idx0 = npr.index_of_nearest(vec,less)
    if more is not None and ref is None:
        idx1 = npr.index_of_nearest(vec,more)
    if ref is not None and less is not None:
        idx0 = npr.index_of_nearest(vec,ref-less)
    if ref is not None and more is not None:
        idx1 = npr.index_of_nearest(vec,ref+more)
    if idx0 > idx1:
        tmp = idx0
        idx0 = idx1
        idx1 = tmp
    return idx0,idx1

def test_time():
    import numpy as np
    tvec = ['2009110601','2009110602','2009110603','2009110604','2009110605','2009110606']
    tvec[:] = [apply_dt_format(x,'%Y%m%d%H') for x in tvec]
    arr = np.array(tvec)
    # If ref only	min=max=index(ref)
    print 'ref in: want 3,3:',\
      values_to_indices_time(tvec,ref='200911060420',reffmt='%Y%m%d%H%M')
    print 'ref dn: want 0,0:',\
      values_to_indices_time(tvec,ref='200811060420',reffmt='%Y%m%d%H%M')
    print 'ref up: want 5,5:',\
      values_to_indices_time(tvec,ref='201211060420',reffmt='%Y%m%d%H%M')
    # If less only	min=index(less), max=len(vec)-1
    print 'less in: want 3,5:',\
      values_to_indices_time(tvec,less='200911060420',lessfmt='%Y%m%d%H%M')
    print 'less dn: want 0,5:',\
      values_to_indices_time(tvec,less='200811060420',lessfmt='%Y%m%d%H%M')
    print 'less up: want 5,5:',\
      values_to_indices_time(tvec,less='201211060420',lessfmt='%Y%m%d%H%M')
    # If more only	min=0, max=index(more)
    print 'more in: want 0,3:',\
      values_to_indices_time(tvec,more='200911060420',morefmt='%Y%m%d%H%M')
    print 'more dn: want 0,0:',\
      values_to_indices_time(tvec,more='200811060420',morefmt='%Y%m%d%H%M')
    print 'more up: want 0,5:',\
      values_to_indices_time(tvec,more='201211060420',morefmt='%Y%m%d%H%M')
    # If ref and less	min=index(ref-less), max=len(vec)-1
    print 'ref-less in: want 1,5:',\
      values_to_indices_time(tvec,ref='200911060420',reffmt='%Y%m%d%H%M',\
                                  less='0240',lessfmt='%H%M')
    # If ref and more	min=0, max=index(ref+more)
    print 'ref+more in: want 0,4:',\
      values_to_indices_time(tvec,ref='200911060420',reffmt='%Y%m%d%H%M',\
                                  more='0020',morefmt='%H%M')
    # If less and more	min=index(less), max=index(more)
    print 'less,more in: want 2,4:',\
      values_to_indices_time(tvec,less='200911060320',lessfmt='%Y%m%d%H%M',\
                                  more='200911060440',morefmt='%Y%m%d%H%M')
    # If ref,less,more	min=index(ref-less), max=index(ref+more)
    print 'ref-less,ref+more in: want 1,4:',\
      values_to_indices_time(tvec,ref='200911060420',reffmt='%Y%m%d%H%M',\
                                  less='0240',lessfmt='%H%M', \
                                  more='0020',morefmt='%H%M')

def test_spatial():
    import numpy_routines as npr
    xvec = npr.create_vector(115,10,5)
    # If ref only	min=max=index(ref)
    print 'ref in: want 3,3:',\
      values_to_indices_spatial(xvec,ref=130.2)
    print 'ref dn: want 0,0:',\
      values_to_indices_spatial(xvec,ref=100.2)
    print 'ref up: want 9,9:',\
      values_to_indices_spatial(xvec,ref=230.2)
    # If less only	min=index(less), max=len(vec)-1
    print 'ref in: want 3,9:',\
      values_to_indices_spatial(xvec,less=130.2)
    print 'ref dn: want 0,9:',\
      values_to_indices_spatial(xvec,less=100.2)
    print 'ref up: want 9,9:',\
      values_to_indices_spatial(xvec,less=230.2)
    # If more only	min=0, max=index(more)
    print 'ref in: want 0,3:',\
      values_to_indices_spatial(xvec,more=130.2)
    print 'ref dn: want 0,0:',\
      values_to_indices_spatial(xvec,more=100.2)
    print 'ref up: want 0,9:',\
      values_to_indices_spatial(xvec,more=230.2)
    # If ref and less	min=index(ref-less), max=len(vec)-1
    print 'ref-less in: want 2,9:',\
      values_to_indices_spatial(xvec,ref=130.2,less=750,lessunit='km')
    # If ref and more	min=0, max=index(ref+more)
    print 'ref+more in: want 0,5:',\
      values_to_indices_spatial(xvec,ref=130.2,more=750,moreunit='km')
    # If less and more	min=index(less), max=index(more)
    print 'less,more in: want 3,7:',\
      values_to_indices_spatial(xvec,less=130.2,more=148)
    # If ref,less,more	min=index(ref-less), max=index(ref+more)
    print 'ref-less,ref+more in: want 2,5:',\
      values_to_indices_spatial(xvec,ref=130.2,less=750,lessunit='km',\
                                               more=750,moreunit='km')


# This function does all three, t,y,x, vectors and mixes both indexes and
# values. This has been split up up into one routine for time values and one
# for spatial values. Indexes can be handled outside of this file.
# Once we're happy that the new routines do what we need then this old routine
# can be deleted.
def make_mask_old(tvec,yvec,xvec,roi):
    """
    Return vectors of the size of tvec,yvec,xvec with Boolean values
    for the elements that define the region of interest.
    ROI, a dict, may include just about any combination of the following to
    define the region of interest:
        idx_ref_x, idx_ref_y, idx_ref_t
        idx_less_x, idx_less_y, idx_less_t
        idx_more_x, idx_more_y, idx_more_t
        val_ref_x, val_ref_y, val_ref_t
        val_less_x, val_less_y, val_less_t
        val_more_x, val_more_y, val_more_t
        unit_less_x, unit_less_y
        unit_more_x, unit_more_y
        fmt_ref_t, fmt_less_t, fmt_more_t
    Where:
        idx is an index into the respective vector.
        val is a value in the respective vector.
        unit relates to the more or less values, e.g. km instead of degrees
        (possibly not fully implemented).
        fmt relates to the val_*_t values for creating a datetime value.
    Vals (and units) are converted to indices, so indices take preference.
    If no roi parameters are given for a vector then all elements of the vector
    are accepted into the roi.
    If a ref is given then less and more are relative to ref, else less and
    more are absolute.
    """
    # Convert vals (and units) to idxs
    clean = {}
    for dim in ['t','y','x']:
        vec = eval(dim+'vec')
        for rel in ['ref','less','more']:
            idx = '_'.join(['idx',rel,dim])
            val = '_'.join(['val',rel,dim])
            unit = '_'.join(['unit',rel,dim])
            if dim == 't': fmt = '_'.join(['fmt',rel,dim])
            if idx in roi:
                #clean[idx] = max(0,roi[idx])
                #clean[idx] = min(clean[idx],len(vec)-1)
                if roi[idx] < 0: clean[idx] = 0
                elif roi[idx] >= len(vec): clean[idx] = len(vec)-1
                else: clean[idx] = roi[idx]
            elif val in roi:

                # handle units
                if dim == 't':
                    # use datetime
                    try:
                        roi[val] = dt.datetime.strptime(roi[val],roi[fmt])
                    except ValueError as e:
                        print 'could not parse datetime:',roi[val],':',roi[fmt]
                        return None
                else:
                    if unit not in roi: roi[unit] = 'deg'
                    if re.search('deg',roi[unit]): pass
                    elif roi[unit] == 'km': roi[val] = roi[val]*0.01
                    elif roi[unit] == 'm':  roi[val] = roi[val]*0.01/1000
                    else:
                        print 'unknown spatial unit:',roi[unit]
                        return None

                # handle ref vs more,less
                idxref = '_'.join(['idx','ref',dim])
                valref = '_'.join(['val','ref',dim])
                if (rel == 'ref') or (idxref not in clean):
                    # ref and absolute less,more
                    clean[idx] = npr.index_of_nearest(vec,roi[val])
                else:
                    if valref not in roi:
                       roi[valref] = vec[idxref]
                    if dim == 't':
                        # use timedelta
                        td0 = dt.datetime.strptime('00','%S')
                        tdelta = roi[val]-td0  # timedelta object
                        if rel=='less': roi[val] = roi[valref]-tdelta
                        if rel=='more': roi[val] = roi[valref]+tdelta
                    else:
                        if rel=='less': roi[val] = roi[valref]-roi[val]
                        if rel=='more': roi[val] = roi[valref]+roi[val]
                    clean[idx] = npr.index_of_nearest(vec,roi[val])
    # should have clean dict with some combination of
    #   idx_ref_x, idx_ref_y, idx_ref_t
    #   idx_less_x, idx_less_y, idx_less_t
    #   idx_more_x, idx_more_y, idx_more_t
    # if none of ref,less,more the whole vector
    # if ref only then corresponding element is True
    # if less only, optionally ref, (=min) then element to end
    # if more only, optionally ref, (=max) then start to element
    # if less and more (=min,max) the element to element



if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print """
Usage: {:s} params
Where:	params		A param.
""".format(sys.argv[0])
        exit()

    
