# Convert and LPDAAC mosiac HDF4 file to NetCDF4

import numpy as np
from pyhdf.SD import SD,SDC
import netcdf_builder as nb
import json_handler as jh
import re
import os
from datetime import datetime
from glob import glob
from collections import OrderedDict
import time
 

defdebug = 0
deftime = 1

def sds2var(sds):
    """
    Rearrange and SDS name to a valid NetCDF variable name
    """
    debug = defdebug
    parts = sds.split('_')
    if re.match('\d',parts[1]):
        ncparts = parts[3:] + parts[1:3] + [parts[0]]
    else:
        ncparts = parts[1:] + [parts[0]]
    # Special case reworking case for MCD43[AB]3
    if ncparts[0] == 'albedo':
        ncparts[0:2] = [ncparts[1],'sky','albedo']
    ncvar = '_'.join(ncparts)
    if debug: print sds,'->',ncvar
    return ncvar

def parse_hdfname(fname):
    debug = defdebug
    labels = ['product','year','doy','region','collection','band','sds']
    parts = os.path.basename(fname).split('.')
    d = dict(zip(labels, parts[:len(labels)]))
    d['datetime'] = datetime.strptime(d['year']+d['doy'],"%Y%j")
    if debug: print d
    return d
    
def get_hdfobj(fname):
    return SD(fname,SDC.READ)

def get_hdfdimvecs(hdfobj):
    debug = defdebug
    latvec = hdfobj.select('latitude')[:]
    lonvec = hdfobj.select('longitude')[:]
    if debug: print latvec,lonvec
    return latvec,lonvec

def get_hdfvar(hdfobj,sds):
    debug = defdebug
    # Data are unscaled (raw)
    data = hdfobj.select(sds)[:]
    if debug: print "max:",data.max(),"min:",data.min()
    return data

def get_hdfvardict(hdfobj,sds):
    debug = defdebug
    #keep = ['_FillValue', 'valid_range', 'add_offset', 'long_name', 'units', 'scale_factor', bitmap or pixel indices]
    ignore = ['HorizontalDatumName', 'scale_factor_err', 'calibrated_nt', 'add_offset_err']
    atts = hdfobj.select(sds).attributes()
    for k in atts.keys():
        if k in ignore: del atts[k]
    if debug: print atts
    return atts

def fix_scaleoffset(d):
    # Translate scale, offset and valid_range to NetCDF format
    #    HDF: phys = (stored - offset) * scale
    #    NC: phys = stored * scale + offset
    #    => scale_nc = scale_hdf
    #    => offset_nc = -scale_hdf * offset_hdf
    # valid_range is stored in raw form, so no translation required
    if 'add_offset' in d:
        d['add_offset'] = str(-float(d['scale_factor'])*float(d['add_offset']))
    else:
        d['add_offset'] = "0"
    return d

def set_ncobj(ncfile,var2miss,var2type,timevec,latvec,lonvec):
    debug = defdebug
    ncobj = nb.nc_open(ncfile,'w',format='NETCDF4')
    nb.nc_set_timelatlon(ncobj,None,len(latvec),len(lonvec))
    for vn,miss in var2miss.items():
        if debug: print vn,miss,var2type[vn]
        nb.nc_set_var(ncobj,vn,dtype=var2type[vn],fill=miss,zlib=True)
    nb.nc_set_var(ncobj,'crs',dims=(),dtype="i4")  # Grid mapping container
    nb.nc_add_time(ncobj,timevec)
    nb.nc_add_data(ncobj,'latitude',latvec)
    nb.nc_add_data(ncobj,'longitude',lonvec)
    return ncobj

def set_ncmeta(timevec,latvec,lonvec):
    meta = {}
    return meta

def update_ncmeta(meta,ncvar,hdfdict):
    return meta

def get_nctype(hdftype):
    typetab = {
           SDC.CHAR:    'CHAR',
           SDC.CHAR8:   'CHAR8',
           SDC.UCHAR8:  'UCHAR8',
           SDC.INT8:    'INT8',
           SDC.UINT8:   'UINT8',
           SDC.INT16:   'INT16',
           SDC.UINT16:  'UINT16',
           SDC.INT32:   'INT32',
           SDC.UINT32:  'UINT32',
           SDC.FLOAT32: 'FLOAT32',
           SDC.FLOAT64: 'FLOAT64'
           }
    nctypes = {
           'FLOAT32':	'f4',
           'FLOAT64':	'f8',
           'INT32':	'i4',
           'INT16':	'i2',
#           'INT64':	'i8',
           'INT8':	'i1',
           'UINT8':	'u1',
           'UINT16':	'u2',
           'UINT32':	'u4',
           'UINT64':	'u8',
           'CHAR':	'S1'
           }
    if typetab[hdftype] not in nctypes:
        raise ValueError("HDF type not in NetCDF type list: "+typetab[hdftype])
    return nctypes[typetab[hdftype]]

def lpdaac_to_nc(hdfpath):
    debug = defdebug
    hdfnames = glob(hdfpath+os.sep+'M*hdf')
    latvec = None
    var2miss = OrderedDict()
    var2type = OrderedDict()
    for hdfname in hdfnames:
        hdfparts = parse_hdfname(hdfname)
        hdfobj = get_hdfobj(hdfname)
        if latvec is None:
            latvec,lonvec = get_hdfdimvecs(hdfobj)
        hdfdict = get_hdfvardict(hdfobj,hdfparts['sds'])
        hdftype = hdfobj.datasets()[hdfparts['sds']][2]
        hdfobj.end()
        ncvar = sds2var(hdfparts['sds'])
        # Need number type as well
        var2type[ncvar] = get_nctype(hdftype)
        var2miss[ncvar] = None
        if '_FillValue' in hdfdict: var2miss[ncvar] = hdfdict['_FillValue']

    stem = '{product}.{collection}.{region}.{year}{doy}'.format(**hdfparts)
    ncfile = stem+'.nc'
    jsonfile = stem+'.json'
    if os.path.exists(ncfile): os.remove(ncfile)
    ncobj = set_ncobj(ncfile,var2miss,var2type,[hdfparts['datetime']],latvec,lonvec)
    meta = set_ncmeta([hdfparts['datetime']],latvec,lonvec)

    for hdfname in hdfnames:
        t0 = time.clock()
        if debug: print hdfname
        hdfparts = parse_hdfname(hdfname)
        hdfobj = get_hdfobj(hdfname)
        hdfvar = get_hdfvar(hdfobj,hdfparts['sds'])
        hdfdict = get_hdfvardict(hdfobj,hdfparts['sds'])
        hdfobj.end()
        ncvar = sds2var(hdfparts['sds'])
        ncobj.variables[ncvar][0,:,:] = hdfvar
        #nb.nc_add_data(ncobj,ncvar,hdfvar,(0,':',':'))
        meta = update_ncmeta(meta,ncvar,hdfdict)
        if deftime: print "Time ({}): {}".format(ncvar,time.clock()-t0)

    #jh.json_dump(jsonfile,meta)
    #nb.nc_set_attributes(ncobj,meta)
    nb.nc_close(ncobj)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print "Usage:"
        print "  ", sys.argv[0], "path/to/epoch/directory"
        exit()
    else:
        hdfpath = sys.argv[1]
        if hdfpath[-1] == os.sep: hdfpath = hdfpath[:-1]
        t0 = time.clock()
        lpdaac_to_nc(hdfpath)
        if deftime: print "Time ({}): {}".format(hdfpath,time.clock()-t0)
