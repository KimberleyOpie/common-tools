#!/usr/bin/python

import netCDF4 as nc4
import numpy as np
from datetime import timedelta
from datetime import datetime
import os
import re

def get_indices(mf,lat,lon):
   ilat = (np.abs(mf.variables['latitude'][:]-lat)).argmin()
   ilon = (np.abs(mf.variables['longitude'][:]-lon)).argmin()
   return ilat,ilon

def get_mfdataset(path_or_files):
    if re.match('http',path_or_files):
        return nc4.Dataset(path_or_files)
    if os.path.isfile(path_or_files):
        return nc4.Dataset(path_or_files)
    return nc4.MFDataset(path_or_files)

def get_timedim(mf):
    time = mf.variables['time']
    dt = nc4.num2date(time[:],units=time.units,calendar=time.calendar)
    return dt
    
def get_timeseries(mf,varnames,ntime,ilat,ilon):
    ts = np.empty([ntime,len(varnames)])
    d = {}   # Make a dict of varnames and units while we're here
    for i in range(len(varnames)):
        ncvar =  mf.variables[varnames[i]]
        print varnames[i],ilat,ilon
        ts[:,i] = ncvar[:,ilat,ilon]
        d[varnames[i]] = {'units':ncvar.units,
                          'tsindex':i,
                          'long_name':ncvar.long_name,
                          'miss':ncvar._FillValue
                         }
    return ts,d

def write_header(fh,d,varnames):
    header = """# Solar direct and global irradiance measured from satellite.
# These data were made available by the Bureau of Meteorology, then converted
# to NetCDF format by CSIRO with an additional interpolated observation time
# grid.
# Year, month, day and hour come from the original BoM filename.
# Minute and second come from the obs_time layer of the NetCDF files, which is
# interpolated by latitude from the lookup table given in the BoM's metadata
# records and is dependent on the date and hour.
#
# Method:	{method}
# Station name:	{name}
# Latitude:	{latitude}
# Longitude:	{longitude}
#
# Latitude grid value:	  {latgrid:.2f}
# Latitude grid index:	  {ilat}
# Longitude grid value:	  {longrid:.2f}
# Longtitude grid index:  {ilon}
#
""".format(**d)
    fh.write(header)
    for k in varnames:
        if re.match('obs',k): continue
        params = """# {}
#    {}
#    Missing value: {}
#
""".format(k,d[k]['long_name'],d[k]['miss'])
        fh.write(params)
    fh.write('#YYYY-MM-DD hh:mm:ss TZ')
    for k in varnames:
        if re.match('obs',k): continue
        fh.write(', {}({})'.format(k,d[k]['units']))
    fh.write('\n')

    

def write_data(fh,dt,ts,d,varnames):
    # 'obs_time' column is minutes after dt
    # Format to yyyy-mm-dd hh:mm:ss TZ,Parameter(units)
    for i in range(len(dt)):
        ymd = dt[i]+timedelta(minutes=ts[i,-1])
        fh.write(' {:%Y-%m-%d %H:%M:%S UT}'.format(ymd))
        for j in range(len(varnames)-1):
            fh.write(', {:{w}.0f}'.format(ts[i,j],w=16))
        fh.write('\n')


def set_stations():
    d = {}
    d['23034'] = {'name':'AdelaideAirport','lat':-34.9524,'lon':138.5204}
    d['15590'] = {'name':'AliceSpringsAirport','lat':-23.7951,'lon':133.889}
    d['3003'] = {'name':'BroomeAirport','lat':-17.9475,'lon':122.2353}
    d['31011'] = {'name':'CairnsAero','lat':-16.8736,'lon':145.7458}
    d['91148'] = {'name':'CapeGrimRadiation','lat':-40.6817,'lon':144.6892}
    d['48237'] = {'name':'CobarAirportAws','lat':-31.5388,'lon':145.7964}
#    d['200284'] = {'name':'CocosIslandAirport','lat':-12.1892,'lon':96.8344}
    d['14015'] = {'name':'DarwinAirport','lat':-12.4239,'lon':130.8925}
    d['8051'] = {'name':'GeraldtonAirportComparison','lat':-28.7953,'lon':114.6975}
    d['12038'] = {'name':'KalgoorlieBoulderAirport','lat':-30.7847,'lon':121.4533}
    d['5007'] = {'name':'LearmonthAirport','lat':-22.2406,'lon':114.0967}
    d['36031'] = {'name':'LongreachAero','lat':-23.4397,'lon':144.2828}
    d['86282'] = {'name':'MelbourneAirport','lat':-37.6655,'lon':144.8321}
    d['76031'] = {'name':'MilduraAirport','lat':-34.2358,'lon':142.0867}
    d['26021'] = {'name':'MountGambierAero','lat':-37.7473,'lon':140.7739}
    d['39083'] = {'name':'RockhamptonAero','lat':-23.3753,'lon':150.4775}
    d['15135'] = {'name':'TennantCreekAirport','lat':-19.6423,'lon':134.1833}
    d['32040'] = {'name':'TownsvilleAero','lat':-19.2483,'lon':146.7661}
    d['72150'] = {'name':'WaggaWaggaAmo','lat':-35.1583,'lon':147.4575}
    d['16001'] = {'name':'WoomeraAerodrome','lat':-31.1558,'lon':136.8054}


    return d



if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print "Usage: {:s} path_regex fileout_prefix [param]".format(sys.argv[0])
        exit()

    params = ['solar_dni','solar_ghi']
    if len(sys.argv) >= 4: params = [sys.argv[3]]
    params = params+['obs_time']
    fileout_stem = sys.argv[2]
    print sys.argv[1]
    mf = get_mfdataset(sys.argv[1])
    dt = get_timedim(mf)
    ndt = len(dt)
    s = set_stations()
    print s.keys()
    for sid in s.keys():
        print sid
        if len(params) > 2:
            fileout = "_".join([fileout_stem,s[sid]['name']])+".txt"
        else:
            fileout = "_".join([fileout_stem,s[sid]['name'],params[0]])+".txt"
        ilat,ilon = get_indices(mf,s[sid]['lat'],s[sid]['lon'])
        ts,d = get_timeseries(mf,params,ndt,ilat,ilon)

        if not os.path.exists(fileout):
            fh = open(fileout,'w')
            d['method'] = 'nearest neighbour'
            d['name'] = s[sid]['name']
            d['latitude'] = s[sid]['lat']
            d['longitude'] = s[sid]['lon']
            d['ilat'] = ilat
            d['ilon'] = ilon
            d['latgrid'] = mf.variables['latitude'][ilat]
            d['longrid'] = mf.variables['longitude'][ilon]
            write_header(fh,d,params)
        else:
            fh = open(fileout,'a')

        write_data(fh,dt,ts,d,params)
        fh.close()
