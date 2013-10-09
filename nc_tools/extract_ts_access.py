#!/usr/bin/python

import cmar_paths
import nc_helper as nch
import nc_subset as ncs
import numpy as np
import postgres_helper as pgh
from datetime import timedelta
from datetime import datetime
import os
import re
import psycopg2

# probably superseded by masks
def get_timeseries(ncobj,varnames,ntime,ilat,ilon):
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

# would be nice to work off nc attributes but still a role for adding the
# subset parameters
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

# would be nice to add to csv_helper as a writer function
def write_data(fh,dt,ts,d,varnames):
    # 'obs_time' column is minutes after dt
    # Format to yyyy-mm-dd hh:mm:ss TZ,Parameter(units)
    for i in range(len(dt)):
        ymd = dt[i]+timedelta(minutes=ts[i,-1])
        fh.write(' {:%Y-%m-%d %H:%M:%S UT}'.format(ymd))
        for j in range(len(varnames)-1):
            fh.write(', {:{w}.0f}'.format(ts[i,j],w=16))
        fh.write('\n')

# could come from a csv file, via csv_helper
def get_stations():
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

def want_vars():
    #exclude=['seg_type','base_date','base_time','valid_date','valid_time','forc_hrs','wrtn_date','wrtn_time','av_sfc_sw_dif','accum_prcp','qsair_scrn','soil_mois_cont','veg_ruff','av_sfc_sw_dir','sfc_temp','sfc_pres','mslp','topog','abl_ht','ttl_cld','uwnd10m','vwnd10m','temp_scrn','dewpt_scrn','lnd_mask','low_cld','mid_cld','hi_cld','av_temp_scrn','tmin_scrn','av_mslp','av_swirrtop','visibility','conv_cldbse_pres','conv_cldtop_pres','vis_prob','av_ttl_cld']
    return ['av_netswsfc','av_netlwsfc','av_lwsfcdown','av_swsfcdown','av_olr']



if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print """
Usage: {:s} files targetlat targetlon output [varnames]
Where:	files	path or regex to set of files.
	targetlat, targetlon
		latitude and longitude for nearest neighbour single pixel
		sampling. Can be extended to fancier subsetting with
		mask_helper routines.
	output	csv or postgres.
	varnames  Optional. List of variables names to extract.
""".format(sys.argv[0])
        exit()

    fnames = glob.glob(files)
    fnames.sort()
    # namedtuple(tname,yname,xname,tvec,yvec,xvec)
    dims = nch.get_dims_tyx(fnames[0])

    # CSV input
    #locs = csv_to_dictoflists('arena_stations.csv')
    #for i in range(locs['_row_count']):
    #locs = set_stations()
    #for sid in locs.keys():
        #ymin,ymax = mh.values_to_indices_spatial(yvec,ref=locs['lat'][i])
        #xmin,xmax = mh.values_to_indices_spatial(xvec,ref=locs['lon'][i])

    # Function input
    locs = get_stations()
    for sid in locs.keys():
        ymin,ymax = mh.values_to_indices_spatial(dims.yvec,ref=locs[sid]['lat'])
        xmin,xmax = mh.values_to_indices_spatial(dims.xvec,ref=locs[sid]['lon'])
        locs['sid']['ymask'] = mh.vector_mask(dims.yvec,ymin,ymax)
        locs['sid']['xmask'] = mh.vector_mask(dims.xvec,xmin,xmax)


    sid = '23034' # probably change to an input parameter
    tmask = True
    ymask = locs[sid]['ymask']
    xmask = locs[sid]['xmask']

    times,data = ncs.extract_ts(fnames,varnames,tmask,ymask,xmask)

    # simple lookup from numpy type to postgres types
    pg_fields = ''
    for i in data.dtype.names:
        if data.dtype[i].type is np.float32:
            pgtype = 'real'
        elif data.dtype[i].type is np.float64:
            pgtype = 'double'
        pg_fields += ', i '+pgtype

    # Postgres output
    tablename = '_'.join('access',sid,'1hr')
    conn = psycopg2.connect(database='solar_isd',user='solar',host='arena-isd-cdc.it.csiro.au',password='dummy')
    curs = conn.cusor()
    curs.execute("drop table if exists {};".format(tablename))
    curs.execute("create table {} (id serial primary key, refdate timestamp{};".format(tablename,pg_fields))
    pgh.time_pgcopy(data,curs,tablname,True)
