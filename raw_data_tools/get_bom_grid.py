# Fetch BoM gridded files from the Bureau server.

import sys
from datetime import datetime
import calendar
import subprocess
import os.path

R = 'http://www.bom.gov.au/web03/ncc/www/awap'
WGET = '/usr/bin/wget'

# Append yyyymmddyyyymmdd.grid.Z to these URLs to get the file
vars = {
  'rain/day' : R+'/rainfall/totals/daily/grid/0.05/history/nat',
  'rain_rmse/day' : R+'/rainfall/rmse/daily/grid/0.05/history/nat',
  'rain_recal/day' : R+'/rainfall/totals/daily/grid/0.05/history/nat_recal',
  'tmax/day' : R+'/temperature/maxave/daily/grid/0.05/history/nat',
  'tmax_rmse/day' : R+'/temperature/maxrmse/daily/grid/0.05/history/nat',
  'tmin/day' : R+'/temperature/minave/daily/grid/0.05/history/nat',
  'tmin_rmse/day' : R+'/temperature/minrmse/daily/grid/0.05/history/nat',
  'rad/day' : R+'/solar/solarave/daily/grid/0.05/history/nat',
  'vph09/day' : R+'/vprp/vprph09/daily/grid/0.05/history/nat',
  'vph15/day' : R+'/vprp/vprph15/daily/grid/0.05/history/nat',
  'rain/month' : R+'/rainfall/totals/month/grid/0.05/history/nat',
  'rain_rmse/month' : R+'/rainfall/rmse/month/grid/0.05/history/nat',
  'tmax/month' : R+'/temperature/maxave/month/grid/0.05/history/nat',
  'tmax_rmse/month' : R+'/temperature/maxrmse/month/grid/0.05/history/nat',
  'tmin/month' : R+'/temperature/minave/month/grid/0.05/history/nat',
  'tmin_rmse/month' : R+'/temperature/minrmse/month/grid/0.05/history/nat',
  'rad/month' : R+'/solar/solarave/month/grid/0.05/history/nat',
  'vph09/month' : R+'/vprp/vprph09/month/grid/0.05/history/nat',
  'vph15/month' : R+'/vprp/vprph15/month/grid/0.05/history/nat'
}

filePrefix = {
  'rain/day' : 'rain_daily',
  'rain_rmse/day' : 'rain_rmse_daily',
  'rain_recal/day' : 'rain_recal_daily',
  'tmax/day' : 'tmax_daily',
  'tmax_rmse/day' : 'tmax_rmse_daily',
  'tmin/day' : 'tmin_daily',
  'tmin_rmse/day' : 'tmin_rmse_daily',
  'rad/day' : 'rad_daily',
  'vph09/day' : 'vph09_daily',
  'vph15/day' : 'vph15_daily',
  'rain/month' : 'rain_month',
  'rain_rmse/month' : 'rain_rmse_month',
  'tmax/month' : 'tmax_month',
  'tmax_rmse/month' : 'tmax_rmse_month',
  'tmin/month' : 'tmin_month',
  'tmin_rmse/month' : 'tmin_rmse_month',
  'rad/month' : 'rad_month',
  'vph09/month' : 'vph09_month',
  'vph15/month' : 'vph15_month',
}

def get_bom_grid(var, yyyymmdd, month):
    """
    Given a variable VAR and date YYYYMMDD, fetch the .grid file from 
    the Bureau server using wget. If MONTH != 0, monthly files are fetched
    rather than daily.  Return the name of the file fetched.
    """

    monthday = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    interval = 'day' if month == 0  else 'month'
    key = var+'/'+interval
    if key not in vars:
        raise ValueError("No such var ("+var+")")

    try:
        s = datetime.strptime(yyyymmdd, "%Y%m%d") 
    except ValueError:
        raise ValueError("Invalid date ("+yyyymmdd+")")

    if month:
        if s.month == 2 and calendar.isleap(s.year): monthday[2] = 29
        fname = "%04d%02d01%04d%02d%02d.grid.Z" % (s.year, s.month, s.year, s.month, monthday[s.month-1])
    else:
        fname = yyyymmdd+yyyymmdd+".grid.Z"
    out = filePrefix[key]+'_'+fname

    if os.path.exists(out): 
        raise RuntimeError("File already exists ("+out+")")

    url = vars[key]+"/"+fname
    if subprocess.call([WGET, '--output-document='+out, url]): 
        os.remove(out)
        raise IOError("Calling %s --output-document=%s %s" % (WGET, out, url))
    return out

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print "USAGE: python get_bom_grid.py VAR YYYYMMDD MONTH(1/0)"
        sys.exit(1)

    try:
        fname = get_bom_grid(sys.argv[1], sys.argv[2], int(sys.argv[3]))
        print "Fetched file: "+fname
    except ValueError as e:
        print "Aborted - erroneous arguments: %s" % e
        sys.exit(1)
    except IOError as e:
        print "URL not found: %s" % e
        sys.exit(1)
    except RuntimeError as e:
        print "File of same name present: %s" % e
        sys.exit(1)
    sys.exit(0)
