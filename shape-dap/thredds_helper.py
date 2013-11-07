#!/usr/bin/python

# Define a simple structure to help navigate different thredds directory
# structures.

import sys

#nci = dict()
nci = dict()

nci['lpdaac_mosaics_cmar'] = {
  'path':'http://thredds0.nci.org.au/thredds/catalog/u39/modis/lpdaac-mosaics-cmar/v1-hdf4/catalog.html',
  'latdim':'latitude',
  'londim':'longitude',
  'timedim':{'fname':'\/(\d{4})\.(\d\d)\.(\d\d)\/'},
  'varnames':{'file':['StructMetadata%2E0', 'OldArchiveMetadata%2E0', 'latitude', 'longitude']}}


#http://thredds0.nci.org.au/thredds/catalog/u39/wirada/cmrset/v20120328-bias-corrected/catalog.html
#http://thredds0.nci.org.au/thredds/catalog/u39/modis/fractionalcover-clw/catalog.html
#http://thredds0.nci.org.au/thredds/catalog/u39/modis/burnt-area-cdu/catalog.html
#http://thredds0.nci.org.au/thredds/catalog/u39/hyperion/surface-reflectance-curtin/catalog.html
#http://thredds0.nci.org.au/thredds/catalog/u39/avhrr/pal-cmar/v1/catalog.html
#http://thredds0.nci.org.au/thredds/catalog/u39/avhrr/fire-frequency-cdu/catalog.html

def thredds_nci(short=None):
    if short is None:
        return nci
    elif short in nci:
        return nci[short]
    else:
        raise KeyError("Key not found in NCI THREDDS helper definition: "+short)

if __name__ == '__main__':
    if len(sys.argv) < 1:
        print "Usage:"
        print "  ", sys.argv[0], "[short_name]"
        exit()
    else:
        try:
            short = None
            if len(sys.argv) >= 2: short = sys.argv[1]
            print thredds_nci(short)
        except:
            raise

