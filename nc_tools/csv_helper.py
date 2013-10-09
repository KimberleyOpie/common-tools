#!/usr/bin/python

import re
import csv

def csv_to_dictoflists(sitescsv):
    header = False
    col2key = {}
    dofl = {}
    rowcount = 0
    with open(sitescsv, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row)==0: continue  # Blank lines
            if re.match('^\s*#',row[0]): continue  # Comment lines
            if not header:
                for i,v in enumerate(row):
                    v = v.strip()
                    col2key[i] = v
                    dofl[v] = []
                header = True
                continue
            for i,v in enumerate(row): dofl[col2key[i]].append(v.strip())
            rowcount += 1
    csvfile.close()
    dofl['_row_count'] = rowcount
    return dofl

def csv_to_listofdicts(sitescsv):
    header = False
    col2key = {}
    lofd = []
    with open(sitescsv,'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row)==0: continue  # Blank lines
            if re.match('^\s*#',row[0]): continue  # Comment lines
            if not header:
                for i,v in enumerate(row): col2key[i] = v.strip()
                header = True
                continue
            d = {}
            for i,v in enumerate(row):
                d[col2key[i]] = v.strip()
            lofd.append(d.copy())
    csvfile.close()
    return lofd


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print """
Usage: {:s} params
Where:	params		A param
""".format(sys.argv[0])
        exit()

    
