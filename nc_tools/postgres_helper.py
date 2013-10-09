#!/usr/bin/python

import psycopg2
import numpy as np
from struct import pack
from io import BytesIO
#import ppygis

def prepare_text(dat):
    """
    Prepare a numpy record array for copying to postgres as text. Text should
    be slower than the binary version, particularly if your data is originally
    numpy arrays. If using text, also consider psycopg2 and postgres's
    copy_expert("COPY test TO STDOUT WITH CSV HEADER") and
    http://www.postgresql.org/docs/current/interactive/sql-copy.html#AEN66712.

    Code borrowed from http://stackoverflow.com/questions/8144002/use-binary-copy-table-from-with-psycopg2/8150329#8150329

    Dat is a numpy record array. Rows and columns are synomynous with the
    destination postgres table but the dtype is a record or structure (name
    and format per column).
    """
    cpy = BytesIO()
    for row in dat:
        cpy.write('\t'.join([repr(x) for x in row]) + '\n')
    return(cpy)

def prepare_binary(dat):
    """
    Prepare a numpy record array for copying to postgres as binary. This should
    be quicker than text, particularly for data already stored in numpy arrays.

    Code borrowed from http://stackoverflow.com/questions/8144002/use-binary-copy-table-from-with-psycopg2/8150329#8150329

    Dat is a numpy record array. Rows and columns are synomynous with the
    destination postgres table but the dtype is a record or structure (name
    and format per column).
    """
    pgcopy_dtype = [('num_fields','>i2')]
    for field, dtype in dat.dtype.descr:
        if isinstance(field,tuple): field = field[1]
        pgcopy_dtype += [(field + '_length', '>i4'),
                         (field, dtype.replace('<', '>'))]
    pgcopy = np.empty(dat.shape, pgcopy_dtype)
    pgcopy['num_fields'] = len(dat.dtype)
    for i in range(len(dat.dtype)):
        field = dat.dtype.names[i]
        pgcopy[field + '_length'] = dat.dtype[i].alignment
        pgcopy[field] = dat[field]
    cpy = BytesIO()
    # Required header for postgres binary copy. See http://www.postgresql.org/docs/current/interactive/sql-copy.html#AEN66756
    cpy.write(pack('!11sii', b'PGCOPY\n\377\r\n\0', 0, 0))
    cpy.write(pgcopy.tostring())  # all rows
    cpy.write(pack('!h', -1))  # file trailer
    return(cpy)

def time_pgcopy(dat, curs, table, binary):
    """
    Copy data in a numpy array directly into a postgres table.
    See http://stackoverflow.com/questions/8144002/use-binary-copy-table-from-with-psycopg2/8150329#8150329

    Use:	time_pgcopy(data, cursor, 'num_data_text', binary=False)
		time_pgcopy(data, cursor, 'num_data_binary', binary=True)
    """
    print('Processing copy object for ' + table)
    tstart = datetime.now()
    if binary:
        cpy = prepare_binary(dat)
    else:  # text
        cpy = prepare_text(dat)
    tendw = datetime.now()
    print('Copy object prepared in ' + str(tendw - tstart) + '; ' +
          str(cpy.tell()) + ' bytes; transfering to database')
    cpy.seek(0)
    if binary:
        curs.copy_expert('COPY ' + table + ' FROM STDIN WITH BINARY', cpy)
    else:  # text
        curs.copy_from(cpy, table)
    conn.commit()
    tend = datetime.now()
    print('Database copy time: ' + str(tend - tendw))
    print('        Total time: ' + str(tend - tstart))
    return

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print """
Usage: {:s} param1
Where:  param1		Some parameter.
""".format(sys.argv[0])
        exit()

