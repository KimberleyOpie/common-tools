import numpy as np



def overlap(mxvec, myvec, txvec, tyvec, verbose=False):
    """
    Find the overlap between two pairs of coordinate variables and return the 
    indices to insert the overlapping region of the second pair into the first.
    If there is no overlap, None is returned instead of the array of indexes.
    NOTE: The indexes returned point to the element numbers of the coordinate
    vars which are the inclusive limits.  This is not the same as the python
    NUMPY indices needed to access the array slices.
    """

    if verbose:
        print " Mosaic-X  %8.3f - %8.3f  %d" % (mxvec[0], mxvec[-1], len(mxvec))
        print "   Test-X  %8.3f - %8.3f  %d" % (txvec[0], txvec[-1], len(txvec))
        print " Mosaic-Y  %8.3f - %8.3f  %d" % (myvec[0], myvec[-1], len(myvec))
        print "   Test-Y  %8.3f - %8.3f  %d" % (tyvec[0], tyvec[-1], len(tyvec))

    # Compute overlap region in test vectors.  If len(either) is zero, then
    # there is no overlap in that dimension and we are done
    # The logic here is as follows: 
    #    Form boolean array    True where txvec >= minumum mxvec
    #                   and    TRUE where txvec <= maximum mxvec
    # Multiply them (effectively logical-AND) to form boolean array True
    # where txvec is within range of mxvec.  
    # Use this to index txvec to produce reduced txvec containing only data in range
    tx = txvec[(txvec >= mxvec[0])*(txvec <= mxvec[-1])]
    if len(tx) == 0:
        if verbose: print "No X overlap"
        return None

    ty = tyvec[(tyvec <= myvec[0])*(tyvec >= myvec[-1])]
    if len(ty) == 0:
        if verbose: print "No Y overlap"
        return None

    # Proceed to extract index ranges by using reduced test vectors to derive
    # indexes for both test and mosaic vectors
    iv = np.arange(len(txvec))[(txvec >= tx[0])*(txvec <= tx[-1])]
    tx1 = iv[0]
    tx2 = iv[-1]
    iv = np.arange(len(mxvec))[(mxvec >= tx[0])*(mxvec <= tx[-1])]
    mx1 = iv[0]
    mx2 = iv[-1]

    iv = np.arange(len(tyvec))[(tyvec <= ty[0])*(tyvec >= ty[-1])]
    ty1 = iv[0]
    ty2 = iv[-1]
    iv = np.arange(len(myvec))[(myvec <= ty[0])*(myvec >= ty[-1])]
    my1 = iv[0]
    my2 = iv[-1]

    if verbose:
        print "Overlap X"
        print "   Mosaic %5d - %5d  (%d)  %8.3f  -  %8.3f  (%8.3f)" % (mx1, mx2, mx2-mx1, mxvec[mx1], mxvec[mx2], mxvec[mx2]-mxvec[mx1])
        print "     Test %5d - %5d  (%d)  %8.3f  -  %8.3f  (%8.3f)" % (tx1, tx2, tx2-tx1, txvec[tx1], txvec[tx2], txvec[tx2]-txvec[tx1])
        print "Overlap Y"
        print "   Mosaic %5d - %5d  (%d)  %8.3f  -  %8.3f  (%8.3f)" % (my1, my2, my2-my1, myvec[my1], myvec[my2], myvec[my1]-myvec[my2])
        print "     Test %5d - %5d  (%d)  %8.3f  -  %8.3f  (%8.3f)" % (ty1, ty2, ty2-ty1, tyvec[ty1], tyvec[ty2], tyvec[ty1]-tyvec[ty2])

    return [tx1, tx2, mx1, mx2, ty1, ty2, my1, my2]


# Support code for testing.


if __name__ == '__main__':

    def makeVec(start, stop, n):
        """
        Return an numpy f64 vector from START to STOP with N elements
        """
        return (np.arange(int(n))*(np.float64(stop)-np.float64(start))/(np.float64(n)-1.0)) + np.float64(start)


    # Set coordinate vars for mosaic window
    mlat = makeVec(10.0, -60.0, 7001)
    mlon = makeVec(80.0, 180.0, 10001)

    print "Region entirely enclosed"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region across RHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(170.0, 190.0, 2001)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region outside RHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(190.0, 215.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region on RHS boundary"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(180.0, 205.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""
    
    print "Region off by one RHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(180.01, 205.01, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""
    
    print "Region across LHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(70.0, 90.0, 2001)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region outside LHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(50.0, 70.0, 2001)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region on LHS boundary"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(60.0, 80.0, 2001)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""
    
    print "Region off by one LHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(59.99, 79.99, 2001)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""


    print "Region across TOP"
    tlat = makeVec(20.0, 0.0, 2001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region outside TOP"
    tlat = makeVec(40.0, 20.0, 2001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region on TOP boundary"
    tlat = makeVec(20.0, 10.0, 1001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region off by one TOP"
    tlat = makeVec(20.01, 10.01, 1001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""


    print "Region across BOTTOM"
    tlat = makeVec(-50.0, -70.0, 2001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region outside BOTTOM"
    tlat = makeVec(-65.0, -85.0, 2001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region on BOTTOM boundary"
    tlat = makeVec(-60.0, -70.0, 1001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region off by one BOTTOM"
    tlat = makeVec(-60.01, -70.01, 1001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""



    print "Region across TOP and BOTTOM"
    tlat = makeVec(20.0, -70.0, 9001)
    tlon = makeVec(110.0, 135.0, 2501)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

    print "Region across LHS and RHS"
    tlat = makeVec(-10.0, -20.0, 1001)
    tlon = makeVec(70.0, 190.0, 12001)
    oindicies = overlap(mlon, mlat, tlon, tlat, True)
    print ""

