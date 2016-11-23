import sys
from math     import degrees
from numpy    import zeros, ones, sqrt, array, mat, append, empty,concatenate
# Forked from PyParaTools

def ZYZRot(A, B, G, scal=1.0):
    """
    Returns a ZYZ rotation matrix when given 3 Euler angles (in degrees).

    See: http://mathworld.wolfram.com/EulerAngles.html

    :param A: The (A)lpha angle (degrees)
    :type  A: float
    :param B: The (B)eta angle (degrees)
    :type  B: float
    :param G: The (G)amma angle (degrees)
    :type  G: float
    :param scal: Scale the rotation matix by a constant [OPTIONAL]
    :type  scal: float

    :rtype: Numpy matrix
    """
    from numpy import zeros
    from math  import cos
    from math  import sin
    from math  import radians

    rot = zeros((3, 3))
    ca = cos(radians(A))
    cb = cos(radians(B))
    cg = cos(radians(G))
    sa = sin(radians(A))
    sb = sin(radians(B))
    sg = sin(radians(G))
    rot[0][0] = ((-sg * sa) + (cb * ca * cg))*scal
    rot[0][1] = (( sg * ca) + (cb * sa * cg))*scal
    rot[0][2] = (( -cg * sb))*scal
    rot[1][0] = ((-cg * sa) - (cb * ca * sg))*scal
    rot[1][1] = ((cg * ca) - (cb * sa * sg))*scal
    rot[1][2] = (( sg * sb))*scal
    rot[2][0] = (( sb * ca))*scal
    rot[2][1] = ((sb * sa))*scal
    rot[2][2] =  (cb)*scal
    return rot

def RotX90():
    """
    Returns the rotation matrix for a 90 degree rotation about X.

    Defined as: [[1,  0, 0], [0,  0, 1], [0, -1, 0]]. Used in UTR determination

    :rtype: Numpy matrix
    """
    from numpy import zeros

    rot = zeros((3, 3))
    rot[0][0] =  1.0
    rot[1][2] =  1.0
    rot[2][1] = -1.0
    return rot

def RotY90():
    """Returns the rotation matrix for a 90 degree rotation about Y.

    Defined as: [[0, 0, -1], [0, 1, 0], [1, 0, 0]]. Used in UTR determination

    :rtype: Numpy matrix
    """
    from numpy import zeros

    rot = zeros((3, 3))
    rot[0][2] = -1.0
    rot[1][1] =  1.0
    rot[2][0] =  1.0
    return rot


def correctAngle(cosv, sinv):
    """Return an angle in correct quadrant given a cosine and sine of the angle.

    Applies sin in [-pi/2, x, pi/2] and cos in [0, x, pi]

    This agrees with Numbat code.

    :param cosv: The value for cos (radians)
    :type  cosv: float
    :param sinv: The value for sin (radians)
    :type  sinv: float

    :rtype: float - the angle in the correct quadrant in radians
    """
    from math import pi

    if (cosv <= pi/2.0):
        if (sinv < 0.0):
            sinv += 2 * pi
            return sinv
        else:
            return sinv
    else:
        if(sinv > 0.0):
            return cosv
        else:
            return -1*(cosv) +2*pi


def ABGFromRotMatrixZYZ(rotMat):
    """
    Return the 3 Euler angles (A, B, G) in ZYZ from a given rotation matrix.

    :param rotMat: The ZYZ rotation matrix
    :type  rotMat: 3x3 Numpy Matrix

    .. warning:: At the present moment this method does not handle cases where
        the generated rotation matrix contains Euler angles at 0,90,180,270 or
        360 degrees respectively. Need to look how Numbat handles this.

    :rtype: tuple - Alpha, Beta, Gamma
    """
    #Rewrite so to remove exception (see docstring)
    from math  import acos
    from math  import asin
    from math  import sin
    from sys   import exit

    try:
        b_c = acos(rotMat[2, 2])
        a_c = acos(rotMat[2, 0]/sin(b_c))
        g_c = acos(-1*rotMat[0, 2]/sin(b_c))
        a_s = asin(rotMat[2, 1]/sin(b_c))
        g_s = asin(rotMat[1, 2]/sin(b_c))
        aE = correctAngle(a_c, a_s)
        bE = b_c
        gE = correctAngle(g_c, g_s)
        return aE, bE, gE
    except ValueError:
        print 80*'-'
        print "ERROR: One of A,B or G is equal to 0, 90, 180, 270 or 360 deg"
        print "This is a known issue and will be fixed in future versions"
        print 80*'-'
        exit(0)


def FromVVU(AxorRh):
    """
    Return the Axial or Rhombic component converted from van Vleck units.

    van Vleck units = m3/3.77 10-35, See:
    http://www.cerm.unifi.it/Downloads/Xplor/PARArestraints.pdf

    :param AxorRh: Axial or Rhombic component of the X-tensor
    :type AxorRh : float

    :rtype: float
    """
    from math import pi

    return AxorRh/(1./((12*pi))*10000.0)

def ToVVU(AxorRh):
    """
    Return the Axial or Rhombic component converted to van Vleck units.

    van Vleck units = m3/3.77 10-35, See:
    http://www.cerm.unifi.it/Downloads/Xplor/PARArestraints.pdf

    :param AxorRh: Axial or Rhombic component of the X-tensor
    :type AxorRh : float

    :rtype: float
    """
    from math import pi

    return AxorRh*(1./((12*pi))*10000.0)


def FixAngle(angle):
    """Return an angle such that it is [0:2pi] bound.

    Useful as the angles determined from an optimization are not [0:2pi] bound

    :param angle: An Euler angle determined from the optimization
    :type angle:  float

    :rtype: float (in degrees)
    """
    return angle % 360.0


def AnglesUTR(p0, ref=0, verbose=False):
                
    """
    Reconfigure a X-tensor or Alignment-tensor into Unique Tensor Rep
    
    Core code has been borrowed from Numbat. Information on UTR can be
    found in JBNMR. 2008 Jul;41(3):179-89. Thanks Christophe Schmitz :-)

    :param ref:     The reference of ParsePara object associated with this
                    fitter. This defaults to 0 (the first ParsePara object)
                        [OPTIONAL]
    :type ref:      int
    :param verbose: Toggle between printing information on the UTR case (0)
                        or being quiet (non-zero). Default is non-zero (not
                        printing)
                        [OPTIONAL]
    :type  verbose: int (0 or non-zero)
    :rtype: a list containing Ax, Rh, A, B and G in UTR.
    """
    # p0=[ax,rh, a,b,g]
    # 0   1  2 3 4
    
    
    a = p0[2]
    b = p0[3]
    g = p0[4]       

    Dx = -ToVVU(p0[0])/3.0 + ToVVU(p0[1])/2.0
    Dy = -ToVVU(p0[0])/3.0 - ToVVU(p0[1])/2.0
    Dz = 2.0/3.0*(ToVVU(p0[0]))
    aDx, aDy, aDz = abs(Dx), abs(Dy), abs(Dz)

    # Determine the UTR case
    if (aDz >= aDy) and (aDy >= aDx):
        if verbose:
            print "UTR Case1"
    if (aDz >= aDx)and (aDx >= aDy):
        g += 90.0
        Dy, Dx = Dx, Dy
        if verbose:
            print "UTR Case2"
    if (aDy >= aDz) and (aDz >= aDx):
        Dy, Dz = Dz, Dy
        rX90 = RotX90()
        rZYZ = ZYZRot(a, b, g)
        nR = mat(rX90) * mat(rZYZ)
        a, b, g = ABGFromRotMatrixZYZ(nR)
        a, b, g = degrees(a), degrees(b), degrees(g)
        if verbose:
            print "UTR Case3"
    if (aDy >= aDx) and (aDx >= aDz):
        g += 90.0
        Dy, Dx = Dx, Dy
        Dz, Dx = Dx, Dz
        rY90 = RotY90()
        rZYZ = ZYZRot(a, b, g)
        nR = mat(rY90) * mat(rZYZ)
        a, b, g = ABGFromRotMatrixZYZ(nR)
        a, b, g = degrees(a), degrees(b), degrees(g)
        if verbose:
            print "UTR Case4"
    if(aDx >= aDz) and (aDz >= aDy):
        g += 90.0
        Dy, Dx = Dx, Dy
        Dy, Dz = Dz, Dy
        rX90 = RotX90()
        rZYZ = ZYZRot(a, b, g)
        nR = mat(rX90) * mat(rZYZ)
        a, b, g = ABGFromRotMatrixZYZ(nR)
        a, b, g = degrees(a), degrees(b), degrees(g)
        if verbose:
            print "UTR Case5"
    if(aDx >= aDy) and (aDy >= aDz):
        Dz, Dx = Dx, Dz
        rY90 = RotY90()
        rZYZ = ZYZRot(a, b, g)
        nR =  mat(rY90)* mat(rZYZ)
        a, b, g = ABGFromRotMatrixZYZ(nR)
        a, b, g = degrees(a), degrees(b), degrees(g)
        if verbose:
            print "UTR Case6"

    #Axial and Rhombic are now in UTR
    Ax = Dz - (Dx + Dy)/2.0
    Rh = Dx - Dy
    Ax, Rh = FromVVU(Ax), FromVVU(Rh)

    # Make Euler angles in 0-360 after manipulation.
    a = FixAngle(a)
    b = FixAngle(b)
    g = FixAngle(g)

    # Do manipulations such that A,B,G in 0-180
    if a >= 0.0 and a < 180.0:
        if b >= 0.0  and  b < 180.0:
            if g >= 0.0  and  g < 180.0:
                pass
            else:
                g += 180.0
        else:
            if g >= 0.0 and g < 180.0:
                b += 180.0
                g = -g +180
            else:
                b += 180.0
                g = -g
    else:
        if b >= 0 and  b < 180.0:
            if g >= 0  and  g < 180.0:
                a += 180.0
                b = -b + 180.0
                g = -g + 180.0
            else:
                a += 180.0
                b = -b + 180.0
                g = -g
        else:
            if g >= 0 and  g < 180.0:
                a += 180.0
                b = -b
                g =  g
            else:
                a += 180.0
                b = -b
                g += 180.0

    # Important. Fix to 0-360 to get in UTR (really 0-180).
    a = FixAngle(a)
    b = FixAngle(b)
    g = FixAngle(g)

    #Update for UTR!
    return [Ax, Rh, a, b, g]
