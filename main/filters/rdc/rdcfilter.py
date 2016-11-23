import ConvertUTR
import math
import numpy as np
from scipy.optimize import leastsq

"""
Module to provide subroutines to parse and fit RDC data
"""


def ZYZRot(A, B, G, scal=1.0):
    """
    *** Adopted from PyParaTools ***
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
    from math import cos
    from math import sin
    from math import radians

    rot = np.zeros((3, 3))
    ca = cos(radians(A))
    cb = cos(radians(B))
    cg = cos(radians(G))
    sa = sin(radians(A))
    sb = sin(radians(B))
    sg = sin(radians(G))
    rot[0][0] = ((-sg * sa) + (cb * ca * cg)) * scal
    rot[0][1] = ((sg * ca) + (cb * sa * cg)) * scal
    rot[0][2] = ((-cg * sb)) * scal
    rot[1][0] = ((-cg * sa) - (cb * ca * sg)) * scal
    rot[1][1] = ((cg * ca) - (cb * sa * sg)) * scal
    rot[1][2] = ((sg * sb)) * scal
    rot[2][0] = ((sb * ca)) * scal
    rot[2][1] = ((sb * sa)) * scal
    rot[2][2] = (cb) * scal
    return rot


def getGMR(spin_type1, spin_type2):
    """
    *** Adopted from PyParaTools ***
    Return the gyromagnetic ratio(s) for the given spin type.
    From: http://nmrwiki.org/wiki/index.php?title=Gyromagnetic_ratio
    The values defined are consistent with those in Xplor-NIH.
    :param spin_type: The atom identifier for the spin (can be H, C, N)
    :type spin_type: String (H, C, N)
    :rtype: float
    """
    from math import pi
    mgr = {'H': ((2 * pi) * 42.576) * 1e6, \
           'C': ((2 * pi) * 10.705) * 1e6, \
           'CA': ((2 * pi) * 10.705) * 1e6, \
           'N': ((2 * pi) * -4.315) * 1e6,}

    return mgr[spin_type1] * mgr[spin_type2]


def rdcScal(S, B0, temp):
    """
    Scaling constant.for RDC calculations
    """
    # TODO: These need to be checked
    hbar = 1.05457148e-34
    kboltz = 1.3806503e-23
    scal = -S * hbar * B0 * B0 / (8 * 15 * math.pi * math.pi * kboltz * temp)
    return scal * 0.01


def getVector(c1, c2):
    """
    return the vector for given two coordinates
    """
    return [c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2]]


def RDC_ZYZ(p0, scal, vec_data):
    """
    The RDC function for ZYZ notation
    :param p0:                     Parameters Dax, Drh, A, B, G
    :vec_data:                     Vector def, GMR*GMR, RDC
    :param scal:                   The RDC scaling constants
    """

    n = len(vec_data)
    err = np.zeros(n)
    rot = ZYZRot(p0[2], p0[3], p0[4])

    for i in xrange(n):
        # [[-0.6490001678466797, -0.3989999294281006, 0.6469993591308594], -7252794860688272.0, -6.593]]
        #             X                  Y                    Z                 gm1*gm2            rdc
        Vec = vec_data[i][0]
        GMRp = vec_data[i][1]
        rdc_measured = vec_data[i][2]
        kscal = scal * GMRp
        X = Vec[0]
        Y = Vec[1]
        Z = Vec[2]
        x_t = rot[0][0] * X + rot[0][1] * Y + rot[0][2] * Z
        y_t = rot[1][0] * X + rot[1][1] * Y + rot[1][2] * Z
        z_t = rot[2][0] * X + rot[2][1] * Y + rot[2][2] * Z
        r2 = (x_t * x_t) + (y_t * y_t) + (z_t * z_t)
        r5 = (r2 * r2) * math.sqrt(r2)
        tmp = 1.0 / r5
        rdc = ((p0[0] * (3.0 * z_t * z_t - r2) + p0[1] * 1.5 * (x_t * x_t - y_t * y_t)) * tmp) * kscal
        err[i] = rdc_measured - rdc
    return err


def getVectorData(s1_def, s2_def, smotif, exp_data):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param exp_data:
    :return:
    """

    rdc_vector = []
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)
    sse_list = ss1_list + ss2_list

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)
    smotif_list = smotif_ss1 + smotif_ss2
    smotif_coors = smotif[1] + smotif[2]

    rdc_data = exp_data['rdc_data']

    for data in rdc_data:
        t_rdc_vector = []
        for res in sse_list:
            res_no = smotif_list[sse_list.index(res)]
            try:
                for entry in data[res_no]:
                    # [58, 'H', 58, 'N', 0.725]
                    for coor in smotif_coors:
                        # [31, 'PHE', 'H', 9.467, -2.327, 10.206]
                        if coor[0] == res_no and coor[2] == 'N':
                            atco1 = coor[3:]
                        if coor[0] == res_no and coor[2] == 'H':
                            atco2 = coor[3:]
                    # print res_no, atco1, atco2
                    t_rdc_vector.append([getVector(atco2, atco1), getGMR('H', 'N'), entry[4]])
            except:
                pass
                """
                res1 = sse_list[smotif_list.index(entry[0])]
                res2 = sse_list[smotif_list.index(entry[2])]
                for coor in smotif_coors:
                    #[31, 'PHE', 'H', 9.467, -2.327, 10.206]
                    if coor[0] == res1 and coor[2] == entry[1]:
                        atco1 = coor[3:]
                    if coor[0] == res2 and coor[2] == entry[1]:
                        atco2 = coor[3:]
                print res1, res2, atco1, atco2
                """
        rdc_vector.append(t_rdc_vector)
    return rdc_vector


def RDCAxRhFit(s1_def, s2_def, smotif, exp_data):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """
    rdc_vectors = getVectorData(s1_def, s2_def, smotif, exp_data)
    pred_axial = exp_data['pred_axial']
    exp_error = exp_data['exp_error']
    total_chisqr = 0
    all_tensors = []
    temp_tensor = []

    for i in range(0, len(rdc_vectors)):
        B0 = 18.8
        TinK = 298.0
        Sorder = 1.0
        scal = rdcScal(Sorder, B0, TinK)
        p0 = [pred_axial[i], pred_axial[i] * (0.33), 10, 20, 30]
        maxcs = 10000  # maximum number of optimization calls
        try:
            soln, cov, info, mesg, success = leastsq(RDC_ZYZ, p0, args=(scal, rdc_vectors[i]), full_output=1,
                                                     maxfev=maxcs)
            chisq = sum(info["fvec"] * info["fvec"]) / len(rdc_vectors[i])
            tensor = ConvertUTR.AnglesUTR(soln)

            if abs(tensor[0]) > pred_axial[i]:
                chisq = 999.999
            if chisq > exp_error[i]:
                chisq = 999.999
        except:
            chisq = 999.999

        if chisq < 999.999:
            temp_tensor.append([chisq, tensor])
    if len(temp_tensor) == len(rdc_vectors):
        return temp_tensor
    else:
        return []


def getNHvectors(coo_arrays, sse_list):
    """
    link to contacts_filter.py

    :param coo_arrays:
    :param sse_list:
    :return:
    """
    import filters.pcs.pcsfilter as pcs
    nh_dict = {}
    for i in range(0, len(sse_list)):
        sse_def = sse_list[i]
        h_array = pcs.getAtomCoo(coo_arrays[i], atom_type='H')
        n_array = pcs.getAtomCoo(coo_arrays[i], atom_type='N')

        sse_range = range(sse_def[4], sse_def[5] + 1)

        for j in range(0, len(sse_range)):
            try:
                nh_dict[sse_range[j]] = [getVector([h_array[0][j], h_array[1][j], h_array[2][j]],
                                                   [n_array[0][j], n_array[1][j], n_array[2][j]])]
            except:
                pass
    return nh_dict


def getVectorData2(transformed_coos, sse_ordered, rdc_data):
    rdc_vector = []
    nh_vectors = getNHvectors(transformed_coos, sse_ordered)
    nh_residues = nh_vectors.keys()
    for data in rdc_data:
        t_rdc_vector = []
        residues = data.keys()
        for res in residues:
            if res in nh_residues:
                t_rdc_vector.append([nh_vectors[res][0], getGMR('H', 'N'), data[res][0][-1]])
        rdc_vector.append(t_rdc_vector)

    return rdc_vector


def RDCAxRhFit2(transformed_coos, sse_ordered, exp_data, stage):
    # print transformed_coos
    # nh_dict = pcs.coorNHdict(transformed_coos, sse_ordered)
    rdc_vectors = getVectorData2(transformed_coos, sse_ordered, exp_data['rdc_data'])
    pred_axial = exp_data['pred_axial']
    exp_error = exp_data['exp_error']
    temp_tensor = []
    for i in range(0, len(rdc_vectors)):
        B0 = 18.8
        TinK = 298.0
        Sorder = 1.0
        scal = rdcScal(Sorder, B0, TinK)
        p0 = [pred_axial[i], pred_axial[i] * (0.33), 10, 20, 30]
        maxcs = 10000  # maximum number of optimization calls
        try:
            soln, cov, info, mesg, success = leastsq(RDC_ZYZ, p0, args=(scal, rdc_vectors[i]), full_output=1,
                                                     maxfev=maxcs)
            chisq = sum(info["fvec"] * info["fvec"]) / len(rdc_vectors[i])
            tensor = ConvertUTR.AnglesUTR(soln)
            if abs(tensor[0]) > pred_axial[i]:
                chisq = 999.999

            if chisq > exp_error[i]:
                chisq = 999.999
            """
            if exp_error[i] - 2 <= chisq <= exp_error[i] + 2:
                pass
            else:
                chisq = 999.999
            """
        except:
            chisq = 999.999

        if chisq < 999.999:
            temp_tensor.append([chisq, tensor])
    if len(temp_tensor) == len(rdc_vectors):
        return temp_tensor
    else:
        return []