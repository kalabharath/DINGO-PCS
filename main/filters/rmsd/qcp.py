#!/usr/bin/env python

"""
Project_Name: rmsd/qcp, File_name: qcp.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 8/05/15 , Time:4:27 AM
"""

import copy

import qcprot


def dumpPDBCoo(coo_array):
    """

    :param coo_array:
    :return:
    """
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]

        pdb_line = "%-6s%5d %2s%1s%s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM', i+1, 'CA', "", 'ALA', 'A', i+1, " ", x, y, z, 1.0, 30.0, ' ', ' ')
        print pdb_line
    print 'TER'
    return True


def dumpPDBCoo2(coo_array):
    """
 
    :param coo_array:
    :return:
    """
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]

        pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM', i+1, atom, res, 'A', res_no, " ", x, y, z, 1.0, 30.0, ' ', ' ')
        print pdb_line
    print 'TER'
    return True


def getCAcoo(frag):
    """
    :param frag:
    :return:
    """
    #print frag
    x, y, z = [], [], []
    for i in range(0, len(frag[0])):
        if frag[3][i] == 'CA':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
    return [x, y, z]

def getCcoo(frag):
    """
    :param frag:
    :return:
    """
    #print frag
    x, y, z = [], [], []
    for i in range(0, len(frag[0])):
        if frag[3][i] == 'C':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
    return [x, y, z]


def getcoo(frag):
    """

    :param frag:
    :return:
    """
    #[41, 'ASP', 'N', 28.117, -17.694, 7.215]
    x, y, z, atom_type, res_no, res = [], [], [], [], [], []
    for i in range(0, len(frag)):
        x.append(frag[i][3])
        y.append(frag[i][4])
        z.append(frag[i][5])
        atom_type.append(frag[i][2])
        res_no.append(frag[i][0])
        res.append(frag[i][1])
    return [x, y, z, atom_type, res_no, res]


def centerCoo(coo_array):
    """

    :param coo_array:
    :return:
    """

    xsum, ysum, zsum = 0, 0, 0

    for i in range(0, len(coo_array[0])):
        xsum += coo_array[0][i]
        ysum += coo_array[1][i]
        zsum += coo_array[2][i]

    xsum /= len(coo_array[0])
    ysum /= len(coo_array[0])
    zsum /= len(coo_array[0])

    for i in range(0, len(coo_array[0])):
        coo_array[0][i] -= xsum
        coo_array[1][i] -= ysum
        coo_array[2][i] -= zsum

    return coo_array, [xsum, ysum, zsum]


def applyTranslation(frag, cen_mass):
    """

    :param frag:
    :param cen_mass:
    :return:
    """

    for i in range(0, len(frag[0])):
        frag[0][i] += cen_mass[0]
        frag[1][i] += cen_mass[1]
        frag[2][i] += cen_mass[2]
    return frag


def translateCM(coo_array, cm):
    """

    :param coo_array:
    :param cm:
    :return:
    """
    xsum, ysum, zsum = cm[0], cm[1], cm[2]
    for i in range(0, len(coo_array[0])):
        coo_array[0][i] -= xsum
        coo_array[1][i] -= ysum
        coo_array[2][i] -= zsum
    return coo_array


def applyRot(frag, rotmat):
    """

    :param frag:
    :param rotmat:
    :return:
    """
    for i in range(0, len(frag[0])):
        x = rotmat[0]*frag[0][i] + rotmat[1]*frag[1][i] + rotmat[2]*frag[2][i]
        y = rotmat[3]*frag[0][i] + rotmat[4]*frag[1][i] + rotmat[5]*frag[2][i]
        z = rotmat[6]*frag[0][i] + rotmat[7]*frag[1][i] + rotmat[8]*frag[2][i]
        frag[0][i] = x
        frag[1][i] = y
        frag[2][i] = z
    return frag


def get_dist(r1, r2):
    import math
    x1, y1, z1 = r1[0], r1[1], r1[2]
    x2, y2, z2 = r2[0], r2[1], r2[2]
    return math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))


def rmsdQCP(psmotif, csmotif, direction):
    """
    use rmsdQCP
    :param psmotif:
    :param csmotif:
    :param direction:
    :return:
    """

    psmotif = copy.copy(psmotif[1])

    if direction == 'left':
        native_fraga = getcoo(psmotif[1])
        frag_a = getcoo(psmotif[1])
        frag_b = getcoo(csmotif[2])
        native_fragb_2ndsse = copy.copy(csmotif[1])
        native_fraga_2ndsse = getcoo(psmotif[2])
    else:

        native_fraga = getcoo(psmotif[2])
        frag_a = getcoo(psmotif[2])
        frag_b = getcoo(csmotif[1])
        native_fragb_2ndsse = copy.copy(csmotif[2])
        native_fraga_2ndsse = getcoo(psmotif[1])

    frag_a, a_cen = centerCoo(frag_a)
    frag_b, b_cen = centerCoo(frag_b)

    frag_aca = getCAcoo(frag_a)
    frag_bca = getCAcoo(frag_b)

    fraglen = len(frag_aca[0])

    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_aca[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_aca[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_aca[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_bca[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_bca[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_bca[2][i])

    rot = qcprot.MakeDvector(9)

    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********

    rotmat = []
    for i in range(0, 9):
        rotmat.append(qcprot.GetDvector(i, rot))

    rotated_fragb = applyRot(frag_b, rotmat)

    # trans_fragb = applyTranslation(rotated_fragb, a_cen)

    # translate the other SSE of the current smotif
    sse_2nd_coos = getcoo(native_fragb_2ndsse)
    # sse_2nd_coos, sec_cm = centerCoo(sse_2nd_coos)
    cm_sse2nd = translateCM(sse_2nd_coos,b_cen)
    rot_sse_2nd = applyRot(cm_sse2nd, rotmat)
    trans_sse2nd = applyTranslation(rot_sse_2nd, a_cen)

    #return 3 arrays of coordinates
    if direction == 'left':
        #transformed_coor = [native_fraga, native_fraga_2ndsse, trans_sse2nd]

        transformed_coor = [trans_sse2nd, native_fraga, native_fraga_2ndsse]
        #
    else:

        transformed_coor = [native_fraga_2ndsse, native_fraga, trans_sse2nd]

    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd, transformed_coor


def rmsdQCP3(previous_smotif, csmotif, direction):
    """

    :param previous_smotif:
    :param csmotif:
    :param direction:
    :return:
    """

    for entry in previous_smotif:
        if 'qcp_rmsd' == entry[0]:
            temp_holder = entry[1]
            presse = temp_holder[:]

    #print csmotif
    if direction == 'left':
        frag_b = getcoo(csmotif[2])
        native_fragb_2ndsse = copy.copy(csmotif[1])
        frag_a = copy.deepcopy(presse[0])

    else:
        try:
            frag_a = copy.deepcopy(presse[-1])
            frag_b = getcoo(csmotif[1])
            native_fragb_2ndsse = copy.copy(csmotif[2])
        except:
            print previous_smotif

    frag_a, a_cen = centerCoo(frag_a)
    frag_b, b_cen = centerCoo(frag_b)

    frag_aca = getCAcoo(frag_a)
    frag_bca = getCAcoo(frag_b)

    fraglen = len(frag_aca[0])

    xyz1 = qcprot.MakeDMatrix(3, fraglen)
    xyz2 = qcprot.MakeDMatrix(3, fraglen)

    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz1, frag_aca[0][i])
        qcprot.SetDArray(1, i, xyz1, frag_aca[1][i])
        qcprot.SetDArray(2, i, xyz1, frag_aca[2][i])
    for i in range(0, fraglen):
        qcprot.SetDArray(0, i, xyz2, frag_bca[0][i])
        qcprot.SetDArray(1, i, xyz2, frag_bca[1][i])
        qcprot.SetDArray(2, i, xyz2, frag_bca[2][i])

    rot = qcprot.MakeDvector(9)

    # *********
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)
    # *********

    # print rmsd

    rotmat = []
    for i in range(0, 9):
        rotmat.append(qcprot.GetDvector(i, rot))

    # translate the other SSE of the current smotif

    sse_2nd_coos = getcoo(native_fragb_2ndsse)
    cm_sse2nd = translateCM(sse_2nd_coos, b_cen)
    rot_sse_2nd = applyRot(cm_sse2nd, rotmat)
    trans_sse2nd = applyTranslation(rot_sse_2nd, a_cen)

    # append the translated coordinates
    temp_holder = copy.deepcopy(presse)

    if direction == 'left':
        temp_holder.insert(0, trans_sse2nd)
    else:
        temp_holder.append(trans_sse2nd)


    qcprot.FreeDMatrix(xyz1)
    qcprot.FreeDMatrix(xyz2)
    qcprot.FreeDArray(rot)

    return rmsd, temp_holder


def clahses(coo_arrays, cdist):
    """
    Find the clashes between all the SSEs
    :param coo_arrays:
    :return:
    """

    for i in range(0, len(coo_arrays)-1):
        sse1 = getCAcoo(coo_arrays[i])
        for p in range(i+1,len(coo_arrays)):
            sse2 = getCAcoo(coo_arrays[p])
            for j in range(0, len(sse1[0])):
                for k in range(0, len(sse2[0])):
                    dist = get_dist([sse1[0][j], sse1[1][j], sse1[2][j]], [sse2[0][k], sse2[1][k], sse2[2][k]])
                    if dist < cdist:
                        return False

    for i in range(0, len(coo_arrays)-1):
        sse1 = getCcoo(coo_arrays[i])
        for p in range(i+1,len(coo_arrays)):
            sse2 = getCcoo(coo_arrays[p])
            for j in range(0, len(sse1[0])):
                for k in range(0, len(sse2[0])):
                    dist = get_dist([sse1[0][j], sse1[1][j], sse1[2][j]], [sse2[0][k], sse2[1][k], sse2[2][k]])
                    if dist < cdist:
                        return False

    return True
