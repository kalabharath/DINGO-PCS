#!/usr/bin/env python

"""
Project_Name: main/filters/pcsfilter, File_name: pcsfilter.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 24/04/15 , Time: 01:30 PM

fit PCS and apply Ax,Rh filters
"""

import math

import fastT1FM


def getHN(ss1_list, ss2_list, smotif, atom_type='H'):
    """
    parse the x,y,z of individual SSE into seperate arrays
    :param smotif:
    :return: arrays of coordinates of individual SSE in seperate arrays
    """
    # TODO return arrays for any given atom type(s)

    rH1, rH2 = [], []
    counter = 0
    for entry in smotif[1]:
        # res_no, aa_type, atom_type, x, y, z
        if entry[2] == atom_type:
            x, y, z = entry[3], entry[4], entry[5]
            res_no = ss1_list[counter]
            counter += 1
            if counter == len(ss1_list):
                break
            rH1.append([x, y, z, res_no])
    counter = 0
    for entry in smotif[2]:
        # res_no, aa_type, atom_type, x, y, z
        if entry[2] == atom_type:
            x, y, z = entry[3], entry[4], entry[5]
            res_no = ss2_list[counter]
            counter += 1
            if counter == len(ss2_list):
                break
            rH1.append([x, y, z, res_no])
    return rH1, rH2


def match_pcss_HN(rh1, rh2, pcs_data):
    """
    :rtype : object
    :param rh1:
    :param rh2:
    :param pcs_data:
    :return:
    """
    smotif_pcs = []

    for entry in rh1:
        res_no = (entry[-1]) - 1  # because python numbering starts from '0'
        smotif_pcs.append(pcs_data[res_no])
    for entry in rh2:
        res_no = (entry[-1]) - 1
        smotif_pcs.append(pcs_data[res_no])
    return smotif_pcs


def PointsOnSpheres(M, N, rMx, rMy, rMz):
    """
    quick way from wikipedia
    :param M:
    :param N:
    :param rMx:
    :param rMy:
    :param rMz:
    :return:
    """

    import math

    node = []
    dlong = math.pi * (3 - math.sqrt(5))  # ~2.39996323
    dz = 2.0 / N
    xlong = 0
    z = 1 - 0.5 * dz
    for k in range(N):
        r = math.sqrt(1 - z * z)
        node.append([math.cos(xlong) * r, math.sin(xlong) * r, z])
        z -= dz
        xlong += dlong

    j = 0
    for i in range(M[0], M[1]):
        for k in range(N):
            fastT1FM.SetDvector(j, rMx, i * node[k][0])
            fastT1FM.SetDvector(j, rMy, i * node[k][1])
            fastT1FM.SetDvector(j, rMz, i * node[k][2])
            j += 1


def fixedPointsOnSpheres(M, sN, rMx, rMy, rMz):
    """
    quick way from wikipedia
    :param M:
    :param sN: num of points on the starting sphere
    :param rMx:
    :param rMy:
    :param rMz:
    :return:
    Hardcoded increments from 200 pts innershell * number of shells
    """

    # import math
    j = 0
    for i in range(M[0], M[1]):
        N = i * sN
        node = []
        # dlong = math.pi * (3 - math.sqrt(5))  # ~2.39996323
        dlong = 2.39996322973
        dz = 2.0 / N
        xlong = 0
        z = 1 - 0.5 * dz

        for l in range(N):
            r = math.sqrt(1 - z * z)
            node.append([math.cos(xlong) * r, math.sin(xlong) * r, z])
            z -= dz
            xlong += dlong

        for k in range(N):
            fastT1FM.SetDvector(j, rMx, i * node[k][0])
            fastT1FM.SetDvector(j, rMy, i * node[k][1])
            fastT1FM.SetDvector(j, rMz, i * node[k][2])
            j += 1
    return True


def usuablePCS(pcs_array):
    """

    :param pcs_array:
    :return:
    """

    if not pcs_array:
        return 0, False
    else:
        total_pcs = 0
        for j in range(0, len(pcs_array[0])):
            counter = 0
            for entry in pcs_array:
                if entry[j] != 999.999:
                    counter += 1
            total_pcs = total_pcs + counter
            if counter <= 5:
                return total_pcs, False

    return float(total_pcs), True


def calcAxRh(saupe_matrices):
    """
    :param saupe_matrices:
    :return:
    """
    import numpy
    from numpy import linalg as LA

    axrh = []
    for t in saupe_matrices:
        w, v = LA.eig(numpy.array([[t[0], t[1], t[2]], [t[1], t[3], t[4]], [t[2], t[4], -t[0] - t[3]]]))
        x = []
        for i in range(3):
            x.append([abs(w[i]), w[i]])
        x.sort()
        for i in range(3):
            w[i] = x[i][1]
        axrh.append([w[2] - 0.5 * (w[0] + w[1]), w[0] - w[1]])
    return axrh


def checkAxRh(axrh, chisqr, total_pcs, stage):
    """

    :param axrh:
    :param chisqr:
    :param total_pcs:
    :param stage:
    :return:
    """
    if stage <=2:

        if (chisqr / total_pcs) > 0.0025:  # no error limit
            return 1.0e+30

    if stage == 3:

        if (chisqr / total_pcs) > (0.01):  # 4 times the standard error limit
            return 1.0e+30
    if stage == 4:
        if (chisqr / total_pcs) > (0.015 ): #6 times the standard limit
            return 1.0e+30

    for metal in axrh:
        for parameter in metal:
            if stage == 1 or stage == 2:
                if abs(parameter) > 80:
                    return 1.0e+30
            if stage == 3 or stage == 4:
                if abs(parameter) > 40:
                    return 1.0e+30

    return chisqr


def checkAxRhCutoffs(axrh, chisqr, total_pcs, exp_data, stage):
    """

    :param axrh:
    :param chisqr:
    :param total_pcs:
    :param stage:
    :param exp_data:
    :return:
    """
    # get the user specified cutoffs for chisqr, axial and rhombic parameters
    chisqr_cutoffs = exp_data['chisqr_cutoff']
    axrh_cutoffs = exp_data['axrh_cutoff']

    # Apply chisqr cutoffs
    if stage == 1:
        if (chisqr / total_pcs) > chisqr_cutoffs[0]:  # no error limit
            return 1.0e+30
    elif stage == 2:
        if (chisqr / total_pcs) > chisqr_cutoffs[1]:  # no error limit
            return 1.0e+30
    elif stage == 3:
        if (chisqr / total_pcs) > chisqr_cutoffs[2]:  # no error limit
            return 1.0e+30
    else:
        if (chisqr / total_pcs) > chisqr_cutoffs[3]:  # no error limit
            return 1.0e+30

    # Apply Axial and Rhombic cutoffs
    for metal in axrh:
        for parameter in metal:
            if stage == 1 :
                if abs(parameter) > axrh_cutoffs[0]:
                    return 1.0e+30
            elif stage == 2:
                if abs(parameter) > axrh_cutoffs[1]:
                    return 1.0e+30
            elif stage == 3:
                if abs(parameter) > axrh_cutoffs[2]:
                    return 1.0e+30
            else:
                if abs(parameter) > axrh_cutoffs[3]:
                    return 1.0e+30

    return chisqr


def PCSAxRhFit(s1_def, s2_def, smotif, exp_data):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """

    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)


    rH1, rH2 = getHN(ss1_list, ss2_list, smotif, atom_type='H')
    pcs_data = exp_data['pcs_data']
    ntags = len(pcs_data)

    # Define Thomas's implementaion of hollow concentric shells

    nM = 200  # 1000 pts in each sphere
    M = [1, 45]  # 40 spheres 10-50 Angstrom

    npts = 198000
    rMx = fastT1FM.MakeDvector(npts)  # allocate memmory
    rMy = fastT1FM.MakeDvector(npts)
    rMz = fastT1FM.MakeDvector(npts)
    # PointsOnSpheres(M, nM, rMx, rMy, rMz)
    fixedPointsOnSpheres(M, nM, rMx, rMy, rMz)

    # Temp storage of tensor values
    temp_tensor = []

    for tag in range(0, ntags):

        smotif_pcs = match_pcss_HN(rH1, rH2, pcs_data[tag])

        total_pcs, pcs_bool = usuablePCS(smotif_pcs)

        if pcs_bool:  # save some time for not running

            # Thomas's fast Tensor calc code init

            frag_len = len(smotif_pcs)
            nsets = len(smotif_pcs[0])
            xyz = fastT1FM.MakeDMatrix(frag_len, 3)
            pcs = fastT1FM.MakeDMatrix(nsets, frag_len)
            xyz_HN = rH1 + rH2

            for k in range(nsets):
                for j in range(frag_len):
                    fastT1FM.SetDArray(k, j, pcs, smotif_pcs[j][k])

            cm = [0.0, 0.0, 0.0]
            for j in range(frag_len):
                cm[0] = cm[0] + xyz_HN[j][0]
                cm[1] = cm[1] + xyz_HN[j][1]
                cm[2] = cm[2] + xyz_HN[j][2]
            cm[0] /= float(frag_len)
            cm[1] /= float(frag_len)
            cm[2] /= float(frag_len)
            for j in range(frag_len):
                fastT1FM.SetDArray(j, 0, xyz, xyz_HN[j][0] - cm[0])
                fastT1FM.SetDArray(j, 1, xyz, xyz_HN[j][1] - cm[1])
                fastT1FM.SetDArray(j, 2, xyz, xyz_HN[j][2] - cm[2])

            tensor = fastT1FM.MakeDMatrix(nsets, 8)
            Xaxrh_range = fastT1FM.MakeDMatrix(nsets, 4)
            for i in range(0, nsets):
                fastT1FM.SetDArray(i, 0, Xaxrh_range, 0.05)
                fastT1FM.SetDArray(i, 1, Xaxrh_range, 200.0)
                fastT1FM.SetDArray(i, 2, Xaxrh_range, 0.05)
                fastT1FM.SetDArray(i, 3, Xaxrh_range, 200.0)
            # ****
            chisqr = fastT1FM.rfastT1FM_multi(npts, rMx, rMy, rMz, nsets, frag_len, xyz, pcs, tensor, Xaxrh_range)
            # ****

            saupe_array = []
            for kk in range(nsets):
                temp_saupe = []
                for j in range(3, 8):
                    temp_saupe.append(fastT1FM.GetDArray(kk, j, tensor))
                saupe_array.append(temp_saupe)

            x = fastT1FM.GetDArray(0, 0, tensor)
            y = fastT1FM.GetDArray(0, 1, tensor)
            z = fastT1FM.GetDArray(0, 2, tensor)
            metal_pos = [x + cm[0], y + cm[1], z + cm[2]]


            # Compute and check Axial and Rhombic parameters
            AxRh = calcAxRh(saupe_array)

            #chisqr = checkAxRh(AxRh, chisqr, total_pcs, stage=1)  # modifies the values of chisqr
            chisqr = checkAxRhCutoffs(AxRh, chisqr, total_pcs, exp_data, stage=1)  # modifies the values of chisqr
            AxRh.append(metal_pos)  # add metal pos

            # Free memory for the variables
            fastT1FM.FreeDMatrix(xyz)
            fastT1FM.FreeDMatrix(pcs)
            fastT1FM.FreeDMatrix(Xaxrh_range)
            fastT1FM.FreeDMatrix(tensor)
        else:
            chisqr = 1.0e+30

        if chisqr < 1.0e+30:

            nchisqr = chisqr / float(total_pcs - (nsets * 5))
            snchisqr = nchisqr / float(math.pow(total_pcs, 1 / 3.0))
            temp_tensor.append([tag, snchisqr, AxRh])

    fastT1FM.FreeDArray(rMx)
    fastT1FM.FreeDArray(rMy)
    fastT1FM.FreeDArray(rMz)

    return temp_tensor


###For stage 2

def getAtomCoo(coo_array, atom_type):
    """
    link to contacts_filter.py

    :param coo_array:
    :return:
    """

    xt, yt, zt = [], [], []
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]
        if atom == atom_type:
            xt.append(x)
            yt.append(y)
            zt.append(z)
    return [xt, yt, zt]


def coorNHdict(coo_arrays, sse_list):
    """
    link to contacts_filter.py

    :param coo_arrays:
    :param sse_list:
    :return:
    """
    nh_dict = {}
    for i in range(0, len(sse_list)):
        sse_def = sse_list[i]
        nh_array = getAtomCoo(coo_arrays[i], atom_type='H')
        sse_range = range(sse_def[4], sse_def[5] + 1)

        for j in range(0, len(sse_range)):
            try:
                nh_dict[sse_range[j]] = [nh_array[0][j], nh_array[1][j], nh_array[2][j]]
            except:
                pass
    return nh_dict


def matchPCS(nh_dict, pcs_data):
    smotif_pcs = []
    xyz_HN = []

    residues = nh_dict.keys()
    residues.sort()
    for res in residues:
        smotif_pcs.append(pcs_data[res - 1])
        xyz = nh_dict[res]
        xyz.append(res)
        xyz_HN.append(xyz)

    return xyz_HN, smotif_pcs


def PCSAxRhFit2(transformed_coos, sse_ordered, exp_data, stage):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """

    nh_dict = coorNHdict(transformed_coos, sse_ordered)

    pcs_data = exp_data['pcs_data']

    # print "pcs_data", pcs_data
    ntags = len(pcs_data)

    # Define Thomas's implementaion of hollow concentric shells

    nM = 200  # 200 pts in starting sphere
    M = [1, 45]  # 40 spheres 10-50 Angstrom

    # npts = (M[1] - M[0]) * nM  # 50 spheres * 1000 pts each
    npts = 198000
    rMx = fastT1FM.MakeDvector(npts)  # allocate memmory
    rMy = fastT1FM.MakeDvector(npts)
    rMz = fastT1FM.MakeDvector(npts)
    # PointsOnSpheres(M, nM, rMx, rMy, rMz)
    fixedPointsOnSpheres(M, nM, rMx, rMy, rMz)

    # Temp storage of tensor values
    temp_tensor = []

    for tag in range(0, ntags):

        xyz_HN, smotif_pcs = matchPCS(nh_dict, pcs_data[tag])

        total_pcs, pcs_bool = usuablePCS(smotif_pcs)

        if pcs_bool:  # save some time for not running

            # Thomas's fast Tensor calc code
            frag_len = len(smotif_pcs)
            nsets = len(smotif_pcs[0])
            xyz = fastT1FM.MakeDMatrix(frag_len, 3)
            pcs = fastT1FM.MakeDMatrix(nsets, frag_len)

            for k in range(nsets):
                for j in range(frag_len):
                    fastT1FM.SetDArray(k, j, pcs, smotif_pcs[j][k])

            cm = [0.0, 0.0, 0.0]
            for j in range(frag_len):
                cm[0] = cm[0] + xyz_HN[j][0]
                cm[1] = cm[1] + xyz_HN[j][1]
                cm[2] = cm[2] + xyz_HN[j][2]
            cm[0] /= float(frag_len)
            cm[1] /= float(frag_len)
            cm[2] /= float(frag_len)
            for j in range(frag_len):
                fastT1FM.SetDArray(j, 0, xyz, xyz_HN[j][0] - cm[0])
                fastT1FM.SetDArray(j, 1, xyz, xyz_HN[j][1] - cm[1])
                fastT1FM.SetDArray(j, 2, xyz, xyz_HN[j][2] - cm[2])

            tensor = fastT1FM.MakeDMatrix(nsets, 8)
            Xaxrh_range = fastT1FM.MakeDMatrix(nsets, 4)
            for i in range(0, nsets):
                fastT1FM.SetDArray(i, 0, Xaxrh_range, 0.05)
                fastT1FM.SetDArray(i, 1, Xaxrh_range, 200.0)
                fastT1FM.SetDArray(i, 2, Xaxrh_range, 0.05)
                fastT1FM.SetDArray(i, 3, Xaxrh_range, 200.0)
            # ****
            chisqr = fastT1FM.rfastT1FM_multi(npts, rMx, rMy, rMz, nsets, frag_len, xyz, pcs, tensor, Xaxrh_range)
            # ****

            saupe_array = []

            for kk in range(nsets):
                temp_saupe = []
                for j in range(3, 8):
                    temp_saupe.append(fastT1FM.GetDArray(kk, j, tensor))
                saupe_array.append(temp_saupe)

            x = fastT1FM.GetDArray(0, 0, tensor)
            y = fastT1FM.GetDArray(0, 1, tensor)
            z = fastT1FM.GetDArray(0, 2, tensor)
            metal_pos = [x + cm[0], y + cm[1], z + cm[2]]

            # Compute and check Axial and Rhombic parameters
            AxRh = calcAxRh(saupe_array)

            #chisqr = checkAxRh(AxRh, chisqr, total_pcs, stage)  # modifies the values of chisqr
            chisqr = checkAxRhCutoffs(AxRh, chisqr, total_pcs, exp_data, stage)  # modifies the values of chisqr
            AxRh.append(metal_pos)  # add metal pos

            # Free memory
            fastT1FM.FreeDMatrix(xyz)
            fastT1FM.FreeDMatrix(pcs)
            fastT1FM.FreeDMatrix(Xaxrh_range)
            fastT1FM.FreeDMatrix(tensor)
        else:
            chisqr = 1.0e+30

        if chisqr < 1.0e+30:

            nchisqr = chisqr / float(total_pcs - (nsets * 5))
            snchisqr = nchisqr / float(math.pow(total_pcs, 1 / 3.0))
            temp_tensor.append([tag, snchisqr, AxRh])


        if stage <= 3:
            if len(temp_tensor) < tag:
                # free memory and return
                fastT1FM.FreeDArray(rMx)
                fastT1FM.FreeDArray(rMy)
                fastT1FM.FreeDArray(rMz)
                tfalse = []
                return tfalse
        if stage == 4:
            if len(temp_tensor) < tag + 1:
                # free memory and return
                fastT1FM.FreeDArray(rMx)
                fastT1FM.FreeDArray(rMy)
                fastT1FM.FreeDArray(rMz)
                tfalse = []
                return tfalse

    fastT1FM.FreeDArray(rMx)
    fastT1FM.FreeDArray(rMy)
    fastT1FM.FreeDArray(rMz)

    return temp_tensor
