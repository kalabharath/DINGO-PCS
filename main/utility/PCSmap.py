#!/usr/bin/env python

"""
Project_Name: main, File_name: PCSmap.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 1/04/15 , Time:10:53 AM
"""
# TODO write comments for each def

def getNspcs(sstart, send, pcs_data):
    nspcs = []

    for dataset in pcs_data:
        tnspcs = 0
        for i in xrange(0, len(dataset)):
            if sstart <= i <= send:
                for entry in dataset[i]:
                    if entry == 999.999:
                        pass
                    else:
                        tnspcs += 1
            else:
                pass
        nspcs.append(tnspcs)
    return nspcs


def SearchPcs(ss_array, pcs_data):
    pcs_array = []
    for i in xrange(0, len(ss_array)):
        nspcs = getNspcs(ss_array[i][3], ss_array[i][4], pcs_data)
        nspcs.sort()
        nspcs.reverse()
        t = 0
        for j in xrange(0, 2):
            t = t + nspcs[j]
        pcs_array.append(t)
    return pcs_array


def get_ij(pcs_array):
    """
        return the index of SS carrying largest PCSs/SS
    """
    ntpcs = 0

    for i in xrange(0, len(pcs_array) - 1):
        if i == 0:
            tpcs = pcs_array[i] + pcs_array[i + 1]
            ti = i
            tj = i + 1
        else:
            ntpcs = pcs_array[i] + pcs_array[i + 1]
        if tpcs >= ntpcs:
            pass
        else:
            tpcs = ntpcs
            ti = i
            tj = i + 1
    return ti, tj


def getRoute(ss_seq, pcsdata):
    """
	returns Smotif search route based PCSs/SS
	"""
    import ss_util as ssutil

    ss_array, lloop = ssutil.genSSDef(ss_seq)

    pcs_array = SearchPcs(ss_array, pcsdata)
    control = 0
    i = 0
    j = 0

    map_route = []
    while not (control == len(ss_array) - 1):
        # [20, 64, 28, 24, 56, 28, 24, 40, 24, 24, 88] 11
        if control == 0:
            control += 1
            i, j = get_ij(pcs_array)  # get the largest PCS smotif
            map_route.append([i, j, 'start'])
        else:
            if i == 0:
                ti = j
                j += 1
                control += 1
                direction = 'right'
                map_route.append([ti, j, direction])
            elif j == len(pcs_array) - 1:
                tj = i
                i -= 1
                control += 1
                direction = 'left'
                map_route.append([i, tj, direction])
            else:
                if pcs_array[i - 1] >= pcs_array[j + 1]:
                    tj = i
                    i -= 1
                    direction = 'left'
                    control += 1
                    map_route.append([i, tj, direction])
                else:
                    ti = j
                    j += 1
                    control += 1
                    direction = 'right'
                    map_route.append([ti, j, direction])
    return map_route