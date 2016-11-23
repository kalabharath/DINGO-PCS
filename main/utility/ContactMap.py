#!/usr/bin/env python

"""
Project_Name: main, File_name: ContactMap.py 
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 1/04/15 , Time:12:54 PM
"""


def getSmotif(ss1, ss2):
    """
	returns end residues of pairs of SSE
	:param ss1:
	:param ss2:
	:return:
	"""
    # ss_element = ['strand', 8, 5, 6, 13],
    s1_start, s1_end = ss1[3], ss1[4]
    s2_start, s2_end = ss2[3], ss2[4]
    return [[s1_start, s1_end], [s2_start, s2_end]]


def getTopSmotif(rank_seq):
    """
	sorts the Smotif indices based on number of contacts
	:param rank_seq:
	:return: top Smotif index
	"""
    keys = rank_seq.keys()
    number_of_contacts = 0
    top_rank = 0
    for key in keys:
        if rank_seq[key] > number_of_contacts:
            number_of_contacts = rank_seq[key]
            top_rank = key
    return top_rank


def getTopRank(ss_def, contacts_def):
    """
	Identifies pairs of SSEs which have largest intra contacts between them
	:param ss_def:
	:param contacts_def:
	:return:  SSE indices
	#TO DO could be written better
	"""
    rank_seq = {}
    contacts_true = contacts_def.keys()
    for i in range(1, len(ss_def)):
        smotif = getSmotif(ss_def[i - 1], ss_def[i])
        no_of_contacts = 0
        for j in range(smotif[0][0], smotif[0][1] + 1):
            if j in contacts_true:
                contacts = contacts_def[j]
                for contact in contacts:
                    if contact in range(smotif[1][0], smotif[1][1] + 1):
                        no_of_contacts += 1

        rank_seq[i - 1] = no_of_contacts
    top_rank = getTopSmotif(rank_seq)

    return top_rank, top_rank + 1


def total_data(map_index, next_sse_index, ss_def, contacts_def):
    """
	Count the total number of contacts that the prospective SSE makes with the rest of selected SSEs
	:param map_index:
	:param next_sse_index:
	:param ss_def:
	:param contacts_def:
	:return: total number of contacts
	"""
    no_of_contacts = 0
    contacts_true = contacts_def.keys()
    start, end = ss_def[next_sse_index][3], ss_def[next_sse_index][4]
    for i in range(start, end + 1):
        if i in contacts_true:
            contacts = contacts_def[i]
            for contact in contacts:
                for index in map_index:
                    tstart, tend = ss_def[index][3], ss_def[index][4]
                    if contact in range(tstart, tend + 1):
                        no_of_contacts += 1
    return no_of_contacts


def nextSS(map_route, ss_def, contacts_def):
    """
	Identifies the next SSE to be in the  map or not
	:param map_route:
	:param ss_def:
	:param contacts_def:
	:return: Smotif SSE indices and direction
	"""
    map_index = []
    for i in range(0, len(map_route)):
        for j in range(0, 2):
            ti = map_route[i][j]
            if ti in map_index:
                pass
            else:
                map_index.append(ti)
    map_index.sort()

    temp_numcontacts = []
    for i in range(-1, 2, 2):  # HARDCODE, next possible SSE two directions
        if i < 0:
            if map_index[0] > 0:
                next_sse_index = map_index[0] + i
                tnum = total_data(map_index, next_sse_index, ss_def, contacts_def)

            else:
                tnum = -1
        else:
            if map_index[-1] <= len(ss_def) - 1:
                next_sse_index = map_index[-1] + i
                tnum = total_data(map_index, next_sse_index, ss_def, contacts_def)

            else:
                tnum = -1

        temp_numcontacts.append(tnum)

    if temp_numcontacts[0] >= temp_numcontacts[1]:
        new_direction = 'left'
        ti = map_index[0] - 1
        tj = map_index[0]
        return ti, tj, new_direction
    else:
        new_direction = 'right'
        ti = map_index[-1]
        tj = map_index[-1] + 1
        return ti, tj, new_direction


def getContactRoute(ss_def, contacts_def):
    """
	Determine the sequnce of SSEs to be assembled using contact information.
	:param ss_def:
	:param contacts_def:
	:return: map_route
	"""
    control = 0  # to control the while loop until the length of all SSEs
    i = 0  # i, j are the indices of the SSE array
    j = 0
    map_route = []
    while control != len(ss_def) - 1:

        # [20, 64, 28, 24, 56, 28, 24, 40, 24, 24, 88] 11
        if control == 0:
            control += 1
            i, j = getTopRank(ss_def, contacts_def)  # get the largest intra contact smotif
            map_route.append([i, j, 'start'])
        else:
            if i == 0:
                ti = j
                j += 1
                control += 1
                direction = 'right'
                map_route.append([ti, j, direction])
            elif j == len(ss_def) - 1:  # -1 to deal with array index
                tj = i
                i -= 1
                control += 1
                direction = 'left'
                map_route.append([i, tj, direction])
            else:
                ti, tj, new_direction = nextSS(map_route, ss_def, contacts_def)
                control += 1
                map_route.append([ti, tj, new_direction])
    return map_route

