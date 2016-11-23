#!/usr/bin/env python

"""
Project_Name: main/filters/contacts, File_name: evfoldContacts.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 16/06/16

Score evolutionary couplings retrieved from EVFold server
"""

import contacts_filter
from filters.constraints.looplengthConstraint import get_dist

def calcFmeasure(gbar, ganoe):
    tp, fp, fn = 0.0, 0.0, 0.0

    for entry in ganoe:
        if entry in gbar:
            tp += 1
        else:
            fn += 1
    for entry in gbar:
        if entry in ganoe:
            pass
        else:
            fp += 1
    if tp:
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        return precision
        #return (2 * precision * recall) / (precision + recall)
    else:
        return False


def calcPLMscore(gbar, plm_score_matrix):
    score = 0.0
    for entry in gbar:
        score = score + plm_score_matrix[entry[0], entry[1]]
    return score


def s1EVcouplings(s1_def, s2_def, smotif, contact_matrix, plm_score_matrix, contacts_cutoff):
    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param contact_matrix:
    :param plm_score_matrix:
    :return:
    """

    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    # print ss1_list, len(ss1_list)
    # print ss2_list, len(ss2_list)
    # print smotif_ss1, len(smotif_ss1)
    # print smotif_ss2, len(smotif_ss2)

    # find out all possible contacts in the Smotif
    contacts_found = []
    for res in smotif_ss1:
        for entry1 in smotif[1]:
            if entry1[2] == 'CA' and entry1[0] == res:
                coo1 = [entry1[3], entry1[4], entry1[5]]
                for entry2 in smotif[2]:
                    if entry2[2] == 'CA':
                        coo2 = [entry2[3], entry2[4], entry2[5]]
                        dist = contacts_filter.get_distance(coo1, coo2)
                        if dist < contacts_cutoff:
                            contacts_found.append(
                                (ss1_list[smotif_ss1.index(entry1[0])], ss2_list[smotif_ss2.index(entry2[0])]))
                            """
                            print ss1_list[smotif_ss1.index(entry1[0])],
                            print ss2_list[smotif_ss2.index(entry2[0])]
                            print entry1[0], entry2[0]
                            """

    # total number of contacts according to the EV coupling matrix is
    # why do i have to do this every other time, i can save and read from the disk
    # or is it much more easier to compute than saving and reading again and again and again
    # or i can supply the data to this function thereby calculating only once

    total_contacts = []
    for entry1 in ss1_list:
        for entry2 in ss2_list:
            if contact_matrix[entry1, entry2]:
                total_contacts.append((entry1, entry2))

    fmeasure = calcFmeasure(contacts_found, total_contacts)

    if fmeasure:
        plm_score = calcPLMscore(contacts_found, plm_score_matrix)
        return fmeasure, plm_score
    else:
        return 0.0, 0.0


def getCAandresi(frag):
    x, y, z = [], [], []
    resi = []
    for i in range(0, len(frag[0])):
        if frag[3][i] == 'CA':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
            resi.append(frag[4][i])
    return resi, [x, y, z]


def s2EVcouplings(transformed_coors, native_sse_order, contact_matrix, plm_score_matrix, contacts_cutoff):
    # print native_sse_order
    import copy
    sse_coors = copy.deepcopy(transformed_coors)

    contacts_found = []
    total_contacts = []
    for i in range(0, len(sse_coors) - 1):

        res_c, ca1 = getCAandresi(sse_coors[i])
        res_n, ca2 = getCAandresi(sse_coors[i + 1])
        ss1_list = range(native_sse_order[i][4], native_sse_order[i][5] + 1)
        ss2_list = range(native_sse_order[i + 1][4], native_sse_order[i + 1][5] + 1)

        for res1 in ss1_list:
            # print ss1_list, res1, res_c, res_c[ss1_list.index(res1)]
            # print ss1_list.index(res1)
            ca_res1 = [ca1[0][ss1_list.index(res1)], ca1[1][ss1_list.index(res1)], ca1[2][ss1_list.index(res1)]]
            for res2 in ss2_list:
                ca_res2 = [ca2[0][ss2_list.index(res2)], ca2[1][ss2_list.index(res2)], ca2[2][ss2_list.index(res2)]]
                dist = get_dist(ca_res1, ca_res2)

                if contact_matrix[res1, res2] or contact_matrix[res2, res1]:
                    total_contacts.append((res1, res2))

                if dist < contacts_cutoff:
                    contacts_found.append((res1, res2))

    # print len(total_contacts), total_contacts
    # print len(contacts_found), contacts_found

    fmeasure = calcFmeasure(contacts_found, total_contacts)
    if fmeasure:
        plm_score = calcPLMscore(contacts_found, plm_score_matrix)
        return fmeasure, plm_score
    else:
        return 0.0, 0.0
