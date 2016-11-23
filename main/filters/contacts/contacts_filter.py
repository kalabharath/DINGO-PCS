#!/usr/bin/env python

"""
Project_Name: main/filters/sequence, File_name: sequence_similarity.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 16/04/15 , Time: 04:52 PM

Whether the Smotifs satisfy the input predicted contacts

*********************
ERROR: Needs review and fixing of directionality
*********************
"""

def get_distance(coo1, coo2):
    import math
    x1,y1,z1=coo1[0],coo1[1],coo1[2]
    x2,y2,z2=coo2[0],coo2[1],coo2[2]
    return math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))


def ContactPredicition(s1_def, s2_def, smotif, exp_data):

    """

    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """
    #     print s1_def, s2_def ['strand', 9, 4, 4, 89, 97] ['strand', 9, 4, 4, 102, 110]

    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    # print ss1_list, ss2_list
    # print smotif_ss1, smotif_ss2
    # print smotif[0][0]

    contacts = exp_data['contacts']

    residue_list = []
    contacts_in_smotif = 0
    for contact in contacts:
        if contact[0] in ss1_list and contact[1] in ss2_list:
            contacts_in_smotif +=1
            res1_index = ss1_list.index(contact[0])
            res2_index = ss2_list.index(contact[1])
            #print contact, res1_index, res2_index

            smotif_res1 = smotif_ss1[res1_index]
            smotif_res2 = smotif_ss2[res2_index]
            #print smotif_res1, smotif_res2

            for res in smotif[1]:
                if res[2] == 'CA' and res[0] == smotif_res1:
                    coo1 = [res[3], res[4], res[5]]

            for res in smotif[2]:
                if res[2] == 'CA' and res[0] == smotif_res2 :
                   coo2 = [res[3], res[4], res[5]]
            if coo1 and coo2 :
                dist = get_distance(coo1,coo2)
                if dist <= contact[2] and dist > 2.0:
                    residue_list.append(True)
                else:
                    residue_list.append(False)
            else:
                print "Error in progressing"
    hits =0
    for entry in residue_list:
        if entry:
            hits +=1

    if hits:
        return contacts_in_smotif, float(hits)/float(contacts_in_smotif)*100.00
    else:
        return contacts_in_smotif, 9999.999

def getCA(coo_array):
    """

    :param coo_array:
    :return:
    """

    xt, yt, zt = [],[],[]
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]
        if atom == 'CA':
            xt.append(x)
            yt.append(y)
            zt.append(z)
    return [xt, yt, zt]

def coorCAdict(coo_arrays, sse_list):
    """

    :param coo_arrays:
    :param sse_list:
    :return:
    """
    ca_dict={}
    for i in range(0,len(sse_list)):
        sse_def = sse_list[i]
        ca_array = getCA(coo_arrays[i])
        sse_range = range(sse_def[4], sse_def[5]+1)
        for j in range(0, len(sse_range)):
            ca_dict[sse_range[j]] = [ca_array[0][j], ca_array[1][j], ca_array[2][j]]
    return ca_dict

def S2ContactPredicition(coo_arrays, sse_list, exp_data):

    """
    :param s1_def:
    :param s2_def:
    :param smotif:
    :param threshold:
    :return:
    """

    contacts = exp_data['contacts']

    ca_dict = coorCAdict(coo_arrays, sse_list)

    ca_res = ca_dict.keys()

    residue_list =[]
    total_contacts = 0
    for contact in contacts:
        res1, res2 = contact[0], contact[1]
        if res1 in ca_res and res2 in ca_res:
            total_contacts +=1
            dist = get_distance(ca_dict[res1], ca_dict[res2])
            if dist <= contact[2] and dist > 2.0 :
                residue_list.append(True)
            else:
                residue_list.append(False)

    hits =0
    for entry in residue_list:
        if entry:
            hits +=1
    if hits:
        return total_contacts, float(hits)/float(total_contacts)*100.00
    else:
        return total_contacts, 9999.999