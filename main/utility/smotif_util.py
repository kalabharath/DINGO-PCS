#!/usr/bin/env python

"""
Project_Name: main/utility, File_name: smotif_util.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:12:54 PM
"""

import utility.io_util as io


def getSmotif(s1, s2):
    """
    Return Smotif based on SS definition

    ('helix-helix', 6, 23, 4)
    new_smotif_file = hh_6_23.db
    ['helix', 6, 9, 4, 146, 151] ['helix', 23, 4, 1, 156, 178]
    """
    print s1, s2
    if s1[0] == 'helix':
        s1_type = 'h'
    if s1[0] == 'strand':
        s1_type = 's'
    if s2[0] == 'helix':
        s2_type = 'h'
    if s2[0] == 'strand':
        s2_type = 's'
    smotif_type = s1_type + s2_type
    s1_len = s1[1]
    s2_len = s2[1]
    smotif = [smotif_type, s1_len, s2_len]
    return smotif


def readSmotifDatabase(smotif, *database_cutoff):
    """

    :param smotif:
    :return:
    """
    # TODO option to parse in database path
    import os

    if database_cutoff:
        # ['', 'home', 'kalabharath', 'projects', 'boss-evo', 'zinr', 'utility']
        cwd = (os.path.dirname(os.path.realpath(__file__))).split("/")
        root_dir = ''
        for entry in cwd[:-2]:
            if entry == '':
                pass
            else:
                root_dir = root_dir + '/' + entry
        smotif_db_path = root_dir + '/databases/database_cutoff_' + database_cutoff[0] + '/'
    else:
        # Backwards compatible with my old code
        # I should retire this soon and keep everything clean

        smotif_db_path = "/short/xc4/kbp502/zinr/main/smotif_cen_db/"
        if os.path.exists(smotif_db_path):
            pass
        else:
            smotif_db_path = "/home/kalabharath/zinr/main/smotif_cen_db/"

    file_name = smotif[0] + "_" + str(smotif[1]) + "_" + str(smotif[2]) + ".db"
    fin = smotif_db_path + file_name
    smotif_data = io.readPickle(fin)
    return smotif_data


def orderSeq(previous_smotif, current_seq, direction):
    """
    :param previous_seq:
    :param current_seq:
    :param direction:
    :return:
    """
    previous_seq = ''

    for entry in previous_smotif:
        if entry[0] == 'seq_filter':
            seq_filter = entry
            previous_seq = seq_filter[1]

    if direction == 'left':
        concat_seq = current_seq + previous_seq
    else:
        concat_seq = previous_seq + current_seq

    return concat_seq


def orderCATH(previous_smotif, current_smotif, direction):
    """

    :param previous_smotif:
    :param current_smotif:
    :param direction:
    :return:
    """
    previous_cath = []
    for entry in previous_smotif:
        if entry[0] == 'cathcodes':
            t_cath = entry[1]
            previous_cath = t_cath[:]

    if direction == 'left':
        previous_cath.insert(0, current_smotif)
    else:
        previous_cath.append(current_smotif)

    return previous_cath
