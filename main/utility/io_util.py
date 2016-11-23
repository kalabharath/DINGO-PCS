#!/usr/bin/env python
import os
import pickle

"""
Project_Name: setup, File_name: util.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 19/03/15 , Time:3:24 PM
Utility file for frequently used functions
"""


def readInputDataFiles(filename):
    lines = readFile(filename)
    data_dict = {}
    for line in lines:
        if line[0] != '#' and line != '\n':
            data_type, data_file = line.split('=')
            #print data_type.strip(), data_file.strip()
            data_dict[data_type.strip()] = data_file.strip()
    if len(data_dict.keys()) == 1:
        return data_dict
    else:
        return data_dict


def readFile(filename):
    """
	read file and return the lines as an array
	:param filename:
	:return:
	"""
    if os.path.isfile(filename):
        with open(filename, 'r') as fin:
            return fin.readlines()
    else:
        return False



def writeFile(array_of_lines, filename):
    """
	write output to a file
	:param array_of_lines:
	:param filename:
	:return:
	"""
    fout = open(filename, 'w')
    for line in array_of_lines:
        fout.write(line)
        fout.write("\n")
    fout.close()
    return True


def dumpPickle(filename, data):
    """
    dump data as a pickle object, usually arrays
    :param filename:
    :param data:
    :return:
    """
    with open(filename, 'w') as f:
        pickle.dump(data, f)
    f.close()
    return True


def readPickle(filename):
    """
    read pickle objects
    :param filename:
    :return:
    """
    if os.path.isfile(filename):
        fin = open(filename, 'r')
        data = pickle.load(fin)
        fin.close()
        return data
    else:
        return False


def readFasta(filename):
    """
	reads in FastA file and returns seq output
	:param filename:
	:return:
	"""
    with open(filename, 'r') as fin:
        lines = fin.readlines()

    header = lines[0].rstrip()
    seq = ''
    for i in range(1, len(lines)):
        seq += lines[i].rstrip()
    return header, seq


def readPsiPred(filename):
    """
	read in ss_seq as predicted by PsiPred
	:param filename:
	:return: ss_seq
	"""
    with open(filename, 'r') as fin:
        lines = fin.readlines()

    ss_seq = ''

    for line in lines:
        if line[0:5] == 'Pred:':
            ss_seq = ss_seq + line[6:-1]
    return ss_seq


def readContacts(filename, probability):
    """
    Parses MetaPSICOV's contact file and returns the
	:param filename:
	:param probability:
	:return contacts array
	"""

    with open(filename, 'r') as fin:
        lines = fin.readlines()
    contacts = []
    contacts_seq = {}
    for line in lines:
        if line != '\n':
            res1, res2, zero, distance, precision = line.split()
            precision = round(float(precision), 1)
            if precision >= probability:
                contacts.append([int(res1), int(res2), float(distance), float(precision)])
                contacts_seq.setdefault(int(res1), []).append(int(res2))
                contacts_seq.setdefault(int(res2), []).append(int(res1))
    return contacts, contacts_seq


def auto_readContacts(probability):
    """
    Parses MetaPSICOV's contact file and returns the
	:param filename:
	:param probability:
	:return contacts array
	"""
    import glob

    filename = glob.glob("*.metapsicov")
    # filename = readInputDataFiles('input_data.txt',['contacts_file'])

    with open(filename[0], 'r') as fin:
        lines = fin.readlines()
    contacts = []
    contacts_seq = {}
    for line in lines:
        if line != '\n':
            res1, res2, zero, distance, precision = line.split()
            precision = round(float(precision), 1)
            if precision >= probability:
                contacts.append([int(res1), int(res2), float(distance), float(precision)])
                contacts_seq.setdefault(int(res1), []).append(int(res2))
                contacts_seq.setdefault(int(res2), []).append(int(res1))

    return contacts
    # return contacts, contacts_seq


def FormatPcs(seqlen, nsets, filenames):
    import string

    pcsdata = []
    for i in range(seqlen):
        pcsdata.append([999.999] * nsets)
    for j in range(len(filenames)):
        f = open(filenames[j], 'r')
        line = f.readline()
        while line:
            xx = string.split(line[:-1])
            resno = int(xx[0]) - 1
            try:
                pcsdata[resno][j] = float(xx[2])
            except:
                pass
            line = f.readline()
        f.close()
    return pcsdata


def getPcsTagInfo(ss_seq, brokerfile):
    brokerin = open(brokerfile, 'r')
    brokerlines = brokerin.readlines()
    brokerin.close()
    ntags = 0
    tsets = []
    npc_files = {}
    for line in brokerlines:
        if line[0:7] == 'CLAIMER':
            ntags += 1
    for i in range(ntags):
        temp = 0
        for l in brokerlines:
            if l[0:18] == 'TS' + str(i + 1) + '_PCS_INPUT_FILE':
                npc_files.setdefault(i, []).append(l[19:-1])
                temp = temp + 1
        tsets.append(temp)
    pcs_data = []

    for j in range(0, ntags):
        pcs_data.append(FormatPcs(len(ss_seq), tsets[j], npc_files[j]))
    return pcs_data