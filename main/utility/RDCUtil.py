#!/usr/bin/env python

"""
Project_Name: main/RDCutil, File_name: RDCUtil.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 8/8/16 , Time:4:34 PM
"""


def matchSeq2SS(aa_seq, ssfile):
    print aa_seq
    print ssfile
    raw_ss = []
    with open(ssfile) as fin:
        lines = fin.readlines()

    for i in range(0, len(lines)):
        if lines[i] == 'FORMAT %4d %1s %2d %2d %8.3f %8.3f %8.3f %4.2f %s\n':
            print lines[i]
            j = i
            for j in range(j, len(lines)):

                content = lines[j].split()
                if len(content) == 9:
                    raw_ss.append(content)
            break

    print len(aa_seq), len(raw_ss)
    diff = len(raw_ss) - len(aa_seq)
    print 'diff', diff
    if diff > 0:
        for i in range(0, diff + 1):
            t_aa = ''
            t_ss = ''
            for j in range(i, i + len(aa_seq)):
                t_aa = t_aa + raw_ss[j][1]
                t_ss = t_ss + raw_ss[j][-1]
            print t_aa, t_ss
            if t_aa == aa_seq:
                print t_aa, len(t_aa)
                # REMARK     h-Helix    e-Strand   c-Coil (Sequence based)
                t_ss = t_ss.replace('c', 'L')
                t_ss = t_ss.replace('e', 'E')
                t_ss = t_ss.replace('h', 'H')
                print t_ss, len(t_ss)
                return t_ss
    else:

        t_aa = ''
        t_ss = ''
        for j in range(0, len(aa_seq)):
            t_aa = t_aa + raw_ss[j][1]
            t_ss = t_ss + raw_ss[j][-1]
        t_ss = t_ss.replace('c', 'L')
        t_ss = t_ss.replace('e', 'E')
        t_ss = t_ss.replace('h', 'H')
        print t_ss, len(t_ss)
        return t_ss

def FormatRdc(seqlen, rdcfile):
    """
    parses rdc from .npc file.
    should be in the format #['179', 'H', '179', 'N', '16.042', '0.0']
    the rdcs are returned as a dict with res_no as key and rdc def as value.
    """

    import io_util as io
    rdc_l = io.readFile(rdcfile)
    rdcs = {}
    for l in rdc_l:
        r1, v1, r2, v2, rdc, tol = l.split()
        rdcs.setdefault(int(r1), []).append([int(r1), v1, int(r2), v2, float(rdc)])
    return rdcs

def getRdcData(rdc_files, ss_seq):
        """
        Parse RDCs from files
        """
        rdc_files = rdc_files.split()
        rdc_data = []
        for j in range(0, len(rdc_files)):
            rdc_data.append(FormatRdc(len(ss_seq), rdc_files[j]))
        return rdc_data
