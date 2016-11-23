#!/usr/bin/env python

"""
Project_Name: rmsd, File_name: test_qcp_rmsd.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 17/03/15 , Time:10:57 AM
"""

import pickle

import qcprot

# extern double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight)
# extern int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore)
#hh_9_37.db
#['2exwB00', '18', '26', '33', '69'] ['3omaC00', '72', '80', '92', '128']
#  RMSD from pymol 33-68 and 92-128 is 3.478


pickle_array = open("/home/kalabharath/zinr/smotif_database/hh_9_37.db", 'r')
smotif_array = pickle.load(pickle_array)
pickle_array.close()
print smotif_array[0][0][0], smotif_array[1][0][0]
smotif_1 = smotif_array[0]

smotif_2 = smotif_array[1]
print smotif_1[0][0], len(smotif_1[0][1]), len(smotif_1[0][2])
frag_ax, frag_ay, frag_az = [], [], []
s2_coords = smotif_1[0][2]
for i in range(0, len(s2_coords)):
    #[69, 'VAL', 'C', 14.469, 26.646, 64.269]
    if s2_coords[i][2] == 'CA':
        frag_ax.append(s2_coords[i][3])
        frag_ay.append(s2_coords[i][4])
        frag_az.append(s2_coords[i][5])
frag_a = [frag_ax, frag_ay, frag_az]


frag_bx, frag_by, frag_bz = [], [], []
s2_coords = smotif_2[0][2]
for i in range(0, len(s2_coords)):
    #[69, 'VAL', 'C', 14.469, 26.646, 64.269]
    if s2_coords[i][2] == 'CA':
        #print s2_coords[i]
        frag_bx.append(s2_coords[i][3])
        frag_by.append(s2_coords[i][4])
        frag_bz.append(s2_coords[i][5])
frag_b = [frag_bx, frag_by, frag_bz]

"""
double          rmsd, x, y, z, euc_dist;
    double        **frag_a, **frag_b;
    int             len = 7;
    double          rotmat[9];
    rmsd = CalcRMSDRotationalMatrix((double **) frag_a, (double **) frag_b, len, rotmat, NULL);
"""

fraglen = 37
xyz1 = qcprot.MakeDMatrix(3, fraglen)
xyz2 = qcprot.MakeDMatrix(3, fraglen)

for i in range(0, fraglen):
    qcprot.SetDArray(0, i, xyz1, frag_a[0][i])
    qcprot.SetDArray(1, i, xyz1, frag_a[1][i])
    qcprot.SetDArray(2, i, xyz1, frag_a[2][i])
for i in range(0, fraglen):
    qcprot.SetDArray(0, i, xyz2, frag_b[0][i])
    qcprot.SetDArray(1, i, xyz2, frag_b[1][i])
    qcprot.SetDArray(2, i, xyz2, frag_b[2][i])

rot_matrix = qcprot.MakeDMatrix(3,3)

rot = qcprot.MakeDvector(9)




rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, fraglen, rot)

for i in range(0,9):
    print qcprot.GetDvector(i,rot)
"""
for i in range(0, fraglen):
    print qcprot.GetDArray(0, i, xyz1), qcprot.GetDArray(1, i, xyz1), qcprot.GetDArray(2, i, xyz1)
    print qcprot.GetDArray(0, i, xyz2), qcprot.GetDArray(1, i, xyz1), qcprot.GetDArray(2, i, xyz1)
"""




print rmsd
