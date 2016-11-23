#!/usr/bin/env python

"""
Project_Name: main, File_name: stage0_rdc.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 08/08/2016

prepare data to run with BOSS-R
"""
import multiprocessing

import utility.RDCUtil    as ru
import utility.io_util    as io
import utility.ss_util    as ss

# Parse the input data text file
data = io.readInputDataFiles('input_data.txt')
datatypes = data.keys()
print datatypes

# Read in the amino acid sequence and the Secondary structure assignment
handle, aa_seq = io.readFasta(data['fasta_file'])
ss_seq = ru.matchSeq2SS(aa_seq, data['ss_file'])
print ss_seq

# Generate fuzzy +/-2 SSE combinations
ss_def, ss_combi = ss.genSSCombinations(ss_seq)
io.dumpPickle("ss_profiles.pickle", ss_combi)
print ss_def
# Read the native pdbs that you can exclude from the smotif search
native_pdbs = data['native_pdbs']
native_pdbs = native_pdbs.lower()
native_pdbs = native_pdbs.split()
print native_pdbs

# Read in the axial, Rhombic and chisqr cutoffs for all of the datasets
axrh_cutoff = data['rdc_axrh_cutoff']
axrh_cutoff = axrh_cutoff.split()
axrh_cutoff = [float(i) for i in axrh_cutoff]

chisqr_cutoff = data['rdc_chisqr_cutoff']
chisqr_cutoff = chisqr_cutoff.split()
chisqr_cutoff = [float(i) for i in chisqr_cutoff]

pred_axial = data['predicted_axial']
pred_axial = pred_axial.split()
pred_axial = [float(i) for i in pred_axial]

exp_error = data['exp_error']
exp_error = exp_error.split()
exp_error = [float(i) for i in exp_error]


rmsd_cutoff = data['rmsd_cutoff']
rmsd_cutoff = rmsd_cutoff.split()
rmsd_cutoff = [float(i) for i in rmsd_cutoff]

print rmsd_cutoff

if (len(chisqr_cutoff) != 4) or (len(axrh_cutoff) != 4):
    print "The number of specified cutoffs should be exactly 4"
clash_distance = float(data['clash_distance'])
print 'clash_distance: ', clash_distance

rdcdata = ru.getRdcData(data['rdc_input_files'], ss_seq)
map_route = [[9, 10, 'start'], [8, 9, 'left'], [7, 8, 'left'], [6, 7, 'left'], [5, 6, 'left'], [4, 5, 'left'],
             [3, 4, 'left'], [2, 3, 'left'], [1, 2, 'left'], [0, 1, 'left']]

io.dumpPickle("rdc_route.pickle", map_route)
database_cutoff = data['database_cutoff']

data_dict = {'ss_seq': ss_seq, 'rdc_data': rdcdata, 'aa_seq': aa_seq, 'natives': native_pdbs, \
             'clash_distance': clash_distance, 'database_cutoff': database_cutoff, \
             'rdc_axrh_cutoff': axrh_cutoff, 'rdc_chisqr_cutoff': chisqr_cutoff, 'rmsd_cutoff': rmsd_cutoff,
             'pred_axial': pred_axial, 'exp_error': exp_error}

io.dumpPickle("exp_data.pickle", data_dict)

fout = open("run.sh", 'w')
fout.write("#!/bin/bash\n")
ncpus = multiprocessing.cpu_count()

for i in range(0, len(map_route)):
    smotif = map_route[i]
    print smotif
    if i == 0:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage1_mpi_run.py\n"
        print run_line
        fout.write(run_line)
    elif i == 1:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage2_mpi_run.py 10\n"
        print run_line
        fout.write(run_line)
    else:
        run_line = "mpirun -np " + str(ncpus) + " python ../../main/stage3x_mpi_run.py 10\n"
        print run_line
        fout.write(run_line)

fout.close()
exit()
