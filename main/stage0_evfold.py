#!/usr/bin/env python

"""
Project_Name: main, File_name: stage0_evfold.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/06/2016

prepare data to run with BOSS-Evo

"""


import multiprocessing
import sys
sys.path.append('/home/kalabharath/projects/boss-evo/zinr/')
import utility.evfold as evfold
import utility.io_util    as io
import utility.ss_util    as ss


data = io.readInputDataFiles('input_data.txt')
aa_seq, ss_seq = evfold.parseSS(data['evfold_ss'])
database_cutoff = data['database_cutoff']
rmsd = data['rmsd']
clash_distance = data['clash_distance']


ss_def, ss_combi = ss.genSSCombinations(ss_seq)
print ss_seq
print ss_def
native_pdbs = evfold.parseNatives(data['native_pdbs'])
contact_matrix, plm_score_matrix, contact_ss_matrix = evfold.parseContacts(data['evfold_plm'], ss_combi, ss_def,
                                                                           len(ss_seq), cutoff_score=0.25)
map_route = evfold.getRoute2(ss_def, contact_matrix)
print map_route
data_dict = {'ss_seq': ss_seq, 'contact_matrix': contact_matrix, 'plm_scores': plm_score_matrix, 'aa_seq': aa_seq,
             'natives': native_pdbs, 'database_cutoff': database_cutoff, 'rmsd': float(rmsd),
             'clash_distance': float(clash_distance)}

"""
Dump all of the data in relavant pickle files
however, have to streamline this to a much better way of doing it
"""
io.dumpPickle("exp_data.pickle", data_dict)
io.dumpPickle("ss_profiles.pickle", ss_combi)
io.dumpPickle("contacts_route.pickle", map_route)

fout = open("run.sh", 'w')
ncpus = multiprocessing.cpu_count()

for i in range(0, len(map_route)):
    smotif = map_route[i]
    print smotif
    if i == 0:
        run_line = "mpirun -np " + str(ncpus) + "python stage1_mpi_run.py"
        print run_line
        fout.write(run_line)
    elif i == 1:
        run_line = "mpirun -np " + str(ncpus) + "python stage2_mpi_run.py"
        print run_line
        fout.write(run_line)
    else:
        run_line = "mpirun -np " + str(ncpus) + "python stage3x_mpi_run.py"
        print run_line
        fout.write(run_line)

fout.close()
exit()
