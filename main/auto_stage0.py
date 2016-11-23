#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_0.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 31/03/15 , Time:10:56 AM

Prepare all the relavant files for stage1 & 2

"""
import os
import sys

sys.path.append("../../main")
# sys.path.append('/short/xc4/kbp502/BOSS/zinr')

import utility.io_util    as io
import utility.ss_util    as ss
import utility.PCSmap     as PCSmap


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


data = io.readInputDataFiles('input_data.txt')

datatypes = data.keys()
handle, aa_seq = io.readFasta(data['fasta_file'])

ss_seq = matchSeq2SS(aa_seq, data['ss_file'])
print ss_seq

# ss_seq = io.readPsiPred(psipred_file)

ss_def, ss_combi = ss.genSSCombinations(ss_seq)
# ss_def, ss_combi = ss.genSSCombinations2(ss_seq)

print ss_combi
io.dumpPickle("ss_profiles.pickle", ss_combi)

# Read in contacts at a given confidence level
if 'contacts_file' in datatypes:
    contacts, contacts_seq = io.readContacts(contactsfile, probability=0.7)

native_pdbs = data['native_pdbs']

native_pdbs = native_pdbs.lower()

native_pdbs = native_pdbs.split()

print native_pdbs

axrh_cutoff = data['axrh_cutoff']
axrh_cutoff = axrh_cutoff.split()
axrh_cutoff = [float(i) for i in axrh_cutoff]

chisqr_cutoff = data['chisqr_cutoff']
chisqr_cutoff = chisqr_cutoff.split()
chisqr_cutoff = [float(i) for i in chisqr_cutoff]

rmsd_cutoff = data['rmsd_cutoff']
rmsd_cutoff = rmsd_cutoff.split()
rmsd_cutoff = [float(i) for i in rmsd_cutoff]

print rmsd_cutoff

if (len(chisqr_cutoff) != 4) or (len(axrh_cutoff) != 4):
    print "The number of specified cutoffs should be exactly 4"

clash_distance = float(data['clash_distance'])

print 'clash_distance: ', clash_distance

# read in PCS data from .npc file from Rosetta's broker file format
pcsdata = io.getPcsTagInfo(ss_seq, data['pcs_broker'])

map_route = PCSmap.getRoute(ss_seq, pcsdata)

print map_route, '1'

io.dumpPickle("pcs_route.pickle", map_route)
# io.dumpPickle("contact_route.pickle", rank_ss)
database_cutoff = data['database_cutoff']

data_dict = {'ss_seq': ss_seq, 'pcs_data': pcsdata, 'aa_seq': aa_seq, 'natives': native_pdbs, \
             'clash_distance': clash_distance, 'database_cutoff': database_cutoff, \
             'axrh_cutoff': axrh_cutoff, 'chisqr_cutoff': chisqr_cutoff, 'rmsd_cutoff': rmsd_cutoff}

io.dumpPickle("exp_data.pickle", data_dict)

res_sofar = 0

total_ss_res = 0

for ss in ss_def:
    tres = ss[1]
    total_ss_res += tres

fout = open("concat.txt", 'w')

for i in range(0, len(map_route)):
    smotif = map_route[i]
    if i == 0:
        print smotif
        s1 = smotif[0]
        s2 = smotif[1]
        print s1, s2
        # print ss_combi[s1][0][1],
        # print ss_combi[s2][0][1]
        res_sofar = ss_combi[s1][0][1] + ss_combi[s2][0][1]
        percent = (res_sofar / float(total_ss_res)) * 100
        print percent
        outline = 'Executable="mpirun -np 128 python stage1_mpi_run.py"\nrun="$Executable"\necho $run\n$run\n'
        print outline
        fout.write(outline)
    else:
        print smotif
        if smotif[-1] == 'right':
            s2 = smotif[1]
            res_sofar = res_sofar + ss_combi[s2][0][1]
            percent = (res_sofar / float(total_ss_res)) * 100
            print percent
            if i == 1:
                num_of_models = (len(map_route) - i) * 125
                print num_of_models
                outline = 'Executable="mpirun -np 128 python stage2_mpi_run.py ' + str(
                    num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                print outline
                fout.write(outline)
                continue

            if 25.0 < percent < 50.0:
                if i >= 2:
                    num_of_models = (len(map_route) - i) * 125
                    print "num_of_models", num_of_models
                    outline = 'Executable="mpirun -np 128 python stage3x_mpi_run.py ' + str(
                        num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                    print outline
                    fout.write(outline)
                else:
                    num_of_models = (len(map_route) - i) * 125
                    print "num_of_models", num_of_models
                    outline = 'Executable="mpirun -np 128 python stage3x_mpi_run.py ' + str(
                        num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                    print outline
                    fout.write(outline)
            if 50.0 < percent < 75.0:
                num_of_models = (len(map_route) - i) * 125
                print "num_of_models", num_of_models
                outline = 'Executable="mpirun -np 128 python stage3x_mpi_run.py ' + str(
                    num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                print outline
                fout.write(outline)
            if percent > 75.0:
                num_of_models = (len(map_route) - i) * 125
                print "num_of_models", num_of_models
                outline = 'Executable="mpirun -np 128 python stage4x_mpi_run.py ' + str(
                    num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                print outline
                fout.write(outline)
        else:
            s1 = smotif[0]
            res_sofar = res_sofar + ss_combi[s1][0][1]
            percent = (res_sofar / float(total_ss_res)) * 100
            print percent
            if i == 1:
                num_of_models = (len(map_route) - i) * 125
                print "num_of_models", num_of_models
                outline = 'Executable="mpirun -np 128 python stage2_mpi_run.py ' + str(
                    num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                print outline
                fout.write(outline)
                continue
            if 25.0 < percent <= 50.0:
                if i >= 2:
                    num_of_models = (len(map_route) - i) * 125
                    print "num_of_models", num_of_models
                    outline = 'Executable="mpirun -np 128 python stage2_mpi_run.py ' + str(
                        num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                    print outline
                    fout.write(outline)
                else:
                    num_of_models = (len(map_route) - i) * 125
                    print "num_of_models", num_of_models
                    outline = 'Executable="mpirun -np 128 python stage2_mpi_run.py ' + str(
                        num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                    print outline
                    fout.write(outline)
            if 50.0 < percent <= 75.0:
                num_of_models = (len(map_route) - i) * 125
                print "num_of_models", num_of_models
                outline = 'Executable="mpirun -np 128 python stage3x_mpi_run.py ' + str(
                    num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                print outline
                fout.write(outline)
            if percent > 75.0:
                num_of_models = (len(map_route) - i) * 125
                print "num_of_models", num_of_models
                outline = 'Executable="mpirun -np 128 python stage4x_mpi_run.py ' + str(
                    num_of_models) + '"\nrun="$Executable"\necho $run\n$run\n'
                print outline
                fout.write(outline)
fout.close()
run = 'cat submit0.sh concat.txt > submitX.sh'
os.system(run)
