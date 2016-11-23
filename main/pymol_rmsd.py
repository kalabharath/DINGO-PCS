import pymol, sys, os, glob
#__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
pymol.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
pymol.finish_launching()
# Disable Pymol's output

#pymol.cmd.feedback("disable","all","actions")
#pymol.cmd.feedback("disable","all","results")
#pymol.cmd.feedback("disable","all","errors")
#pymol.cmd.feedback("disable","all","warnings")

def get_seqs(filename):

    with open(filename) as fin:
        lines = fin.readlines()

    for line in lines:
        if line[2:12] == 'seq_filter':
            #print line
            pass
        if line[2:12] == 'PCS_filter':
            print line
    return True



target2 = glob.glob("./setup/ideal*.pdb")

target2 = target2[0]

pdbfiles = []

seqs = get_seqs('final.log')


for i in range(0,50):
    pdbfiles.append("model_"+str(i)+".pdb")

for pdbfile in pdbfiles:
    target1 = pdbfile

    pymol.cmd.load(target1, 'model0')
    pymol.cmd.load(target2, 'ideal')


    #rmsd = pymol.cmd.align('model', 'ideal')

    rmsd = pymol.cmd.align('model0', 'ideal', cycles=10, transform=0, object = 'aln')
    #print rmsd
    print rmsd, target1, target2
    pymol.cmd.refresh()
    pymol.cmd.reinitialize()
