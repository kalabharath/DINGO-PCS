import  sys, os, glob, cPickle
import utility.io_util as io
import filters.rmsd.qcp as rms
import  numpy as np

def conc_sss(ss1,ss2):
    #x,y,z, atom_type, res_no, res = [], [], [], [], [], []
    x = np.concatenate((ss1[0], ss2[0]))
    y = np.concatenate((ss1[1], ss2[1]))
    z = np.concatenate((ss1[2], ss2[2]))
    atom_type =ss1[3]+ss2[3]
    res_no= ss1[4]+ss2[4]
    res = ss1[5]+ss2[5]
    # x = ss1[0]+ss2[0]
    # y = ss1[1]+ss2[1]
    # z = ss1[2]+ss2[2]
    #return [x,y,z]
    return [x,y,z, atom_type, res_no, res]
def translateSSE(cm, sse):
    for i in range(0,len(sse)):
        sse[i][-3] -= cm[0]
        sse[i][-2] -= cm[1]
        sse[i][-1] -= cm[2]

        sse[i][-3] = round(sse[i][-3], 3)
        sse[i][-2] = round(sse[i][-2], 3)
        sse[i][-1] = round(sse[i][-1], 3)

        #print sse[i]

    return sse

def dumpPDBCoo2(coo_array, i):
    outfile = open(str(i)+"_.pdb", 'w')
    for i in range(0, len(coo_array[0])):
        x = coo_array[0][i]
        y = coo_array[1][i]
        z = coo_array[2][i]
        atom = coo_array[3][i]
        res_no = coo_array[4][i]
        res = coo_array[5][i]

        pdb_line = "%-6s%5d  %-2s%5s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"\
                   %('ATOM',i+1,atom,res,'A',res_no," ",x, y, z,1.0,30.0,' ',' \n')
        outfile.write(pdb_line)
        #print pdb_line
    print 'TER'
    outfile.close()
    return True

smotif_db = glob.glob("/home/kalabharath/zinr/final_smotif_database/*.db")
print len(smotif_db)

for i in range(0, len(smotif_db)):
#for i in range(0, 1):

    temp_smotif=[]
    smotif = io.readPickle(smotif_db[i])
    print smotif_db[i]
    print len(smotif)
    print smotif_db[i][45:]
    for j in range(0, len(smotif)):
    #for j in range(0, 1):
        #print smotif[j][0]
        #print smotif[j][0][0] #['1pp9P00', '172', '203', '222', '245']
        #print smotif[j][0][1] #172- [203, 'THR', 'H', 49.366, 87.846, 102.07]]
        #print smotif[j][0][2] #222- [245, 'PHE', 'H', 40.016, 57.964, 78.166]]

        ss1 = rms.getcoo(smotif[j][0][1]) #x,y,z, atom_type, res_no, res = [], [], [], [], [], []
        ss2 = rms.getcoo(smotif[j][0][2]) #return [x,y,z, atom_type, res_no, res]

        #dumpPDBCoo2(ss1, 5)
        #dumpPDBCoo2(ss2, 6)

        #concatenate coordinate arrys of both SSEs to identify their center of mass
        ss1_plus_ss2 = conc_sss(ss1,ss2)

        cen_ss1, cm_ss1 = rms.centerCoo(ss1_plus_ss2)
        #dumpPDBCoo2(cen_ss1, 7)

        # Translate the original SSEs to the new center of mass
        tr_ss1 = translateSSE(cm_ss1,smotif[j][0][1])
        tr_ss2 = translateSSE(cm_ss1,smotif[j][0][2])
        smotif_def = smotif[j][0][0]
        #add the modified Smotifs to temp array
        temp_smotif.append([smotif_def, tr_ss1, tr_ss2])
    io.dumpPickle('./smotif_cen_db/'+smotif_db[i][45:], temp_smotif)
