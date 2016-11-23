from filters.contacts.contacts_filter import get_distance
from filters.contacts.evfoldContacts import calcFmeasure
from filters.constraints.looplengthConstraint import get_dist

def calcPrecision(gbar, ganoe):
    tp, fp, fn = 0.0, 0.0, 0.0

    for entry in ganoe:
        if entry in gbar:
            tp += 1
        else:
            fn += 1
    for entry in gbar:
        if entry in ganoe:
            pass
        else:
            fp += 1
    if tp:
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        return precision
        #return (2 * precision * recall) / (precision + recall)
    else:
        return False



def s1NOEfit(s1_def, s2_def, smotif, exp_data):


    noe_cutoff = 5.0
    noe_matrix = exp_data['noe_data']
    ss1_list = range(s1_def[4], s1_def[5] + 1)
    ss2_list = range(s2_def[4], s2_def[5] + 1)

    smotif_ss1 = range(int(smotif[0][1]), int(smotif[0][2]) + 1)
    smotif_ss2 = range(int(smotif[0][3]), int(smotif[0][4]) + 1)

    noes_found =[]
    for res in smotif_ss1:
        for entry1 in smotif[1]:
            if entry1[2] == 'H' and entry1[0] == res:
                coo1 = [entry1[3], entry1[4], entry1[5]]
                for entry2 in smotif[2]:
                    if entry2[2] == 'H':
                        coo2 = [entry2[3], entry2[4], entry2[5]]
                        dist = get_distance(coo1, coo2)
                        if dist < noe_cutoff:
                            noes_found.append(
                                (ss1_list[smotif_ss1.index(entry1[0])], ss2_list[smotif_ss2.index(entry2[0])]))

    if len(noes_found) == 0:
        return []
    noes_total = []
    for entry1 in ss1_list:
        for entry2 in ss2_list:
            if noe_matrix[entry1, entry2]:
                noes_total.append((entry1, entry2))

    fmeasure = calcFmeasure(noes_found, noes_total)
    #precision = calcPrecision(noes_found, total_noes)
    #print len(total_noes), total_noes
    #print len(noes_found), noes_found

    return fmeasure

def getNHandresi(frag):
    x, y, z = [], [], []
    resi = []
    for i in range(0, len(frag[0])):
        if frag[3][i] == 'H':
            x.append(frag[0][i])
            y.append(frag[1][i])
            z.append(frag[2][i])
            resi.append(frag[4][i])
    return resi, [x, y, z]

def s2NOEfit(transformed_coors, native_sse_order, exp_data):

    noe_cutoff = 5.0
    import copy
    sse_coors = copy.deepcopy(transformed_coors)
    noe_matrix = exp_data['noe_data']

    noes_found = []
    noes_total = []
    for i in range(0, len(sse_coors) - 1):
        res_c, ca1 = getNHandresi(sse_coors[i])
        res_n, ca2 = getNHandresi(sse_coors[i + 1])
        ss1_list = range(native_sse_order[i][4], native_sse_order[i][5] + 1)
        ss2_list = range(native_sse_order[i + 1][4], native_sse_order[i + 1][5] + 1)
        try:
            for res1 in ss1_list:
                # print ss1_list, res1, res_c, res_c[ss1_list.index(res1)]
                # print ss1_list.index(res1)
                ca_res1 = [ca1[0][ss1_list.index(res1)], ca1[1][ss1_list.index(res1)], ca1[2][ss1_list.index(res1)]]
                for res2 in ss2_list:
                    ca_res2 = [ca2[0][ss2_list.index(res2)], ca2[1][ss2_list.index(res2)], ca2[2][ss2_list.index(res2)]]
                    dist = get_dist(ca_res1, ca_res2)

                    if noe_matrix[res1, res2] or noe_matrix[res2, res1]:
                        noes_total.append((res1, res2))

                    if dist < noe_cutoff:
                        noes_found.append((res1, res2))
        except:
            return []
    if len(noes_found) == 0:
        return []

    fmeasure = calcFmeasure(noes_found, noes_total)


    return fmeasure