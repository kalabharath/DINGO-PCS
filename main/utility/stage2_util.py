import collections, glob, os, math
import io_util as io

def enum(*sequential, **named):
    """

    :param sequential:
    :param named:
    :return:
    """
    # fake an enumerated type in Python

    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def checkFile(i):
    regex = str(i) + "_*.pickle"
    file_list = glob.glob(regex)
    if len(file_list) > 0:
        return True
    else:
        return False


def getNextSmotif(map_route):
    # [[6, 7, 'start'], [5, 6, 'left'], [4, 5, 'left'], [3, 4, 'left'], [2, 3, 'left'], [1, 2, 'left'], [0, 1, 'left']
    for i in range(0, len(map_route)):
        if not checkFile(i):
            return i, map_route[i]


def scoreCombination(score_list):
    """

    :param score_list:
    :return:
    """

    import itertools
    min_score = 999
    combi_list = list(itertools.combinations(score_list, 2))

    for combi in combi_list:
        c1 = combi[0]
        c2 = combi[1]
        if c1 and c2:
            if c1 + c2 < min_score:
                min_score = c1 + c2
    return min_score


def scoreCombination4t(score_list):
    """

    :param score_list:
    :return:
    """

    import itertools
    min_score = 999
    combi_list = list(itertools.combinations(score_list, 3))

    for combi in combi_list:
        c1 = combi[0]
        c2 = combi[1]
        c3 = combi[2]
        if c1 and c2 and c3:
            if c1 + c2 + c3 < min_score:
                min_score = c1 + c2 + c3
    return min_score


def getNchiSum(pcs_filter, stage):
    """

    :param pcs_filter:
    :return:
    """
    snchi = 999.999
    tensors = pcs_filter[1]

    if len(tensors) == 1:
        # Discourage single tag scoring by returning high score
        return 999.999  # discourage double tag score only for 4 tags

    if len(tensors) == 2 and stage == 2:  # stage 2
        snchi = 0
        for tensor in tensors:
            nchi = tensor[1]
            snchi += nchi
            # return 999.999 #discourage double tag score only for 4 tags

    if len(tensors) == 3 and stage <= 3:  # stage 2 & 3
        # Scoring three tags, get lowest Nchi for 2
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination(score_list)

    if len(tensors) >= 4 and stage <= 3:  # stage 2,3 & 4
        # For 4 tags, get lowest Nchi for 3
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination4t(score_list)
        snchi /= 100.0  # artificially increase the priority

    if len(tensors) >= 4 and stage == 4:  # stage 4
        # For 4 tags, get lowest Nchi for 3
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = scoreCombination4t(score_list)
    """
    if len(tensors) ==3 and stage == 4:
        score_list = []
        for tensor in tensors:
            score_list.append(tensor[1])
        snchi = score_list[0]+score_list[1]+score_list[2]
    """
    return snchi

def rdcSumChi(rdc_data, stage):
    snchi = 999.999
    tensors = rdc_data[1]
    if len(tensors) == 1:
        for tensor in tensors:
            return tensor[0]

    if len(tensors) == 2 and stage == 2:
        snchi = 0
        for tensor in tensors:
            nchi = tensor[0]
            snchi += nchi
    if len(tensors) == 2:
        snchi = 0
        for tensor in tensors:
            nchi = tensor[0]
            snchi += nchi

    return snchi
def makeTopPickle(previous_smotif_index, num_hits, stage):
    """
    Concatenate data from all of the threads, organize, remove redundancies, rank
     and extract top hits as defined
    :param previous_smotif_index:
    :param num_hits:
    :param stage:
    :return:
    """
    hits = []
    regex = str(previous_smotif_index) + "_*_*.pickle"
    file_list = glob.glob(regex)
    for f in file_list:
        t_hits = io.readPickle(f)
        for t_hit in t_hits:
            hits.append(t_hit)
    """
    identifiers: smotif, smotif_def, seq_filter, contacts_filter, PCS_filter, qcp_rmsd, Evofilter
                 RDC_filter, NOE_filter
    """

    new_dict = collections.defaultdict(list)
    pcs_filter = False
    contact_filter = False
    rdc_filter = False
    noe_filter = False
    for hit in hits:
        # thread_data contains data from each search and filter thread.
        for data_filter in hit:
            if data_filter[0] == 'PCS_filter':
                pcs_filter = True
                pcs_data = data_filter
                Nchi = getNchiSum(pcs_data, stage)
                # new_dict.setdefault(Nchi, []).append(entry)
                new_dict[Nchi].append(hit)

            if data_filter[0] == 'Evofilter':
                contact_filter = True
                new_dict[data_filter[1]].append(hit)

            if data_filter[0] == 'RDC_filter':
                rdc_filter = True
                rdc_data = data_filter
                Nchi = rdcSumChi(rdc_data, stage)
                for filter in hit:
                    if filter[0] == 'NOE_filter':
                        noe_filter = True
                        noe_fmeasure = filter[1]
                        Nchi = Nchi / math.pow(10, noe_fmeasure * 10)
                        new_dict[Nchi].append(hit)
                if not noe_filter:
                    new_dict[Nchi].append(hit)

    # ************************************************
    # Exclude the redundant entries and rank top hits
    # ************************************************

    keys = new_dict.keys()
    keys.sort()
    if contact_filter and not pcs_filter:
        # Contact filter data should be as high as possible
        keys.reverse()

    # Exclude the redundant data.

    # non_redundant = {}
    non_redundant = collections.defaultdict(list)
    seqs = []
    smotif_seq = ''
    Nchi = 0.0
    for i in range(0, len(keys)):
        entries = new_dict[keys[i]]
        for entry in entries:
            for ent in entry:
                if ent[0] == 'smotif':
                    name = ent[1][0]
                if ent[0] == 'seq_filter':
                    seq_filter = ent
                    smotif_seq = seq_filter[1]
                if ent[0] == 'PCS_filter':
                    pcs_data = ent
                    Nchi = getNchiSum(pcs_data, stage)
                if ent[0] == 'Evofilter':
                    Nchi = ent[1]
                if ent[0] == 'RDC_filter':
                    rdc_data = ent
                    Nchi = rdcSumChi(rdc_data, stage)
                    if noe_filter:
                        for ent in entry:
                            if ent[0] == 'NOE_filter':
                                noe_fmeasure = ent[1]
                                Nchi = Nchi /math.pow(10, noe_fmeasure * 10)
                    else:
                        Nchi = rdcSumChi(rdc_data, stage)

            if smotif_seq not in seqs:
                seqs.append(smotif_seq)
                # non_redundant.setdefault(Nchi, []).append(entry)
                non_redundant[Nchi].append(entry)

    # Rank top hits and dump the data
    keys = non_redundant.keys()
    keys.sort()
    if contact_filter and not pcs_filter:
        keys.reverse()
    dump_pickle = []
    print "Dumping data to disk"
    count_top_hits = 0
    while (True):
        for key in keys:
            if key ==  999.999:
                # Do not work on these entries
                continue
            entries = non_redundant[key]
            for entry in entries:
                dump_pickle.append(entry)
                print "final sele", entry[0][1][0][0], key
                count_top_hits += 1
            if count_top_hits >= num_hits:
                break
        if count_top_hits >= num_hits:
            break
        else:
            print "could only extract ", count_top_hits
            break

    io.dumpPickle(str(previous_smotif_index) + "_tophits.pickle", dump_pickle)
    print "actual number in top hits ", len(dump_pickle)
    return range(count_top_hits)

def getRunSeq(num_hits, stage):
    """
    generate run seq, a seq list of pairs of
    indexes of profiles for job scheduling
    """

    ss_profiles = io.readPickle("ss_profiles.pickle")
    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    try:
        next_index, next_smotif = getNextSmotif(map_route)
        print next_index, next_smotif
    except TypeError:
        return [999], 999

    direction = next_smotif[-1]
    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]
    # get and make a list of top 10(n) of the previous run
    top_hits = makeTopPickle(next_index - 1, num_hits, stage)  # send the previous Smotif index

    # delete two stages down pickled files
    check_pickle = str(next_index - 2) + str("_*_*.pickle")
    file_list = glob.glob(check_pickle)

    if len(file_list) > 10:
        remove = "rm "+check_pickle
        os.system(remove)

    if top_hits:
        run_seq = []
        for i in range(len(top_hits)):
            for j in range(len(next_ss_list)):
                run_seq.append([i, j])
        return run_seq, next_index


def getPreviousSmotif(index):

    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    next_index, next_smotif = getNextSmotif(map_route)
    top_hits = io.readPickle(str(next_index - 1) + "_tophits.pickle")  # Read in previous index hits
    # print len(top_hits)
    return top_hits[index]


def getSS2(index):


    if os.path.isfile("contacts_route.pickle"):
        map_route = io.readPickle("contacts_route.pickle")
    elif os.path.isfile("pcs_route.pickle"):
        map_route = io.readPickle("pcs_route.pickle")
    elif os.path.isfile("rdc_route.pickle"):
        map_route = io.readPickle("rdc_route.pickle")

    ss_profiles = io.readPickle("ss_profiles.pickle")
    
    next_index, next_smotif = getNextSmotif(map_route)
    direction = next_smotif[-1]

    if direction == 'left':
        next_ss_list = ss_profiles[next_smotif[0]]
    else:
        next_ss_list = ss_profiles[next_smotif[1]]

    return next_ss_list[index], direction


def rename_pickle(index):
    import glob, os

    file_list = glob.glob("tx_*.pickle")
    for file in file_list:
        mv_cmd = "mv " + file + " " + str(index) + file[2:]
        os.system(mv_cmd)
    return True
