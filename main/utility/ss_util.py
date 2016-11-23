def genSSDef(ss_seq):
    """
	returns an array of smotifs derived from the ss_seq with extended
	functionality ie information about left_loop length
	return [[ss_type,len_ss,l_loop,start,end]['strand', 5, 4, 5, 9]]
	:param ss_seq:
	:return:
	"""
    helix, strand, loop, l_loop = 1, 1, 1, 0
    return_array = []
    seq = []
    for i in range(0, len(ss_seq)):
        if i == 0:
            continue
        else:
            if (ss_seq[i] == 'C') or (ss_seq[i] == 'L'):
                l_loop += 1
            if ss_seq[i] == ss_seq[i - 1]:
                if ss_seq[i] == 'H':
                    seq.append(i)
                    helix += 1
                if ss_seq[i] == 'E':
                    seq.append(i)
                    strand += 1
                if (ss_seq[i] == 'C') or (ss_seq[i] == 'L'):
                    loop += 1
            else:
                if (helix > 4) and ( (strand == 1) ):
                    end = i
                    return_array.append(['helix', helix, l_loop, seq[0], end])
                    seq = []
                    l_loop = 0
                elif (strand > 4) and ( helix == 1):
                    end = i
                    return_array.append(['strand', strand, l_loop, seq[0], end])
                    seq = []
                    l_loop = 0
                elif (loop > 4) and ( (strand == 1) and (helix == 1)):
                    end = i
                    seq = []
                helix, strand, loop = 1, 1, 1
                seq = []
    return return_array, l_loop + 1


def genSSDef_BACK(ss_seq):
    """
    BACKED up comy that uses a min of 4 residues per SSE


	returns an array of smotifs derived from the ss_seq with extended
	functionality ie information about left_loop length
	return [[ss_type,len_ss,l_loop,start,end]['strand', 5, 4, 5, 9]]
	:param ss_seq:
	:return:
	"""
    helix, strand, loop, l_loop = 1, 1, 1, 0
    return_array = []
    seq = []
    for i in range(0, len(ss_seq)):
        if i == 0:
            continue
        else:
            if (ss_seq[i] == 'C') or (ss_seq[i] == 'L'):
                l_loop += 1
            if ss_seq[i] == ss_seq[i - 1]:
                if ss_seq[i] == 'H':
                    seq.append(i)
                    helix += 1
                if ss_seq[i] == 'E':
                    seq.append(i)
                    strand += 1
                if (ss_seq[i] == 'C') or (ss_seq[i] == 'L'):
                    loop += 1
            else:
                if (helix > 3) and ( (strand == 1) ):
                    end = i
                    return_array.append(['helix', helix, l_loop, seq[0], end])
                    seq = []
                    l_loop = 0
                elif (strand > 3) and ( helix == 1):
                    end = i
                    return_array.append(['strand', strand, l_loop, seq[0], end])
                    seq = []
                    l_loop = 0
                elif (loop > 3) and ( (strand == 1) and (helix == 1)):
                    end = i
                    seq = []
                helix, strand, loop = 1, 1, 1
                seq = []
    return return_array, l_loop + 1

def check_profile(ss_type, tlen_ss, tl_loop, tr_loop, tstart, tend, X, Y):
    """
	the def checks whether each SS is of len > 5 residues
	and whether an smotif is possible with at least 1 residue inbetween
	X denotes either left or right of the SS.
	Y denotes the possible extension of +/- 2.
	return [ss_type,len_ss,l_loop,r_loop,start,end]
	"""
    if X:
        if Y == 0:
            tlen_ss -= 2
            tstart += 1
            tend = tend - 1
            tr_loop += 1
            tl_loop += 1
            if (tlen_ss > 4) and (tl_loop >= 0 and tr_loop >= 0 ):
                if (tend - tstart + 1) == tlen_ss:
                    return [ss_type, tlen_ss, tl_loop, tr_loop, tstart, tend]
            else:
                return False
        else:
            tlen_ss = tlen_ss + Y
            # tstart
            tend = tend + Y
            tr_loop = tr_loop + (-1 * Y)
            if (tlen_ss > 4) and (tl_loop >= 0 and tr_loop >= 0 ):
                if (tend - tstart + 1) == tlen_ss:
                    return [ss_type, tlen_ss, tl_loop, tr_loop, tstart, tend]
            else:
                return False
    else:
        if Y == 0:
            tlen_ss += 2
            tstart -= 1
            tend = tend + 1
            tr_loop -= 1
            tl_loop -= 1
            if (tlen_ss > 4) and (tl_loop >= 0 and tr_loop >= 0 ):
                if (tend - tstart + 1) == tlen_ss:
                    return [ss_type, tlen_ss, tl_loop, tr_loop, tstart, tend]
            else:
                return False
        else:
            tlen_ss = tlen_ss + Y
            tstart = tstart + (-1 * Y)
            # tend
            tl_loop = tl_loop + (-1 * Y)
            if (tlen_ss > 4) and (tl_loop >= 0 and tr_loop >= 0 ):
                if (tend - tstart + 1) == tlen_ss:
                    return [ss_type, tlen_ss, tl_loop, tr_loop, tstart, tend]
            else:
                return False


def genSSCombinations(ss_seq):
    """
	:param ss_seq:
	:return exn_ss, ss_combi
	"""
    exn_ss, last_loop = genSSDef(ss_seq)
    ss_combi = {}
    for k in range(0, len(exn_ss)):
        if k != len(exn_ss) - 1:
            ss_type = exn_ss[k][0]
            len_ss = exn_ss[k][1]
            l_loop = exn_ss[k][2]
            r_loop = exn_ss[k + 1][2]
            start = exn_ss[k][3]
            end = exn_ss[k][4]
            arr_0 = [ss_type, len_ss, l_loop, r_loop, start, end]
            ss_combi.setdefault(k, []).append(arr_0)
            for i in range(0, 2):
                if i == 0:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
                else:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
        else:
            ss_type = exn_ss[k][0]
            len_ss = exn_ss[k][1]
            l_loop = exn_ss[k][2]
            r_loop = last_loop
            start = exn_ss[k][3]
            end = exn_ss[k][4]
            arr_0 = [ss_type, len_ss, l_loop, r_loop, start, end]
            ss_combi.setdefault(k, []).append(arr_0)
            for i in range(0, 2):
                if i == 0:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
                else:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
    return exn_ss, ss_combi
def genSSCombinations2(ss_seq):
    """
	:param ss_seq:
	:return exn_ss, ss_combi
	"""
    exn_ss, last_loop = genSSDef(ss_seq)
    ss_combi = {}
    for k in range(0, len(exn_ss)):
        if k != len(exn_ss) - 1:
            ss_type = exn_ss[k][0]
            len_ss = exn_ss[k][1]
            l_loop = exn_ss[k][2]
            r_loop = exn_ss[k + 1][2]
            start = exn_ss[k][3]
            end = exn_ss[k][4]
            arr_0 = [ss_type, len_ss, l_loop, r_loop, start, end]
            ss_combi.setdefault(k, []).append(arr_0)
            for i in range(0, 2):
                if i == 0:
                    for j in range(0, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
                else:
                    for j in range(0, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
        else:
            ss_type = exn_ss[k][0]
            len_ss = exn_ss[k][1]
            l_loop = exn_ss[k][2]
            r_loop = last_loop
            start = exn_ss[k][3]
            end = exn_ss[k][4]
            arr_0 = [ss_type, len_ss, l_loop, r_loop, start, end]
            ss_combi.setdefault(k, []).append(arr_0)
            for i in range(0, 2):
                if i == 0:
                    for j in range(0, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
                else:
                    for j in range(0, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
    return exn_ss, ss_combi

def genSSCombinationsshort(ss_seq):
    """
	:param ss_seq:
	:return exn_ss, ss_combi
	"""
    exn_ss, last_loop = genSSDef(ss_seq)
    print exn_ss
    ss_combi = {}
    for k in range(0, len(exn_ss)):
        if k != len(exn_ss) - 1:
            ss_type = exn_ss[k][0]
            len_ss = exn_ss[k][1]
            l_loop = exn_ss[k][2]
            r_loop = exn_ss[k + 1][2]
            start = exn_ss[k][3]
            end = exn_ss[k][4]
            arr_0 = [ss_type, len_ss, l_loop, r_loop, start, end]
            ss_combi.setdefault(k, []).append(arr_0)
            for i in range(0, 2):
                if i == 0:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
                else:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
        else:
            ss_type = exn_ss[k][0]
            len_ss = exn_ss[k][1]
            l_loop = exn_ss[k][2]
            r_loop = last_loop
            start = exn_ss[k][3]
            end = exn_ss[k][4]
            arr_0 = [ss_type, len_ss, l_loop, r_loop, start, end]
            ss_combi.setdefault(k, []).append(arr_0)
            for i in range(0, 2):
                if i == 0:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
                else:
                    for j in range(-2, 3):
                        arr = check_profile(ss_type, len_ss, l_loop, r_loop, start, end, i, j)
                        if arr:
                            ss_combi.setdefault(k, []).append(arr)
    return exn_ss, ss_combi
