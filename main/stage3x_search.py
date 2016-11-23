#!/usr/bin/env python

"""
Project_Name: main, File_name: stage3x_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 3x in parallel
"""
import time

import filters.constraints.looplengthConstraint as llc
import filters.pcs.pcsfilter as Pfilter
import filters.rdc.rdcfilter as Rfilter
import filters.noe.noefilter as Nfilter
import filters.rmsd.qcp as qcp
import filters.sequence.sequence_similarity as Sfilter
import utility.io_util as io
import utility.smotif_util as sm
import utility.stage2_util as uts2


def getfromDB(previous_smotif, current_ss, direction, database_cutoff):
    # print "previous_smotif: ", previous_smotif

    searched_smotifs = []
    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            searched_smotifs = entry[-1]

    # ['smotif_def', [['helix', 6, 7, 5, 145, 150], ['helix', 23, 5, 1, 156, 178], ['strand', 5, 7, 8, 133, 137]]]

    if direction == 'left':
        previous_ss = searched_smotifs[0]
    else:
        previous_ss = searched_smotifs[-1]

    # print "previous_ss: ", previous_ss
    # print "current_ss : ", current_ss

    if direction == 'left':  # double check this implementation

        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def, database_cutoff), smotif_def


def orderSSE(previous_smotif, current_sse, direction):
    """

    :param previous_smotif:
    :param current_sse:
    :return:
    """

    for entry in previous_smotif:
        # ['qcp_rmsd', transformed_coos, sse_ordered, rmsd]
        if 'qcp_rmsd' == entry[0]:
            previous_sse = entry[2]

            if direction == 'left':
                previous_sse.insert(0,current_sse)
            else:
                previous_sse.append(current_sse)
            return previous_sse


def SmotifSearch(index_array):
    """
    Main()
    :param index_array:
    :return:
    """


    preSSE = uts2.getPreviousSmotif(index_array[0])
    current_ss, direction = uts2.getSS2(index_array[1])
    print current_ss, direction

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts', 'natives']
    csmotif_data, smotif_def = getfromDB(preSSE, current_ss, direction, exp_data['database_cutoff'])

    if not csmotif_data:
        # If the smotif library doesn't exist.
        # Terminate further execution by return value.
        return True

    sse_ordered = orderSSE(preSSE, current_ss, direction)
    #print sse_ordered
    dump_log = []
    no_clashes = False

    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************

    for i in range(0, len(csmotif_data)):
        # ************************************************
        # Applying different filters for the Smotif assembly
        # ************************************************

        # Exclude the natives, if present.
        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = csmotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            #if pdbid not in ['2z2i']:
            if pdbid in natives:
                # Stop further execution, but resume iteration
                continue

        # ************************************************
        # RMSD filter using QCP method
        # quickly filters non-overlapping Smotifs
        # ************************************************

        rmsd, transformed_coos = qcp.rmsdQCP3(preSSE, csmotif_data[i], direction)


        if rmsd <= exp_data['rmsd_cutoff'][2]:

            # Loop constraint restricts the overlapping smotifs is not drifted far away.
            loop_constraint = llc.loopConstraint(transformed_coos, sse_ordered, direction, smotif_def)
            if loop_constraint:
                # Check whether the SSEs with in the assembled smotifs are clashing to one another
                no_clashes = qcp.clahses(transformed_coos, exp_data['clash_distance'])
            else:
                no_clashes = False


        if rmsd <= exp_data['rmsd_cutoff'][2] and no_clashes:
            # Prepare temp log array to save data at the end
            tlog, noe_fmeasure, pcs_tensor_fits, rdc_tensor_fits = [], [], [], []
            tlog.append(['smotif', csmotif_data[i]])
            tlog.append(['smotif_def', sse_ordered])
            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])

            cathcodes = sm.orderCATH(preSSE, csmotif_data[i][0], direction)
            #print cathcodes
            tlog.append(['cathcodes', cathcodes])

            # ************************************************
            # Sequence filter
            # Aligns the smotif seq to target seq and calculates
            # sequence identity and the alignment score
            # ************************************************

            csse_seq, seq_identity, blosum62_score = Sfilter.S2SequenceSimilarity(current_ss, csmotif_data[i],
                                                                                  direction, exp_data)

            concat_seq = sm.orderSeq(preSSE, csse_seq, direction)

            tlog.append(['seq_filter', concat_seq, csse_seq, seq_identity, blosum62_score])

            # ************************************************
            # Pseudocontact Shift filter
            # uses experimental PCS data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'pcs_data' in exp_data_types:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data, stage = 3)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            # ************************************************
            # Ambiguous NOE score filter
            # uses experimental ambiguous noe data to filter Smotifs
            # scoring based on f-measure?
            # ************************************************

            if 'noe_data' in exp_data_types:
                noe_fmeasure = Nfilter.s2NOEfit(transformed_coos, sse_ordered, exp_data)
                tlog.append(['NOE_filter', noe_fmeasure])

            # ************************************************
            # Residual dipolar coupling filter
            # uses experimental RDC data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'rdc_data' in exp_data_types:
                if noe_fmeasure and noe_fmeasure > 0.5:
                    rdc_tensor_fits = Rfilter.RDCAxRhFit2(transformed_coos, sse_ordered, exp_data, stage=2)
                    tlog.append(['RDC_filter', rdc_tensor_fits])


            if pcs_tensor_fits or rdc_tensor_fits:
                #print csmotif_data[i][0],"seq_id", seq_identity, "rmsd=", rmsd, cathcodes
                dump_log.append(tlog)

    if len(dump_log) > 0:
        print "num of hits", len(dump_log)
        io.dumpPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)
    return True
