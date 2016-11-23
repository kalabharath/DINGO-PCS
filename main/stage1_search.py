#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time:01:05 PM

stage 1 in parallel
"""

import filters.pcs.pcsfilter as Pfilter
import filters.rdc.rdcfilter as Rfilter
import filters.noe.noefilter as Nfilter
import filters.sequence.sequence_similarity as Sfilter
import utility.io_util as io
import utility.smotif_util as sm
import utility.stage1_util as uts1


def getSSdef(index_array):
    """

    :param index_array:
    :return:
    """
    s1_list, s2_list = uts1.getSSlist()
    return s1_list[index_array[0]], s2_list[index_array[1]]


def SmotifSearch(index_array):
    """
    Main ()
    :param index_array:
    :return:
    """

    # print index_array
    s1_def, s2_def = getSSdef(index_array)
    smotif_def = sm.getSmotif(s1_def, s2_def)
    # print s1_def, s2_def

    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    smotif_data = sm.readSmotifDatabase(smotif_def, exp_data['database_cutoff'])

    if not smotif_data:
        # If the smotif library doesn't exist, terminate further execution.
        return True

    dump_log = []
    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************

    for i in range(0, len(smotif_data)):
        # loop over for all of the entries in the smotif_db file

        # ************************************************
        # Excluding the natives
        # ************************************************

        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            tpdbid = smotif_data[i][0][0]
            pdbid = tpdbid[0:4]
            if pdbid in natives:
            #if pdbid not in ['2z2i']:
                # Stop further execution, but, iterate.
                continue

        # ************************************************
        # Applying different filters to Smotifs
        # Prepare temp log array to save data at the end
        # ************************************************

        tlog, pcs_tensor_fits, rdc_tensor_fits, noe_fmeasure = [], [], [], []
        tlog.append(['smotif', smotif_data[i]])
        tlog.append(['smotif_def', [s1_def, s2_def]])
        tlog.append(['cathcodes', [smotif_data[i][0]]])

        # ************************************************
        # Sequence filter
        # Aligns the smotif seq to target seq and calculates
        # sequence identity and the alignment score
        # ************************************************

        smotif_seq, seq_identity, blosum62_score = \
            Sfilter.SequenceSimilarity(s1_def, s2_def, smotif_data[i], exp_data)
        tlog.append(['seq_filter', smotif_seq, seq_identity, blosum62_score])

        # ************************************************
        # Pseudocontact Shift filter
        # uses experimental PCS data to filter Smotifs
        # scoring based on normalised chisqr
        # ************************************************

        if 'pcs_data' in exp_data_types:
            pcs_tensor_fits = Pfilter.PCSAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['PCS_filter', pcs_tensor_fits])


        # ************************************************
        # Ambiguous NOE score filter
        # uses experimental ambiguous noe data to filter Smotifs
        # scoring based on f-measure?
        # ************************************************

        if 'noe_data' in exp_data_types:
            noe_fmeasure = Nfilter.s1NOEfit(s1_def, s2_def, smotif_data[i], exp_data)
            tlog.append(['NOE_filter', noe_fmeasure ])



        # ************************************************
        # Residual dipolar coupling filter
        # uses experimental RDC data to filter Smotifs
        # scoring based on normalised chisqr
        # ************************************************

        if 'rdc_data' in exp_data_types:

            if noe_fmeasure and noe_fmeasure > 0.5:
                rdc_tensor_fits = Rfilter.RDCAxRhFit(s1_def, s2_def, smotif_data[i], exp_data)
                tlog.append(['RDC_filter', rdc_tensor_fits])


        # Dump the data to the disk
        if pcs_tensor_fits or rdc_tensor_fits:
            # print smotif_data[i][0][0], "seq_id", seq_identity, "i=", i, "/", len(smotif_data)
            print tpdbid, noe_fmeasure, rdc_tensor_fits
            dump_log.append(tlog)

    # Save all of the hits in pickled arrays
    if dump_log:
        print "num of hits", len(dump_log)
        io.dumpPickle('0_' + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)

    return True
