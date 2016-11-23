#!/usr/bin/env python

"""
Project_Name: main, File_name: stage_1_mpi_run.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 13/04/15 , Time:10:05 AM

Perform stage 2 in parallel
"""
import sys
sys.path.append('../../main/')
import time
from   mpi4py import  MPI

import stage2_search as S2search
import utility.stage2_util as util

# Define MPI messaage tags
tags = util.enum('READY', 'DONE', 'EXIT', 'START')

# Init MPI comm

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
status = MPI.Status()

if rank == 0:
    num_hits = int(sys.argv[1])
    print num_hits
    try:
        tasks, sse_index = util.getRunSeq(num_hits, stage=2)
    except:
        print "Couldn't extract top hits within the specified cutoffs: Exiting..."
        for i in range(0, size - 1):
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == tags.READY:
                comm.send(None, dest=source, tag=tags.EXIT)
        exit()

    if sse_index == 999 :
        # kill all slaves if there is there is EOL
        # only makes sense for self submitting jobs
        for i in range(0, size-1):
            source = status.Get_source()
            # comm.send(None, dest = source, tag = tags.EXIT)
        exit()


    stime = time.time()


    task_index =0 # control the number of processes with this index number
    finished_task = 0
    num_workers = size -1 # 1 processor is reserved for master.
    closed_workers = 0 # control the workers with no more work that can be assigned

    # print ("Master starting with {} workers".format(num_workers))
    while closed_workers < num_workers:
        # Manage/distribute all processes in this while loop
        data = comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # worker is ready, send her something to do
            if task_index < len(tasks):
                comm.send(tasks[task_index], dest = source, tag = tags.START)
                # print ("Sending task {} to worker {}".format(task_index, source))
                task_index +=1 # increment its
            else:
                # everything is done, lets grant freedom to all
                comm.send(None, dest = source, tag = tags.EXIT)
        elif tag == tags.DONE:
            # take the result from the worker
            results = data
            ctime = time.time()
            elapsed = ctime-stime
            finished_task += 1
            print "Finishing..", finished_task, "of", len(tasks), "Smotifs, Elapsed", round((elapsed)/(60), 1), "mins"
            # print ("Got data from  worker {}".format(source))
        elif tag == tags.EXIT:
            # print ("Worker {} exited".format(source))
            closed_workers += 1

    # print "All Done, Master exiting"
    util.rename_pickle(sse_index)
    exit()

# On the worker processes
else:
    # print ("I am a worker with rank {} on {}".format(rank, name))
    while True:  # initiaite infinite loop
        comm.send(None, dest= 0, tag= tags.READY)
        # tell the master that you are ready and waiting for new assignment
        task = comm.recv(source = 0, tag = MPI.ANY_SOURCE, status = status)
        tag = status.Get_tag()

        if tag == tags.START:
            # TODO this is where you actually do something
            result = S2search.SmotifSearch(task)

            comm.send(result, dest=0, tag= tags.DONE)

        elif tag == tags.EXIT:
            # break the (while) infinite loop because there is no more work that can be assigned
            break

    # Tell the master respectfully that you are exiting
    comm.send(None, dest= 0, tag= tags.EXIT)
