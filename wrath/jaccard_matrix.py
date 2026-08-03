from multiprocessing import Process
from multiprocessing.queues import SimpleQueue
import sys
from threading import Thread
from time import sleep
import time

import numpy as np
import polars as pl
import pysam

start_time = time.time()

resultsReceived = 0
windowQueued = 0
resultsWritten = 0


def freqs_wrapper(inQueue, resultQueue, number_win, inFile, windowfile):
    while True:
        windowNumber, windowLine = inQueue.get()
        if windowNumber == -1:
            resultQueue.put((-1, None,))
            break
        array_i = np.zeros((number_win))
        array_u = np.ones((number_win))
        bedfile1 = inFile.fetch(windowLine[0], windowLine[1], windowLine[2], parser=pysam.asBed(), multiple_iterators=True)
        barcodes1 = [rowbed1.name for rowbed1 in bedfile1]
        for index2, row2 in enumerate(windowfile[windowNumber:].iter_rows(), start=windowNumber):
            bedfile2 = inFile.fetch(row2[0], row2[1], row2[2], parser=pysam.asBed(), multiple_iterators=True)
            barcodes2 = [rowbed2.name for rowbed2 in bedfile2]
            intersect = np.intersect1d(barcodes1, barcodes2)
            union = np.union1d(barcodes1, barcodes2)
            array_i[index2] = intersect.size
            if union.size > 1:
                array_u[index2] = union.size
            elif union.size < 1:
                array_u[index2] = 1
                array_i[index2] = 0
            array_u[index2] = union.size
        outArray = np.divide(array_i, array_u)
        resultQueue.put((windowNumber, outArray,))


def sorter(doneQueue, writeQueue, verbose, nWorkerThreads):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    threadsComplete = 0
    while True:
        windowNumber, results = doneQueue.get()
        if windowNumber == -1:
            threadsComplete += 1
        if threadsComplete == nWorkerThreads:
            writeQueue.put((-1, None,))
            break
        resultsReceived += 1
        if verbose:
            sys.stderr.write("Sorter received window {}\n".format(windowNumber))
        if windowNumber == expect:
            writeQueue.put((windowNumber, results))
            if verbose:
                sys.stderr.write("window {} sent to writer\n".format(windowNumber))
            expect += 1
            while True:
                try:
                    results = sortBuffer.pop(str(expect))
                    writeQueue.put((expect, results))
                    if verbose:
                        sys.stderr.write("window {} sent to writer\n".format(expect))
                    expect += 1
                except KeyError:
                    break
        else:
            sortBuffer[str(windowNumber)] = results


def writer(writeQueue, out, verbose):
    global resultsWritten
    while True:
        windowNumber, results = writeQueue.get()
        if windowNumber == -1:
            break
        if verbose:
            sys.stderr.write("Writer received window {}\n".format(windowNumber))
        np.savetxt(out, results, fmt='%.10f', newline=',')
        out.write("\n")
        resultsWritten += 1


def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("{} windows queued | {} windows analysed | {} windows written\n".format(
            windowQueued, resultsReceived, resultsWritten))

def jaccard_matrix(winfile: str, barcodefile: str, outfile: str | None = None, threads: int = 1, verbose: bool = False):
    global windowQueued, resultsReceived, resultsWritten

    outfile = open(outfile, "wt") if outfile else sys.stdout

    windowfile = pl.read_csv(winfile, separator='\t', has_header=False)
    num_win = windowfile.height

    tbx = pysam.TabixFile(barcodefile)

    windowQueued = 0
    resultsReceived = 0
    resultsWritten = 0

    inQueue = SimpleQueue()
    resultQueue = SimpleQueue()
    writeQueue = SimpleQueue()

    workerThreads = []
    sys.stderr.write("\nStarting {} worker threads\n".format(threads))
    for _ in range(threads):
        workerThread = Process(target=freqs_wrapper, args=(inQueue, resultQueue, num_win, tbx, windowfile,))
        workerThread.daemon = True
        workerThread.start()
        workerThreads.append(workerThread)

    sorterThread = Thread(target=sorter, args=(resultQueue, writeQueue, verbose, threads,))
    sorterThread.daemon = True
    sorterThread.start()

    writerThread = Thread(target=writer, args=(writeQueue, outfile, verbose,))
    writerThread.daemon = True
    writerThread.start()

    checkerThread = Thread(target=checkStats)
    checkerThread.daemon = True
    checkerThread.start()

    for windowIdx, windowLine in enumerate(windowfile.iter_rows()):
        inQueue.put((windowQueued, windowLine))
        windowQueued += 1

    for _ in range(threads):
        inQueue.put((-1, None,))

    sys.stderr.write("\nClosing worker threads\n")
    for workerThread in workerThreads:
        workerThread.join()

    sorterThread.join()
    writerThread.join()

    outfile.close()
    sys.stderr.write(f"\nDone in {time.time() - start_time}\n")