'''
alignment utility functions
'''
import os
import sys
import subprocess
import logging
import mmap
from operator import itemgetter
import numpy as np

## public functions ##
def create_idx(asm_fasta, index_file):
    """     make bowtie2 index
    Parameters:
    -----------
        asm_fasta           : str
        index_file          : str
    """

    # run the command.
    subprocess.call(['bowtie2-build', '-f', asm_fasta, index_file])

def create_aln(size_file, index_file, fastq_file, sam_file, ant_dir, num_cpu):
    """     make bowtie2 alignment and
    pull out multimappers/
    Parameters:
    -----------
        index_file          : str
        fastq_file          : str
    """

    # create sizes.
    sizes = dict()
    with open(size_file, "rb") as fin:
        lines = fin.readlines()
    for line in lines:
        sz, name = line.strip().split()
        sz = int(sz)
        sizes[name] = sz

    # create the annotation arrays.
    annotes = dict()
    for ref in sizes:
        annotes[ref] = np.zeros(sizes[ref], dtype=np.int)

    # create alignment command.
    cmd = ['bowtie2','--reorder', '-k', '10', '-q','-p',str(num_cpu), '-x', index_file, '-U', fastq_file]

    # call the command and pipe output.
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    # open single-sam file.
    sam_out = open(sam_file, "wb")

    # loop over each alignment.
    for status, line in _extract_sam(output):

        # dump good ones.
        if status == True:
            sam_out.write(line)
            continue

        # record location of baduns.
        tokens = line.strip().split()
        start = int(tokens[3])
        stop = start + len(tokens[9])
        rname = tokens[2]

        # annotate that.
        annotes[rname][start:stop] += 1

    # close output.
    sam_out.close()

    # serialize multimap.
    for ref in annotes:

        # create name.
        fname = '%s/%s.npy' % (ant_dir, ref)

        # look for existing.
        if os.path.isfile(fname):
            tmp = np.load(fname)
            annotes[ref] = annotes[ref] + tmp

        # serialize it.
        np.save(fname, annotes[ref])

def pair_sam(sam_in_1, sam_in_2, sam_out_1, sam_out_2, key_size):
    """ pairs SAM files """

    # memory map the SAM files1.
    fin1 = open(sam_in_1, "r+")
    fin2 = open(sam_in_2, "r+")

    map1 = mmap.mmap(fin1.fileno(), 0, access=mmap.ACCESS_COPY)
    map2 = mmap.mmap(fin2.fileno(), 0, access=mmap.ACCESS_COPY)

    # create lists from data.
    hitlist1 = list()
    hitlist2 = list()
    for p1, p2 in _sam_gen(map1, map2, key_size):
        hitlist1.append(p1)
        hitlist2.append(p2)

    # seek files bake to begining.
    map1.seek(0)
    map2.seek(0)

    # sort lists by name, reverse so we can pop from end.
    hitlist1.sort(key=itemgetter(1), reverse=True)
    hitlist2.sort(key=itemgetter(1), reverse=True)

    # open output files.
    fout1 = open(sam_out_1, "wb")
    fout2 = open(sam_out_2, "wb")

    # generator of pairs.
    for p1, p2 in _pair_gen(hitlist1, hitlist2):

        # load sam info from map.
        map1.seek(p1[0])
        map2.seek(p2[0])

        # write out info.
        fout1.write(map1.readline())
        fout2.write(map2.readline())

    # close output files.
    fout1.close()
    fout2.close()

    # close memmory mapped files.
    map1.close()
    map2.close()

    fin1.close()
    fin2.close()


def pair_sam2(sam_in_1, sam_in_2, sam_out_1, sam_out_2, key_size):
    """ pairs SAM files """

    # memory map the first SAM file.
    logging.info("opening memory map files")

    # build name arrays.
    logging.info("extracting name array 1")
    id1 = _extract_names(sam_in_1, key_size)
    logging.info("extracting name array 2")
    id2 = _extract_names(sam_in_2, key_size)

    # sort the names and copy.
    logging.info("sorting name array 1")
    srt1 = np.sort(id1[:]['name'])
    logging.info("sorting name array 2")
    srt2 = np.sort(id2[:]['name'])

    # compute unique in each pair.
    logging.info("unique 1")
    uq1 = _numpy_unique(srt1)
    logging.info("unique 2")
    uq2 = _numpy_unique(srt2)

    # compute the intersection of unique.
    logging.info("intersection")
    valid_list = np.intersect1d(uq1, uq2, assume_unique=True)

    # sanity check.
    assert len(valid_list) != 0, 'cant have no valid stuff'

    # create a set.
    logging.info("set")
    valid = set(list(valid_list))

    # write the entries.
    logging.info("writing")
    _write_valid(sam_in_1, id1, valid, sam_out_1)
    _write_valid(sam_in_2, id2, valid, sam_out_2)
    logging.info("done")


## internal functions ##

def _write_valid(sam_in_1, id1, valid, sam_out_1):
    """ writes entries from valid set"""

    # open output.
    fout1 = open(sam_out_1, "wb")
    fin1 = open(sam_in_1, "rb")

    # generator of pairs.
    idx = 0
    for line in fin1:

        # operate.
        if id1[idx]['name'] in valid:
            fout1.write(line)

        # udpate
        idx += 1

    # close em.
    fout1.close()
    fin1.close()

def _numpy_unique(srt1):
    """ return unique subset"""

    # create mask.
    good = np.zeros(srt1.shape[0], dtype=np.bool)
    good[:] = False

    # iterate over non-boundry cases.
    for i in range(1, srt1.shape[0]-1):

        # must not match its neighbors.
        if srt1[i-1] != srt1[i] and srt1[i+1] != srt1[i]:
            good[i] = True

    # check the first one.
    if srt1[0] != srt1[1]:
        good[0] = True

    # check the last one.
    if srt1[-1] != srt1[-2]:
        good[-1] = True

    # return the subset slice.
    return srt1[good]


def _extract_names(file_name, key_size):
    """ builds numpy array of name hits"""

    # count lines.
    with open(file_name, "rb") as fin:
        line_cnt1 = 0
        for line in fin:
            line_cnt1 += 1

    # allocate array.
    id1 = np.zeros(line_cnt1, dtype=np.dtype([('name','S25'),('row',np.int)]))

    # copy data into array.
    with open(file_name, "rb") as fin:

        idx = 0
        for line1 in fin:
            # operate.
            if key_size == 0:
                id1[idx]['name'] = line1.split("\t")[0]
            else:
                id1[idx]['name'] = line1.split("\t")[0][0:-key_size]
            id1[idx]['row'] = idx

            # reset.
            idx += 1

    # return the array.
    return id1

def _extract_sam(output):
    ''' extracts output form SAM'''

    # extract unique to file, save multimap annotations.
    for line in iter(output.stdout.readline,''):

        # skip header.
        if line[0] == '@': continue

        # split.
        tokens = line.strip().split()

        # check for no align.
        if tokens[2] == '*':
            continue

        # check for MAPQ > 2:
        if int(tokens[4]) < 2:
            yield False, line
        else:
            # its good, yield it.
            yield True, line


def _sam_gen(map1, map2, key_size):
    '''yields the SAM name and the line index'''

    # loop till end of file.
    line1 = map1.readline()
    line2 = map2.readline()
    pos1 = 0
    pos2 = 0
    while line1 != '' and line2 != '':

        # process it.
        tok1 = line1.strip().split()
        tok2 = line2.strip().split()

        # remove to key.
        if key_size != 0:
            key1 = tok1[0][0:-key_size]
            key2 = tok2[0][0:-key_size]
        else:
            key1 = tok1[0]
            key2 = tok2[0]

        # yield the name and line number.
        yield (pos1, key1), (pos2, key2)

        # update info.
        pos1 += len(line1)
        pos2 += len(line2)
        line1 = map1.readline()
        line2 = map2.readline()

def _pair_gen(hitlist1, hitlist2):
    ''' does an in-order walk to find pairs '''

    # loop till each list is empty.
    while len(hitlist1) > 0 and len(hitlist2) > 0:

        # peek for a match.
        if hitlist1[-1][1] == hitlist2[-1][1]:

            # yield it.
            yield hitlist1[-1], hitlist2[-1]

            # change left.
            hitlist1.pop()

        else:

            # pop smaller.
            if hitlist1[-1][1] < hitlist2[-1][1]:
                hitlist1.pop()
            else:
                hitlist2.pop()
