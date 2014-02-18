'''
drivers to prepare the data for scaffolding
'''
import sys
import logging
import string
import subprocess
import os

from utils.misc import *

## private functions ##

## public functions ##
def prep_silp(work_dir, ant_dir, prep_sh, ctg_fasta, read1_sam, read2_sam, ins_size, std_dev, bundle_size, pair_mode, args):
    """ SILP scaffolder """

    # create command base.
    cmd_base = ['python', '/home/jrl03001/code/SILP2/silp.py']

    # preprocess.
    job = list()
    job.append('#!/bin/bash')
    job.append('# preprocess')
    job.append('%s nodes -w %s -c %s' % (' '.join(cmd_base), work_dir, ctg_fasta))
    job.append('%s edges -w %s -i %i -s %i -%s -s1 %s -s2 %s' %
        (' '.join(cmd_base), work_dir, ins_size, std_dev, pair_mode, read1_sam, read2_sam))
    job.append('%s bundles -w %s -b %i -p 90 -bup 1 -r %s -i %i -s %i' % (' '.join(cmd_base), work_dir, bundle_size, ant_dir, ins_size, std_dev))
    job.append('%s decompose -w %s -m 2500' % (' '.join(cmd_base), work_dir))
    job.append('')
    runit(job, prep_sh)

def prep_opera(work_dir, ant_dir, prep_sh, ctg_fasta, read1_sam, read2_sam, ins_size, std_dev, bundle_size, pair_mode, args):
    """ opera scaffolder """

    # opera specific directories.
    input_dir = '%s/input' % work_dir
    results_dir = '%s/results' % work_dir
    script_dir = '%s/script' % work_dir
    output_dir = '%s/output' % work_dir
    log_dir = '%s/log' % work_dir
    for x in [input_dir, results_dir, script_dir, output_dir, log_dir]:
        create_dir(x)

    # opera sepecific files.
    edge_file = '%s/opera_raw.txt' % input_dir
    config_file = '%s/opera.cfg' % input_dir
    script_file = '%s/opera.sh' % input_dir
    tout_file = '%s/scaffolds.scaf' % results_dir
    log_file = '%s/opera.log' % log_dir

    # open pair of sam files.
    fin1 = open(read1_sam, "rb")
    fin2 = open(read2_sam, "rb")
    fout = open(edge_file, 'wb')

    # create the edge file from paired sam.
    idx = 1
    for line1 in fin1:
        line2 = fin2.readline()

        # skip headers.
        if line1[0] == "@" or line2[0] == "@":
            continue

        # tokenize.
        tokens1 = line1.strip().split("\t")
        tokens2 = line2.strip().split("\t")

        # get data.
        rname1 = tokens1[2]
        rname2 = tokens2[2]

        qname1 = tokens1[0]
        qname2 = tokens2[0]

        pos1 = int(tokens1[3])
        pos2 = int(tokens2[3])

        szseq1 = len(tokens1[9])
        szseq2 = len(tokens2[9])

        seq1 = tokens1[9]
        seq2 = tokens2[9]

        # prepare orientation.
        if tokens1[1] == "0":
            orien1 = "+"
        else:
            orien1 = "-"
        if tokens2[1] == "0":
            orien2 = "+"
        else:
            orien2 = "-"

        # save names.
        nm1 = "%i.1" % idx
        nm2 = "%i.2" % idx
        idx += 1

        # create txt.
        line1 = "%s\t%s\t%s\t%i\t%s\t%s\t%s\n" % (nm1, orien1, rname1, pos1, tokens1[9], tokens1[10], "0")
        line2 = "%s\t%s\t%s\t%i\t%s\t%s\t%s\n" % (nm2, orien2, rname2, pos2, tokens2[9], tokens2[10], "0")

        # write to file.
        fout.write(line1)
        fout.write(line2)

    # close up shop.
    fin1.close()
    fin2.close()
    fout.close()

    # prepare config script
    txt = """#
# Essential Parameters
#

# Please always supply absolute path of each file,
# Because relative path may not work all the time.

# Output folder for final results
output_folder=${opera_results_dir}

# Contig file
contig_file=${opera_node_file}

#----------------------------------------------------------------------------------------------
# Advanced Parameters

#
# Scaffolding related parameters
#

# Scaffold name in result file
scaffold_name=scaffold

# PET cluster threshold (default=5) (Opera will discard all clusters

# Should Opera abort when running time for specific subgraph is longer
# than 30 minutes (true or false, default=true)
abort=false

#----------------------------------------------------------------------------------------------
#
# Contig file related parameters
#

# Format of contig file (fasta or statistic, default=fasta)
file_format=fasta

# Program name generating contig file (velvet or soap, default=velvet)
file_type=velvet

# Should the repeat contigs be filtered (yes or no, default=true)
filter_repeat=yes

# Repeat threshold (default=1.5): If the coverage of a contig is higher than
# (repeat threshold * average coverage), then it is considered as repeat
repeat_threshold=1.5

# Contig size threshold (default=500): Opera will not use the contigs whose length
# is shorter than this value
contig_size_threshold=10

#----------------------------------------------------------------------------------------------
#
# Library parameters.
#

[LIB]
calculate_ori=no
read_ori=in
map_type=bowtie
calculate_lib=no
lib_mean=${lib_mean}
lib_std=${lib_std}
cluster_threshold=${pet_size}
map_file=${opera_edge_file}
"""
    txt = string.Template(txt)
    txt = txt.substitute(opera_results_dir=results_dir, opera_node_file=ctg_fasta, opera_edge_file=edge_file, pet_size=bundle_size, lib_mean=ins_size, lib_std=std_dev)

    # write to file.
    fout = open(config_file, "wb")
    fout.write(txt)
    fout.close()


def prep_mip(work_dir, ant_dir, prep_sh, ctg_fasta, read1_sam, read2_sam, ins_size, std_dev, bundle_size, pair_mode, args):
    """ mip scaffolder """

    # simplify files.
    ant_dir = '/dev/null'
    pre_sh = prep_sh

    contig_file =  ctg_fasta
    test1 = read1_sam
    test2 = read2_sam
    
    # sanity.
    if os.path.isfile(test1) == False or os.path.isfile(test2) == False:
        logging.error("missing file:")
        logging.error(test1)
        logging.error(test2)
    
    # mip specific files.
    SAM_1_MOD = "%s/sam_1.sam" % work_dir
    SAM_2_MOD = "%s/sam_2.sam" % work_dir
    MERGED_FILE = "%s/merged" % work_dir
    MERGED1_FILE = "%s/merged.sorted1" % work_dir
    MERGED2_FILE = "%s/merged.sorted2" % work_dir
    FILTERED_FILE = "%s/filtered.txt" % work_dir
    COV_FILE = "%s/coverage.txt" % work_dir
    PARAM_FILE = "%s/parameters.txt" % work_dir
    MIP_FILE = "%s/scaffolds2.txt" % work_dir

    # modify SAM file to be happy with MIP.
    fin1 = open(test1, "rb")
    fin2 = open(test2, "rb")
    fout1 = open(SAM_1_MOD, "wb")
    fout2 = open(SAM_2_MOD, "wb")
    for line1 in fin1:
        line2 = fin2.readline()
        
        # tokenize.
        tokens1 = line1.strip().split("\t")
        tokens2 = line2.strip().split("\t")
        
        # change naming convention.
        tokens1[0] = tokens1[0].replace("/1","_R3")
        tokens2[0] = tokens2[0].replace("/2","_F3")
        
        # modify orientation.
        if pair_mode == 'rf':
            if tokens1[1] == "0":
                tokens1[1] = "16"
            else:
                tokens1[1] = "0"
        else:
            print pair_mode
            logging.error('un-written')
            sys.exit(1)
        
        
        # write it back out.
        fout1.write('\t'.join(tokens1) + '\n')
        fout2.write('\t'.join(tokens2) + '\n')

    fin1.close()    
    fin2.close()    
    fout1.close()   
    fout2.close()   

    PROG_1 = "/opt/mip/scripts/merge.sh"
    PROG_2 = "/opt/mip/scripts/filter-mappings.sh"
    PROG_3 = "/home/jrl03001/code/ScafValidate/mip/contig_coverage.py"
    #'''
    # merge data.
    with open(pre_sh, 'wb') as fout:
        
        # prepare.
        fout.write('#!/bin/bash\n')
        cmd = [PROG_1, SAM_2_MOD, SAM_1_MOD, MERGED_FILE]
        fout.write(' '.join(cmd) + '\n')
        cmd = [PROG_2, MERGED1_FILE, MERGED2_FILE, FILTERED_FILE]
        fout.write(' '.join(cmd) + '\n')
        cmd = ["python", PROG_3, contig_file, COV_FILE, SAM_1_MOD, SAM_2_MOD]
        fout.write(' '.join(cmd) + '\n')
        
    # run it.
    subprocess.call(["chmod", "u+x", pre_sh])
    if subprocess.call([pre_sh], cwd=work_dir) != 0:
        logging.error("couldn't prepare MIP")
        sys.exit(1)

    #'''
    # write the config.
    txt = '''# Upper bound for genome length (required)
genome_length=5000000000

#parameter specifications for the first stage
[STAGE]
# Maximum biconnected component size. (optional)
#maximum_biconnected_component=50
# Maximum allowed degree in scaffolding graph. (optional)
maximum_degree=50
# Maximum coverage for nonrepetitive contig. (optional)
#maximum_coverage=100
# The maximum overlap between contigs that is allowed without checking for
# sequence similarity. By default this is set based on the variablility in
# insert size lengths of each library. (optional)
#maximum_overlap=100
# The minimum support for an edge. (optional)
minimum_support=%i
# Should edges with negative estimated distance be checked for sequence
# similarity or removed automatically? (optional)
check_negative_edges=1

# library specification for the first stage
[LIBRARY]
# File in SAM format containing mappings for the mate pair reads
# to the contigs
mappings=%s
# Orientation of the mate pairs (in current version must be SOLID)
orientation=SOLID
# Insert length
insert_length=%i
# Minimum insert length
min_insert_length=%i
# Maximum insert length
max_insert_length=%i
''' % (bundle_size, FILTERED_FILE, ins_size, (ins_size - (3*std_dev)), (ins_size + (3*std_dev)))

    
    # write to file.
    fout = open(PARAM_FILE, "wb")
    fout.write(txt)
    fout.close()
