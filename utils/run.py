'''
drivers to prepare the data for scaffolding
'''
import subprocess
import os
import time
import numpy as np

from utils.misc import *

## private functions ##

## public functions ##
def run_silp(work_dir, ant_dir, run_sh, ctg_fasta, scf_fasta, ctg_agp, read1_sam, read2_sam, ins_size, std_dev, bundle_size, pair_mode, args):
    """ SILP scaffolder """

    # create command base.
    cmd_base = ['python', '/home/jrl03001/code/SILP2/silp.py']

    # scaffolding/
    job = list()
    job.append('#!/bin/bash')
    job.append('# scaffolding')
    job.append('# start time')
    job.append('start=$(date +%s)')
    job.append('%s orient -w %s -z %s ' % (' '.join(cmd_base), work_dir, args.weight_mode))
    job.append('%s order -w %s ' % (' '.join(cmd_base), work_dir))
    job.append('%s gap -w %s ' % (' '.join(cmd_base), work_dir))
    job.append('%s write -w %s -a %s' % (' '.join(cmd_base), work_dir, ctg_agp))
    job.append('%s fasta -w %s -a %s -c %s -f %s' % (' '.join(cmd_base), work_dir, ctg_agp, ctg_fasta, scf_fasta))
    job.append('# stop time')
    job.append('stop=$(date +%s)')
    job.append('echo RUNTIME: $(expr $stop - $start) >> %s' % ctg_agp)
    job.append('')
    runit(job, run_sh)


def run_opera(work_dir, ant_dir, run_sh, ctg_fasta, scf_fasta, ctg_agp, read1_sam, read2_sam, ins_size, std_dev, bundle_size, pair_mode, args):
    """ OPERA scaffolder """

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
    tscf_file = '%s/scaffoldSeq.fasta' % results_dir
    log_file = '%s/opera.log' % log_dir

    # write execution script.
    job = list()
    job.append('#!/bin/bash')
    job.append('# scaffolding')
    job.append('# start time')
    job.append('start=$(date +%s)')
    job.append('# run the job')
    job.append('/opt/opera/bin/opera %s &> %s' % (config_file, log_file))
    job.append('# stop time')
    job.append('stop=$(date +%s)')
    job.append('echo RUNTIME: $(expr $stop - $start) >> %s' % ctg_agp)
    job.append('')
    runit(job, run_sh)

    # load results.
    fin = open(tout_file, "rb")
    lines = fin.readlines()
    fin.close()

    # load the runtime string.
    with open(ctg_agp, 'rb') as fin:
        time_info = fin.readline().strip()

    # create an AGP.
    fout = open(ctg_agp, "wb")
    for lidx in range(len(lines)):

        # setup line.
        line = lines[lidx]

        # check header.
        if line[0] == ">":
            header = line.strip().replace(">","").split()[0]
            #header = line.strip().replace(">","")
            part = 1
            scaf_start = 1
            scaf_stop = 0
            continue

        # tokenize.
        tmp = line.strip().split()

        # get ctgid.
        ctgid = tmp[0]

        if tmp[1] == "BE":
            orien = "+"
        else:
            orien = "-"

        ctglen = int(tmp[2])
        gaplen = int(tmp[3])

        # make indexs.
        scaf_stop = scaf_start + ctglen

        # write out AGP.
        fout.write("%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t%s\n" % \
            (header, scaf_start, scaf_stop, part, "W", ctgid, 1, ctglen, orien))
        #print "write"
        # increment pointers.
        scaf_start = scaf_stop + 1
        part += 1

        # add gap if size is not 0
        if gaplen != 0:

            # make indexs.
            scaf_stop = scaf_start + gaplen

            # add gap.
            fout.write("%s\t%i\t%i\t%i\t%s\t%i\t%s\t%s\n" % \
                (header, scaf_start, scaf_stop, part, "N", gaplen, "fragment", "no"))
            #print "write"
            # increment pointers.
            scaf_start = scaf_stop + 1
            part += 1

        # otherwise.
        else:
            # check if we can peak.
            if lidx+1 >= len(lines):
                continue


            # check if new scaffold.
            if lines[lidx+1][0] != ">":
                # make fake gap.
                scaf_stop = scaf_start + 10

                # add gap.
                fout.write("%s\t%i\t%i\t%i\t%s\t%i\t%s\t%s\n" % \
                    (header, scaf_start, scaf_stop, part, "N", gaplen, "fragment", "no"))
                #print "write"
                # increment pointers.
                scaf_start = scaf_stop + 1
                part += 1
                continue

    # write the timer.
    fout.write(time_info + '\n')

    # close AGP file.
    fout.close()

    # move the output files to standard locations.
    subprocess.call(['cp', tscf_file, scf_fasta])


def run_mip(work_dir, ant_dir, run_sh, ctg_fasta, scf_fasta, ctg_agp, read1_sam, read2_sam, ins_size, std_dev, bundle_size, pair_mode, args):
    
    # simplify files.
    ant_dir = '/dev/null'

    contig_file =  ctg_fasta
    test1 = read1_sam
    test2 = read2_sam
    
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
    OUT_FILE = "%s/scaffolds.fasta" % work_dir

    
    # execute mip.
    tstart = time.time()
    cmd = list()
    cmd.append("/opt/mip/scripts/mip-scaffolder.pl")
    cmd.append(PARAM_FILE)
    cmd.append(contig_file)
    cmd.append(COV_FILE)
    cmd.append(work_dir)
    
    with open(run_sh, "wb") as fout:
        fout.write('#!/bin/bash\n')
        fout.write(' '.join(cmd) + '\n')
        
    subprocess.call(["chmod", "u+x", run_sh])
    if subprocess.call([run_sh], cwd=work_dir) != 0:
        logging.warning("couldn't run mip")
        sys.exit(1)
    tstop = time.time()
    
    # compute time.
    trun = tstop - tstart
    
    # load MIP results into dictionary.
    fin = open(MIP_FILE, "rb")
    lines = fin.readlines()
    fin.close()
    mip = {}
    for line in lines:
        
        # tokenize.
        tmp = line.split()
        
        # set current scaffold.
        if tmp[0].count("scaffold") > 0:
            # set current.
            curs = tmp[0]
            
            # set default.
            mip[curs] = []
            continue
            
        # append.
        mip[curs].append(tmp)
        
    # count the number of rows for agp file.
    size = 0
    for scafid in mip:
        size += len(mip[scafid]) + len(mip[scafid]) - 1

    # allocate the array.
    agps = np.zeros(size, dtype=agp_dt)

    # begin copying data.
    idx = 0
    for scafid in mip:
        scaf_start = 1
        scaf_idx = 1
        for i in range(len(mip[scafid])):
            
            # tokenize.
            entry = mip[scafid][i]
            ctg_name = entry[1]
            if entry[2] == "F":
                orien = 0
            else:
                orien = 1
            start = int(entry[3])
            stop = int(entry[4])
            
            # add to agp.
            agps[idx]['scaf_name'] = scafid
            agps[idx]['scaf_start'] = scaf_start
            agps[idx]['scaf_stop'] = scaf_start + abs(stop - start)
            agps[idx]['scaf_idx'] = scaf_idx
            agps[idx]['comp_type'] = "W"
            agps[idx]['comp_name'] = ctg_name
            agps[idx]['comp_start'] = 1
            agps[idx]['comp_stop'] = abs(stop - start)
            agps[idx]['comp_orien'] = orien
            agps[idx]['comp_linkage'] = 0

            # move up counts.
            scaf_start = agps[idx]['scaf_stop'] + 1
            scaf_idx += 1
            idx += 1
            
            # do gap.
            if entry[1] != mip[scafid][-1][1]:
                    
                # add gap.
                agps[idx]['scaf_name'] = scafid
                agps[idx]['scaf_start'] = scaf_start
                agps[idx]['scaf_stop'] = int(mip[scafid][i+1][3]) - 1
                agps[idx]['scaf_idx'] = scaf_idx
                agps[idx]['comp_type'] = "N"
                agps[idx]['comp_name'] = "fragment"
                agps[idx]['comp_start'] = 1
                agps[idx]['comp_stop'] = abs(agps[idx]['scaf_stop'] - agps[idx]['scaf_start'])
                agps[idx]['comp_orien'] = 0
                agps[idx]['comp_linkage'] = 0   
                
                # move up counts.
                scaf_start = int(mip[scafid][i+1][3])
                scaf_idx += 1
                idx += 1

    # copy results.
    subprocess.call(['cp', OUT_FILE, scf_fasta])

    # save the agp.
    save_agps(ctg_agp, agps)

    # add time.
    with open(ctg_agp, 'rb') as fin:
        lines = fin.read()
        lines += 'RUNTIME: %.2f' % trun
    with open(ctg_agp, 'wb') as fout:
        fout.write(lines)
        
    # move the results.
