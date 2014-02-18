#!/usr/bin/python
'''
standalone adapter for scripts outside of the directory structure.
'''
### imports ###

# system
import time
import subprocess
import warnings
import argparse
import logging
import time
import sys
import os
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s', )

# local
from utils.align import create_idx, create_aln, pair_sam, pair_sam2
from utils.prep import prep_silp, prep_opera, prep_mip
from utils.run import run_silp, run_opera, run_mip
from utils.eval import simulation_evaluation, quast_evaluation

# hack to silence argparser.
warnings.filterwarnings('ignore', category=DeprecationWarning)

### classes ###

### functions ###

def align(args):
    ''' aligns reads against reference '''

    # validate parameters.
    assert os.path.isdir(args.base_dir), 'base_dir'
    assert os.path.isfile(args.ctg_fasta), 'ctg_fasta'
    assert os.path.isfile(args.read1_fastq), 'read1_fastq'
    assert os.path.isfile(args.read2_fastq), 'read2_fastq'
    assert os.path.isfile(args.size_file), 'size_file'

    # relavent files.
    base_dir = os.path.abspath(args.base_dir)

    size_file = os.path.abspath(args.size_file)
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    read1_fastq = os.path.abspath(args.read1_fastq)
    read2_fastq = os.path.abspath(args.read2_fastq)

    tmp1_sam = os.path.abspath('%s/tmp1.sam' % base_dir)
    tmp2_sam = os.path.abspath('%s/tmp2.sam' % base_dir)

    read1_sam = os.path.abspath('%s/read1.sam' % base_dir)
    read2_sam = os.path.abspath('%s/read2.sam' % base_dir)

    ant_dir = '%s/ant' % base_dir
    idx_dir = '%s/index' % base_dir
    idx_file = '%s/index' % idx_dir

    # build index if not present.
    if os.path.isdir(idx_dir) == False:
        subprocess.call(["mkdir", "-p", idx_dir])
        create_idx(ctg_fasta, idx_file)

    # remove annotation dir if present.
    if os.path.isdir(ant_dir) == True:
        subprocess.call(["rm", "-rf", ant_dir])
    subprocess.call(["mkdir", "-p", ant_dir])

    # perform alignment.
    create_aln(size_file, idx_file, read1_fastq, tmp1_sam, ant_dir, args.num_cpu)
    create_aln(size_file, idx_file, read2_fastq, tmp2_sam, ant_dir, args.num_cpu)

    # pair the alignment.
    pair_sam2(tmp1_sam, tmp2_sam, read1_sam, read2_sam, args.key_size)

def pair(args):
    ''' pairs two same files '''

    # validate parameters.
    assert os.path.isdir(args.base_dir), 'base_dir'
    assert os.path.isfile(args.tmp1_sam), 'tmp1_sam'
    assert os.path.isfile(args.tmp2_sam), 'tmp2_sam'

    # relavent files.
    base_dir = os.path.abspath(args.base_dir)
    tmp1_sam = os.path.abspath(args.tmp1_sam)
    tmp2_sam = os.path.abspath(args.tmp2_sam)
    read1_sam = os.path.abspath('%s/read1.sam' % base_dir)
    read2_sam = os.path.abspath('%s/read2.sam' % base_dir)

    # pair the alignment.
    pair_sam2(tmp1_sam, tmp2_sam, read1_sam, read2_sam, args.key_size)

def meta_combine(args):
    """ combines two seperate alignments"""

    # validate parameters.
    assert os.path.isdir(args.base_dir), 'base_dir'
    assert os.path.isdir(args.work1_dir), 'w1_dir'
    assert os.path.isdir(args.work2_dir), 'w2_dir'

    # switch based on subset.
    if args.subsample == None:

        # combine reads.
        with open("%s/read1.sam" % args.base_dir, "wb") as fout:
            subprocess.call(["cat","%s/read1.sam" % args.work1_dir, "%s/read1.sam" % args.work2_dir], stdout=fout)

        with open("%s/read2.sam" % args.base_dir, "wb") as fout:
            subprocess.call(["cat","%s/read2.sam" % args.work1_dir, "%s/read2.sam" % args.work2_dir], stdout=fout)

    else:

        # count number of lines.
        cnt = 0
        with open("%s/read1.sam" % args.work1_dir, "r") as fin:
            for line in fin:
                cnt += 1
                               
        # only output a percentage of them.
        outp = float(cnt) * float(args.subsample)
        
        # first one.
        with open("%s/read1.sam" % args.base_dir, "wb") as fout:
            cnt = 0
            with open("%s/read1.sam" % args.work1_dir, "r") as fin:
                for line in fin:
                    if cnt > outp: break
                    fout.write(line)
                    cnt += 1
            with open("%s/read1.sam" % args.work2_dir, "r") as fin:
                for line in fin:
                    fout.write(line)
                    
        # second one.
        with open("%s/read2.sam" % args.base_dir, "wb") as fout:
            cnt = 0
            with open("%s/read2.sam" % args.work1_dir, "r") as fin:
                for line in fin:
                    if cnt > outp: break
                    fout.write(line)
                    cnt += 1
            with open("%s/read2.sam" % args.work2_dir, "r") as fin:
                for line in fin:
                    fout.write(line)

    # copy one annotation.
    subprocess.call(["cp", "%s/ant" % args.work1_dir, "%s/ant" % args.base_dir, "-R"])
    subprocess.call(["chmod", "u+w", "%s/ant" % args.base_dir, "-R"])

    # combine the rest.
    tmp1_dir = "%s/ant" % args.base_dir
    tmp2_dir = "%s/ant" % args.work2_dir
    for x in os.listdir(tmp2_dir):

        # simplify.
        tmp1_x = "%s/%s" % (tmp1_dir,x)
        tmp2_x = "%s/%s" % (tmp2_dir,x)
        
        # copy it if no problem.
        if os.path.isfile(tmp1_x) == False:
            subprocess.call(["cp", tmp2_x, tmp1_x])

        else:
            # combine it.
            tmp1_n = np.load(tmp1_x)
            tmp2_n = np.load(tmp2_x)
            np.save(tmp1_x, tmp1_n+tmp2_n)

def prepare(args):
    """ prepares alignment given scaffolding method"""

    # extract information.
    work_dir = os.path.abspath(args.base_dir)
    align_dir = os.path.abspath(args.align_dir)
    ant_dir = '%s/ant' % align_dir

    read1_sam = '%s/read1.sam' % align_dir
    read2_sam = '%s/read2.sam' % align_dir
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    prep_sh = '%s/prep.sh' % work_dir

    # translate pairmode.
    if args.pair_mode == 0:
        pair_mode = 'ff'
    elif args.pair_mode == 1:
        pair_mode = 'fr'
    elif args.pair_mode == 2:
        pair_mode = 'rf'

    # run the right prepare.
    cmd_args = (work_dir, ant_dir, prep_sh, ctg_fasta, read1_sam, read2_sam, args.ins_size, args.std_dev, args.bundle_size, pair_mode, args)
    if args.meth_mode == 0:
        prep_silp(*cmd_args)
    elif args.meth_mode == 1:
        prep_opera(*cmd_args)
    elif args.meth_mode == 2:
        prep_mip(*cmd_args)


def run(args):
    """ extracts information for scaffolding"""

    # extract information.
    work_dir = os.path.abspath(args.base_dir)
    align_dir = os.path.abspath(args.align_dir)
    ant_dir = '%s/ant' % align_dir

    read1_sam = '%s/read1.sam' % align_dir
    read2_sam = '%s/read2.sam' % align_dir
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    ctg_agp = '%s/scf.agp' % work_dir
    scf_fasta = '%s/scf.fasta' % work_dir
    run_sh = '%s/run.sh' % work_dir

    # translate pairmode.
    if args.pair_mode == 0:
        pair_mode = 'ff'
    elif args.pair_mode == 1:
        pair_mode = 'fr'
    elif args.pair_mode == 2:
        pair_mode = 'rf'

    # run the right prepare.
    cmd_args = (work_dir, ant_dir, run_sh, ctg_fasta, scf_fasta, ctg_agp, read1_sam, read2_sam, args.ins_size, args.std_dev, args.bundle_size, pair_mode, args)
    if args.meth_mode == 0:
        run_silp(*cmd_args)
    elif args.meth_mode == 1:
        run_opera(*cmd_args)
    elif args.meth_mode == 2:
        run_mip(*cmd_args)

def sim_eval(args):
    """evaluates scaffold"""

    # simplify data.
    ref_fasta = os.path.abspath(args.ref_fasta)
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    scf_fasta = os.path.abspath(args.scf_fasta)
    ref_agp = os.path.abspath(args.ref_agp)
    scf_agp = os.path.abspath(args.scf_agp)

    # call evaluation code.
    A, B, C, D, ppv, mcc, gap_dev, runtime, ref_N50, ctg_N50, scf_N50, tp_N50 = simulation_evaluation(ref_fasta, ctg_fasta, scf_fasta, ref_agp, scf_agp)

    # format and return.
    txt = list()
    txt.append('%d %d %d %d' % (A, B, C, D))
    txt.append('%.2f %.2f' % (ppv, mcc))
    txt.append('%.2f %.2f' % (gap_dev, runtime))
    txt.append('%d %d %d %d' % (ref_N50, ctg_N50, scf_N50, tp_N50))
    print ' '.join(txt)

def real_eval(args):
    """evaluates scaffold using alignme"""

    # simplify data.
    ref_fasta = os.path.abspath(args.ref_fasta)
    ctg_fasta = os.path.abspath(args.ctg_fasta)
    scf_fasta = os.path.abspath(args.scf_fasta)
    #ref_agp = os.path.abspath(args.ref_agp)
    #scf_agp = os.path.abspath(args.scf_agp)

    # call quast
    res = quast_evaluation(ref_fasta, ctg_fasta, scf_fasta, args.wd, args.threads)


    # quast specific.
    try: N50 = int(float(res[1,np.where(res[0,:] == "N50")[0][0]]))
    except: N50 = 0
    try: NG50 = int(float(res[1,np.where(res[0,:] == "NG50")[0][0]]))
    except: NG50 = 0
    try: GFRAC = float(res[1,np.where(res[0,:] == "Genome fraction (%)")[0][0]])
    except: GFRAC = 0.0
    try: NA50 = int(float(res[1,np.where(res[0,:] == "NA50")[0][0]]))
    except: NA50 = 0
    try: NGA50 = int(float(res[1,np.where(res[0,:] == "NGA50")[0][0]]))
    except: NGA50 = 0

    # format and return.
    txt = list()
    txt.append("%d %d %.2f %d %d" % (N50, NG50, GFRAC, NA50, NGA50))

    print ' '.join(txt)

### script ###

if __name__ == '__main__':

    # mode parser.
    main_p = argparse.ArgumentParser()
    subp = main_p.add_subparsers(help='sub-command help')

    ### data preparation ###
    # import reference into working directory.
    subp_p = subp.add_parser('align', help='aligns reads against contigs.')
    subp_p.add_argument('-w', dest='base_dir', required=True, help='working directory')
    subp_p.add_argument('-p', dest='num_cpu', type=int, required=True, help='number of threads for bowtie2')
    subp_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    subp_p.add_argument('-q1', dest='read1_fastq', required=True, help='read first file')
    subp_p.add_argument('-q2', dest='read2_fastq', required=True, help='read second file')
    subp_p.add_argument('-s', dest='size_file', required=True, help='size file')
    subp_p.add_argument('-k', dest='key_size', type=int, required=True, help='size of PE key at end of each read')
    subp_p.set_defaults(func=align)

    subp_p = subp.add_parser('pair', help='pairs 2 SAM files')
    subp_p.add_argument('-w', dest='base_dir', required=True, help='working directory')
    subp_p.add_argument('-s1', dest='tmp1_sam', required=True, help='read first file')
    subp_p.add_argument('-s2', dest='tmp2_sam', required=True, help='read second file')
    subp_p.add_argument('-k', dest='key_size', type=int, required=True, help='size of PE key at end of each read')
    subp_p.set_defaults(func=pair)

    # combines two seperate alignments.
    subp_p = subp.add_parser('meta_combine', help='combines two seperate alignments. Used for meta scaffold testing.')
    subp_p.add_argument('-w', dest='base_dir', required=True, help='working directory')
    subp_p.add_argument('-w1', dest='work1_dir', required=True, help='working directory')
    subp_p.add_argument('-w2', dest='work2_dir', required=True, help='working directory')
    subp_p.add_argument('-s', dest='subsample', type=float, required=False, help='subsample first directory')
    subp_p.set_defaults(func=meta_combine)

    # prepare the reference for scaffolding.
    subp_p = subp.add_parser('prep', help='runs pre-processing stage.')
    subp_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    subp_p.add_argument('-w', dest='base_dir', required=True, help='scaffolding directory')
    subp_p.add_argument('-a', dest='align_dir', required=True, help='alignment directory')
    subp_p.add_argument('-i', dest='ins_size', type=int, required=True, help='insert size')
    subp_p.add_argument('-s', dest='std_dev', type=int, required=True, help='std_devition')
    subp_p.add_argument('-b', dest='bundle_size', type=int, required=True, help='bundle size')
    me_g = subp_p.add_mutually_exclusive_group(required=True)
    me_g.add_argument('-ff', dest='pair_mode', action='store_const', const=0, help='SOLiD style -> ->')
    me_g.add_argument('-fr', dest='pair_mode', action='store_const', const=1, help='innie style -> <-')
    me_g.add_argument('-rf', dest='pair_mode', action='store_const', const=2, help='outtie style <- ->')
    me_g = subp_p.add_mutually_exclusive_group(required=True)
    me_g.add_argument('-silp', dest='meth_mode', action='store_const', const=0, help='SILP')
    me_g.add_argument('-opera', dest='meth_mode', action='store_const', const=1, help='OPERA')
    me_g.add_argument('-mip', dest='meth_mode', action='store_const', const=2, help='MIP')
    subp_p.set_defaults(func=prepare)

    # run the scaffolding.
    subp_p = subp.add_parser('run', help='runs the actual scaffolding')
    subp_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    subp_p.add_argument('-w', dest='base_dir', required=True, help='scaffolding directory')
    subp_p.add_argument('-a', dest='align_dir', required=True, help='alignment directory')
    subp_p.add_argument('-i', dest='ins_size', type=int, required=True, help='insert size')
    subp_p.add_argument('-s', dest='std_dev', type=int, required=True, help='std_devition')
    subp_p.add_argument('-b', dest='bundle_size', type=int, required=True, help='bundle size')
    me_g = subp_p.add_mutually_exclusive_group(required=True)
    me_g.add_argument('-ff', dest='pair_mode', action='store_const', const=0, help='SOLiD style -> ->')
    me_g.add_argument('-fr', dest='pair_mode', action='store_const', const=1, help='innie style -> <-')
    me_g.add_argument('-rf', dest='pair_mode', action='store_const', const=2, help='outtie style <- ->')
    me_g = subp_p.add_mutually_exclusive_group(required=True)
    me_g.add_argument('-silp', dest='meth_mode', action='store_const', const=0, help='SILP')
    me_g.add_argument('-opera', dest='meth_mode', action='store_const', const=1, help='OPERA')
    me_g.add_argument('-mip', dest='meth_mode', action='store_const', const=2, help='MIP')
    # SILP2 arguments.
    subp_p.add_argument('-z', dest='weight_mode', type=int, default=0, help='SILP2: weight mode')
    subp_p.set_defaults(func=run)

    # evaluate the scaffolding.
    subp_p = subp.add_parser('sim_eval', help='evaluates scaffolding if there is a true AGP')
    subp_p.add_argument('-r', dest='ref_fasta', required=True, help='reference fasta')
    subp_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    subp_p.add_argument('-s', dest='scf_fasta', required=True, help='scaffold fasta')
    subp_p.add_argument('-a', dest='ref_agp', required=True, help='reference agp')
    subp_p.add_argument('-y', dest='scf_agp', required=True, help='predicted agp')
    subp_p.set_defaults(func=sim_eval)

    # evaluate the scaffolding.
    subp_p = subp.add_parser('real_eval', help='alignment based scaffold evaluation.')
    subp_p.add_argument('-r', dest='ref_fasta', required=True, help='reference fasta')
    subp_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    subp_p.add_argument('-s', dest='scf_fasta', required=True, help='scaffold fasta')
    subp_p.add_argument('-y', dest='scf_agp', required=True, help='predicted agp')
    subp_p.add_argument('-wd', dest='wd', required=True, help='alignment working directory')
    subp_p.add_argument('-p', dest='threads', type=int, required=True, help='number of threads to align')
    subp_p.set_defaults(func=real_eval)

    args = main_p.parse_args()
    args.func(args)
