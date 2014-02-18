import os
import sys
import math
import numpy as np
import logging
import subprocess
import networkx as nx

### definitions ###
agp_dt = np.dtype([\
    ('scaf_name', 'S255'),\
    ('scaf_start', np.long),\
    ('scaf_stop', np.long),\
    ('scaf_idx', np.long),\
    ('comp_type', 'S50'),\
    ('comp_name', 'S255'),\
    ('comp_start', np.long),\
    ('comp_stop', np.long),\
    ('comp_orien', np.long),\
    ('comp_linkage', np.long),\
])

### private functions ###

def _load_agp(fpath):
    ''' read agp file into array.'''

    # read in agp.
    fin = open(fpath, "rb")
    lines = fin.readlines()
    fin.close()

    # count number of lines minus comments.
    cnt = 0
    for line in lines:
        if line[0] != "#" and len(line) != 0 and line.strip().split()[0] != "RUNTIME:":
            cnt += 1

    # instantiate array.
    agp_edges = np.zeros(cnt, dtype=agp_dt)

    # parse agp.
    idx = 0
    for line in lines:
        # tokenize.
        if line[0] == "#": continue
        tmp = line.strip().split()
        if len(tmp) == 0: continue
        if tmp[0] == "RUNTIME:": continue

        # get general tokenize.
        agp_edges[idx]['scaf_name'] = tmp[0]
        agp_edges[idx]['scaf_start'] = int(float(tmp[1]))
        agp_edges[idx]['scaf_stop'] = int(float(tmp[2]))
        agp_edges[idx]['scaf_idx'] = int(float(tmp[3]))
        agp_edges[idx]['comp_type'] = tmp[4]

        # contig.
        if tmp[4] == "W":
            # get parts.
            agp_edges[idx]['comp_name'] = tmp[5]
            agp_edges[idx]['comp_start'] = int(tmp[6])
            agp_edges[idx]['comp_stop'] = int(tmp[7])
            if tmp[8] == "+":
                agp_edges[idx]['comp_orien'] = 0
            else:
                agp_edges[idx]['comp_orien'] = 1

        else:

            # save entry.
            agp_edges[idx]['comp_name'] = tmp[6]
            agp_edges[idx]['comp_start'] = 1
            agp_edges[idx]['comp_stop'] = int(tmp[5])
            if tmp[7] != "yes":
                agp_edges[idx]['comp_linkage'] = 0
            else:
                agp_edges[idx]['comp_linkage'] = 1


        # update index.
        idx += 1

    # shirnk array.
    agp_edges.resize(idx)

    return agp_edges

def _agp_graph(fpath, RG=None):
    ''' returns agp graph '''

    # load agp array.
    agp_edges = _load_agp(fpath)

    # make digraph.
    G = nx.DiGraph()

    # add nodes.
    for i in range(agp_edges.size):

        # skip contigs themselves.
        if agp_edges[i]['comp_type'] != 'W': continue

        # add node info.
        name = agp_edges[i]['comp_name']
        width = agp_edges[i]['comp_stop']
        orien = agp_edges[i]['comp_orien']
        start = agp_edges[i]['scaf_start']
        stop = agp_edges[i]['scaf_stop']
        scaf = agp_edges[i]['scaf_name']
        G.add_node(agp_edges[i]['comp_name'], {'width':width, 'orien':orien, 'start':start, 'stop':stop,'scaf':scaf})

    # add edges.
    for i in range(agp_edges.size):

        # skip contigs themselves.
        if agp_edges[i]['comp_type'] != 'N': continue
        if i == agp_edges.shape[0] - 1: continue

        # add sorted edges.
        ctg1 = agp_edges[i-1]['comp_name']
        ctg2 = agp_edges[i+1]['comp_name']
        gap = agp_edges[i]['comp_stop']

        G.add_edge(ctg1, ctg2, {'gap':gap})

    # if RG is supplied ensure all contigs exist.
    if RG != None:
        for n in RG.nodes():
            if G.has_node(n) == False:
                G.add_node(n, RG.node[n])

    # done.
    return G

def _uconn_N50(RG, TG):
    """ calcultes UCONN N50s"""

    # flip TG to be most consistent with RG.
    TG = _dag_flip(RG, TG)

    # calculate teh TPN50
    min_size = 1
    WG_REFN50 = _calculate_n50(RG, gap=True, min_size=min_size)
    NG_REFN50 = _calculate_n50(RG, gap=False, min_size=min_size)

    # calculate reported n50
    WG_N50 = _calculate_n50(TG, gap=True, min_size=min_size)
    NG_N50 = _calculate_n50(TG, gap=False, min_size=min_size)

    # remove edges not in TG.
    TG = _remove_difference(RG, TG)

    # calculate teh TPN50
    WG_TPN50 = _calculate_n50(TG, gap=True, min_size=min_size)
    NG_TPN50 = _calculate_n50(TG, gap=False, min_size=min_size)

    # yield results.
    return WG_N50, WG_TPN50, TG

def _compute_formula(true_adj, test_adj, n):
    ''' computes statistics when TG is already oriented. This
    method computes A and B exactly, but uses formulas
    to find C and D.'''

    # set vars.
    m = len(true_adj)
    N = n * (n-1)

    # compute the A statistic set.
    A = true_adj.intersection(test_adj)

    # compute the B statistic set.
    B = test_adj.difference(true_adj)

    # compute the C statistic set.
    C = m - len(A)

    # compute the D statistic set.
    D = N - m - len(B)

    # return categories.
    return len(A), len(B), C, D, A


def _gap_deviation(RG, TG):
    ''' calculate MAPE '''

    # find true edges.
    radj = set(RG.edges())
    tadj = set(TG.edges())
    edges = tadj.intersection(radj)

    # calculate MAPE.
    s = 0.0
    for e0, e1 in edges:
        At = RG[e0][e1]['gap']
        Ft = TG[e0][e1]['gap']
        if At != 0.0:
            s += abs(float(At - Ft) / float(At))

    # return it.
    try:
        return (100.0 / float(len(edges))) * s
    except:
        return -1.0

def _get_runtime(agp_file):

    # read in agp.
    fin = open(agp_file, "rb")
    lines = fin.readlines()
    fin.close()

    if len(lines) == 0:
        return -1.0

    if lines[-1].count("RUNTIME") == 0:
        return -1.0

    # parse last line.
    tmp = lines[-1].strip().split()

    return float(tmp[1])


def _sanity_check(RG, TG, A, B, C, D):
    ''' make sure these parameter make sense '''

    # set vars.
    n = RG.number_of_nodes()
    m = RG.number_of_edges()
    N = n * (n-1)

    # run check 1: m = A + C.
    check1 = m == A + C

    # run check 2: N - m = B + D.
    check2 = N - m == B + D

    # run report.
    fail = False

    if check1 == False:
        logging.error("failed check 1: m = A + C, %i = %i + %i: %i" % (m, A, C, m == A + C))
        fail = True

    if check2 == False:
        logging.error("failed check 2: N-m = B + D, %i - %i = %i + %i: %i" % (N, m, B, D, N-m==B+D))
        fail = True

    if fail == True:
        sys.exit()


def _dag_flip(RG, TG):
    ''' flip to make most consistant '''

    # check it.
    _graph_check(RG)
    _graph_check(TG)

    # make set of reference.
    rset = set(RG.edges())

    # loop over each component.
    NG = nx.DiGraph()
    for comp in nx.weakly_connected_components(TG):

        # turn to subgraph.
        subg = TG.subgraph(comp)

        # make sets.
        tset1 = set(subg.edges())
        tset2 = set([(e1, e0) for e0, e1 in tset1])

        # check version.
        s1 = len(rset.intersection(tset1))
        s2 = len(rset.intersection(tset2))

        # add to nodes to new graph.
        for n in subg.nodes():
            NG.add_node(n, subg.node[n])

        # add the best matching orientation.
        for e0, e1 in subg.edges():
            if s1 >= s2:
                NG.add_edge(e0, e1, subg[e0][e1])
            else:
                NG.add_edge(e1, e0, subg[e0][e1])

    # return it.
    return NG

def _graph_check(G):
    ''' makes sure its linear and a DAG '''

    # path check.
    for n in G.nodes():
        if G.number_of_edges(n) > 2:
            logging.error("bad scaffold graph 1 ")
            sys.exit(1)

    # DAG check.
    if nx.is_directed_acyclic_graph(G) == False:
        logging.error("bad scaffold graph 2")
        sys.exit(1)

def _calculate_n50(G, min_size=False, gap=False):
    ''' calculates scaffold N50 given scaffold graph'''

    # compute the scaffold sizes.
    sizes = _calc_sizes(G, min_size,gap)

    # calculate n50.
    sizes.sort(reverse = True)
    s = sum(sizes)
    limit = s * 0.5
    for l in sizes:
        s -= l
        if s <= limit:
            return int(l)
            
def _fasta_n50(fasta_file, min_size=False, gap=False):
    ''' calculates scaffold N50 given fasta file'''

    # load fasta.
    seqs = _load_fasta(fasta_file)

    # compute the scaffold sizes.
    sizes = [len(seqs[x]) for x in seqs]

    # calculate n50.
    sizes.sort(reverse = True)
    s = sum(sizes)
    limit = s * 0.5
    for l in sizes:
        s -= l
        if s <= limit:
            return int(l)


def _calc_sizes(G, min_size, gap, scaf_only=False):
    sizes = list()
    for comp in nx.weakly_connected_components(G):

        # skip non scaffolds.
        if scaf_only == True:
            if len(comp) < 2:
                continue

        # add contig size.
        size = 0
        for n in comp:
            size += G.node[n]['width']

        # add gap size.
        if gap != False:
            for p,q in G.edges(comp):
                size += G[p][q]['gap']

        # skip this.
        if min_size != False and size < min_size:
            continue

        # save the size.
        sizes.append(size)
    return sizes


def _calc_counts(G, min_size, gap, scaf_only=False):
    sizes = list()
    for comp in nx.weakly_connected_components(G):

        # skip non scaffolds.
        if scaf_only == True:
            if len(comp) < 2:
                continue

        # add contig size.
        size = 0
        for n in comp:
            size += G.node[n]['width']

        # add gap size.
        if gap != False:
            for p,q in G.edges(comp):
                size += G[p][q]['gap']

        # skip this.
        if min_size != False and size < min_size:
            continue

        # save the size.
        sizes.append(len(comp))

    return sizes

def _remove_difference(RG, TG):
    ''' removes edges not in RG'''

    # make edge sets.
    rset = set(RG.edges())
    tset = set(TG.edges())

    # identify edges in test that shouldn't be there
    to_remove = tset.difference(rset)

    # remove them.
    TG.remove_edges_from(to_remove)

    # return modified graph.
    return TG

def _load_fasta(file_path):
    ''' loads fasta file into dictionary'''

    # read file into memory.
    fin = open(file_path)
    lines = fin.readlines()
    fin.close()

    # build dictionary.
    data = dict()
    seq = ""
    for line in lines:

        # Skip blanks.
        if len(line) < 2: continue
        if line[0] == "#": continue

        # remove blanks.
        line = line.strip()

        # Check for ids.
        if line.count(">") > 0:

            # Check if ending seq.
            if len(seq) > 0:

                # save.
                data[head] = seq.upper()

            # reset head.
            head = line.replace(">","")
            seq = ""

            # skip to next line.
            continue

        # Filter chars.
        seq += line

    # save the last one.
    data[head] = seq.upper()

    # return dictionary.
    return data

### public functions ###

def simulation_evaluation(ref_fasta, ctg_fasta, scf_fasta, ref_agp, scf_agp):
    """ evaluates the simulated dataset """

    # load the reference and test graph from AGP files.
    RG = _agp_graph(ref_agp)
    TG = _agp_graph(scf_agp, RG=RG)

    # calculate referece and ctg n50/
    ref_N50 = _fasta_n50(ref_fasta, gap=True)
    ctg_N50 = _fasta_n50(ctg_fasta, gap=True)
    scf_N50 = _fasta_n50(scf_fasta, gap=True)

    # ensure we have all contigs accounted for.
    assert set(TG.nodes()) == set(RG.nodes()),\
        'not all contigs accounted for'

    # compute the exact N50
    scf_N50_2, tp_N50, LG = _uconn_N50(RG, TG)

    # make edge sets.
    radj = set(RG.edges())
    tadj = set(TG.edges())

    # calculate 4 parameters.
    n = RG.number_of_nodes()
    A, B, C, D, A_SET = _compute_formula(radj, tadj, n)

    # compute derivativ stats.
    sensitivity = (100*(float(A) / float(A+C)))
    ppv = (100*(float(A) / float(A+B)))
    if float(math.sqrt((A+B)*(A+C)*(D+B)*(D+C))) != 0.0:
        mcc = (100*(float((A*D)-(B*C)) / float(math.sqrt((A+B)*(A+C)*(D+B)*(D+C)))))
    else:
        mcc = 0.0

    # get gap deviation.
    gap_dev = _gap_deviation(RG, TG)

    # get the runtime.
    runtime = _get_runtime(scf_agp)

    # print results.
    # A B C D gapMape runtime sensitivity ppv mcc REFN50 TPN50 REPN50
    return A, B, C, D, ppv, mcc, gap_dev, runtime, ref_N50, ctg_N50, scf_N50, tp_N50



def quast_evaluation(ref, ctgs, scf, sdir, threads):

    # note that we are running it.
    logging.info("running quast %s %s %s" % (ref, ctgs, scf))

    # run the comparison script.
    x = ref.split("/")
    #z = '/'.join(x[0:-1] + ["genes.gff"])

    cmd = [\
        'python',\
        '/opt/quast/quast.py',\
        "-e",
        "-f",
     #   "-G", z,\
        "-o", "./",\
        "-R", ref,\
        "-T", str(threads),\
        scf\
    ]

    '''
    # remove dir.
    if os.path.isdir(sdir) == True:
        subprocess.call(["rm", "-rf", sdir])

    # make dir.
    if os.path.isdir(sdir) == False:
        subprocess.call(['mkdir', '-p', sdir])

    # make log file.
    lfile = "%s/log.txt" % sdir

    # capture output.
    with open(lfile, 'w') as fout:
        #if subprocess.call(cmd, stdout=fout, stderr=fout, cwd=sdir) != 0:
        if subprocess.call(cmd, cwd=sdir) != 0:
            logging.error('error running quast')
            return -1, -1
    '''
    # open the output.
    with open("%s/report.tsv" % sdir) as fin:
        results = fin.readlines()
        tmp = list()
        for result in results:
            tmp.append(result.strip())
        results = tmp

    # extract results.
    lbls = list()
    res1 = list()
    res2 = list()
    cnt =  0
    for x in results:
        tmp = x.split("\t")
        lbls.append(tmp[0])
        res1.append(tmp[1])
        #res2.append(tmp[2])
        res2.append("")
        #print cnt, tmp
        cnt += 1

    res = np.array([lbls,res1,res2])

    return res
