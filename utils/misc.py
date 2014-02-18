import subprocess
import os
import numpy as np

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


def runit(job, path):
    """ runs the job """
    with open(path, 'wb') as fout:
        fout.write('\n'.join(job))
    subprocess.call(['chmod', 'a+x', path])
    subprocess.call(['bash', path])

def create_dir(dir_path):
    ''' creates directory if necessary'''
    if os.path.isdir(dir_path) == False:
        if subprocess.call(["mkdir", dir_path]) != 0:
            logging.error("couldn't make dir")
            sys.exit(1)



def save_agps(agp_file, agp):
    ''' saves agp to disk.'''

    # write to file.
    fout = open(agp_file, "w")

    # write each entry.
    z = len(agp_dt.names)
    for i in range(agp.size):

        # sanity skip.
        if agp[i]['scaf_name'] == "":
            continue

        # format result.
        tmp = agp[i]
        if tmp['comp_type'] == "W":
            # get orientation.
            if tmp["comp_orien"] == 0:
                o = "+"
            else:
                o = "-"

            # write contig.
            txt = str(tmp['scaf_name']) + "\t"
            txt += str(tmp['scaf_start']) + "\t"
            txt += str(tmp['scaf_stop']) + "\t"
            txt += str(tmp['scaf_idx']) + "\t"
            txt += str(tmp['comp_type']) + "\t"
            txt += str(tmp['comp_name']) + "\t"
            txt += str(tmp['comp_start']) + "\t"
            txt += str(tmp['comp_stop']) + "\t"
            txt += o + "\n"

        else:
            # get linkage.
            if tmp['comp_linkage'] == 0:
                o = "no"
            else:
                o = "yes"

            # write gap.
            txt = str(tmp['scaf_name']) + "\t"
            txt += str(tmp['scaf_start']) + "\t"
            txt += str(tmp['scaf_stop']) + "\t"
            txt += str(tmp['scaf_idx']) + "\t"
            txt += str(tmp['comp_type']) + "\t"
            txt += str(tmp['comp_stop'] - tmp['comp_start']) + "\t"
            txt += str(tmp['comp_name']) + "\t"
            txt += o + "\n"

        fout.write(txt)

    # close file.
    fout.close()
