"""
  Reads the output of messpf
"""


def read_pf(pf_file_name):
    """ reads the pf.dat file.
    """

    temps, logq, dq_dt, dq2_dt2 = [], [], [], []
    with open(pf_file_name, 'r') as pf_file:
        for i, line in enumerate(pf_file):
            if i != 0:
                tmp = line.strip().split()
                temps.append(tmp[0])
                logq.append(tmp[1])
                dq_dt.append(tmp[2])
                dq2_dt2.append(tmp[3])

    return temps, logq, dq_dt, dq2_dt2
