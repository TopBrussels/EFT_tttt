import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM

import os
import sys
import time
import argparse
import logging
import json
from pprint import pprint

def progress(current, total, status=''):
    fullBarLength = 80
    doneBarLength = int(round(fullBarLength * current / float(total)))

    percents = round(100.0 * current / float(total), 1)
    bar = '>' * doneBarLength + ' ' * (fullBarLength - doneBarLength)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()

sig_SM = 0.009308 * 1000. # fb from MG

def main():
    np.set_printoptions(edgeitems=3)
    np.core.arrayprint._line_width = 120

    eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)

    xs_limit = 0.0177+sig_SM/1000.

    min_val, max_val = -15., 15
    marginalized_step = 3
    wilson_of_interest_step = 0.1
    marginal_min_wilson = +100000000.
    marginal_max_wilson = -100000000.

    total = ((max_val - min_val)/marginalized_step)**4 * (max_val - min_val)/wilson_of_interest_step
    print "Total number of evaluations: {}".format(total)
    i = 0
    for C1 in np.arange(min_val, max_val, marginalized_step):
        for C2 in np.arange(min_val, max_val, marginalized_step):
            for C3 in np.arange(min_val, max_val, marginalized_step):
                for C4 in np.arange(min_val, max_val, marginalized_step):
                    progress(i, total, 'Progress'); i += ((max_val - min_val)/wilson_of_interest_step)
                    for C5 in np.arange(min_val, max_val, wilson_of_interest_step):
                        xs = eft.gen_eft_xs([C5, C1, C2, C3, C4])
                        if xs < xs_limit:
                            if xs > marginal_max_wilson:
                                marginal_max_wilson = xs
                            if xs < marginal_min_wilson:
                                marginal_min_wilson = xs

    print '\nMarginal limits: [{}; {}]'.format(marginal_min_wilson, marginal_max_wilson)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    start_time = time.time()
    logging.info(time.asctime())
    exitcode = main()
    logging.info(time.asctime())
    logging.info('TOTAL TIME IN MINUTES:' + str((time.time() - start_time) / 60.0))