import time, os
import numpy as np
from acados_settings_dev import *
from plotFcn import *
from tracks.readDataFcn import getTrack
import matplotlib.pyplot as plt

import centerline
import numpy as np


# get coords and track widths out of array
def import_track(file_path: str) -> np.ndarray:

    # load data from csv file
    csv_data_temp = np.loadtxt(file_path, delimiter=';')

    # get coords and track widths out of array
    refline = csv_data_temp[:, 0:2]
    w_tr_r = csv_data_temp[:, 2]
    w_tr_l = csv_data_temp[:, 3]

    # assemble to a single array
    reftrack_imp = np.column_stack((refline, w_tr_r, w_tr_l))

    return reftrack_imp

def main():
    track_imp = import_track("trackdrive.csv")
    cl = centerline.Centerline(track_imp)
    reftrack, kappa, normvec = cl.discretize(0, cl.end(), 100)

    """
    Example of the frc_racecars in simulation without obstacle avoidance:
    This example is for the optimal racing of the frc race cars.
    The model is a simple bicycle model and the lateral acceleration is
    constraint in order to validate the model assumptions.
    The simulation starts at s=-2m until one round is completed(s=8.71m).
    The beginning is cut in the final plots to simulate a 'warm start'.
    """

    track = "LMS_Track.txt"
    [Sref, _, _, _, _] = getTrack(track)

    Tf = 1.0  # prediction horizon
    N = 50  # number of discretization steps
    T = 10.00  # maximum simulation time[s]
    sref_N = 3  # reference for final reference progress


    # load model
    constraint, model, acados_solver = acados_settings(Tf, N, track)



main()
