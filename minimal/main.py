#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

# author: Daniel Kloeser

import time, os
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import cl

# from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
import acados_template as ac
import scipy.linalg
import numpy as np
import types

import casadi as ca


def bicycle_model(s, kappa):
    # define structs
    constraint  = types.SimpleNamespace()
    model       = types.SimpleNamespace()

    model_name = "Spatialbicycle_model"
    length = len(s)
    pathlength = s[-1]

    # compute spline interpolations
    kapparef_s = ca.interpolant("kapparef_s", "bspline", [s], kappa)

    ## Race car parameters
    m = 160.0
    C1 = 0.5
    C2 = 15.5
    Cm1 = 0.28  # power
    Cm2 = 0.0   # powerlim term
    Cr0 = 0.000 # air friction term
    Cr2 = 0.000 # rolling friction term

    ## CasADi Model
    # set up states & controls
    s           = ca.MX.sym("s")
    n           = ca.MX.sym("n")
    alpha       = ca.MX.sym("alpha")
    v           = ca.MX.sym("v")
    D           = ca.MX.sym("D")
    delta       = ca.MX.sym("delta")
    x           = ca.vertcat(s, n, alpha, v, D, delta)

    # controls
    derD        = ca.MX.sym("derD")
    derDelta    = ca.MX.sym("derDelta")
    u           = ca.vertcat(derD, derDelta)

    # xdot
    sdot        = ca.MX.sym("sdot")
    ndot        = ca.MX.sym("ndot")
    alphadot    = ca.MX.sym("alphadot")
    vdot        = ca.MX.sym("vdot")
    Ddot        = ca.MX.sym("Ddot")
    deltadot    = ca.MX.sym("deltadot")
    xdot        = ca.vertcat(sdot, ndot, alphadot, vdot, Ddot, deltadot)

    # algebraic variables
    z = ca.vertcat([])

    # parameters
    p = ca.vertcat([])

    # dynamics
    Fxd = (Cm1 - Cm2 * v) * D - Cr2 * v * v - Cr0 * ca.tanh(5 * v)
    sdota = (v * ca.cos(alpha + C1 * delta)) / (1 - kapparef_s(s) * n)
    f_expl = ca.vertcat(
        sdota,
        v * ca.sin(alpha + C1 * delta),
        v * C2 * delta - kapparef_s(s) * sdota,
        Fxd / m * ca.cos(C1 * delta),
        derD,
        derDelta,
    )

    # constraint on forces
    a_lat = C2 * v * v * delta + Fxd * ca.sin(C1 * delta) / m
    a_long = Fxd / m

    # Model bounds
    model.n_min = -0.12  # width of the track [m]
    model.n_max = 0.12  # width of the track [m]

    # state bounds
    model.throttle_min = -1.0
    model.throttle_max = 1.0

    model.delta_min = -0.40  # minimum steering angle [rad]
    model.delta_max = 0.40  # maximum steering angle [rad]

    # input bounds
    model.ddelta_min = -2.0  # minimum change rate of stering angle [rad/s]
    model.ddelta_max = 2.0  # maximum change rate of steering angle [rad/s]
    model.dthrottle_min = -10  # -10.0  # minimum throttle change rate
    model.dthrottle_max = 10  # 10.0  # maximum throttle change rate

    # nonlinear constraint
    constraint.alat_min = -4  # maximum lateral force [m/s^2]
    constraint.alat_max = 4  # maximum lateral force [m/s^1]

    constraint.along_min = -4  # maximum lateral force [m/s^2]
    constraint.along_max = 4  # maximum lateral force [m/s^2]

    # Define initial conditions
    model.x0 = np.array([-2, 0, 0, 0, 0, 0])

    # define constraints struct
    constraint.alat = ca.Function("a_lat", [x, u], [a_lat])
    constraint.pathlength = pathlength
    constraint.expr = ca.vertcat(a_long, a_lat, n, D, delta)

    # Define model struct
    params = types.SimpleNamespace()
    params.C1 = C1
    params.C2 = C2
    params.Cm1 = Cm1
    params.Cm2 = Cm2
    params.Cr0 = Cr0
    params.Cr2 = Cr2
    model.f_impl_expr = xdot - f_expl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.p = p
    model.name = model_name
    model.params = params

    return model, constraint


def acados_settings(s, kappa):
    Tf          = 1.0       # prediction horizon
    N           = 50        # number of discretization steps
    T           = 10.00     # maximum simulation time[s]
    sref_N      = 3         # reference for final reference progress

    # create render arguments
    ocp = ac.AcadosOcp()

    # export model
    model, constraint = bicycle_model(s, kappa)

    # # define acados ODE
    # model_ac = ac.AcadosModel()
    # model_ac.f_impl_expr = model.f_impl_expr
    # model_ac.f_expl_expr = model.f_expl_expr
    # model_ac.x = model.x
    # model_ac.xdot = model.xdot
    # model_ac.u = model.u
    # model_ac.z = model.z
    # model_ac.p = model.p
    # model_ac.name = model.name
    # ocp.model = model_ac

    # # define constraint
    # model_ac.con_h_expr = constraint.expr

    # # set dimensions
    # nx = model.x.size()[0]
    # nu = model.u.size()[0]
    # ny = nx + nu
    # ny_e = nx

    # ocp.dims.N = N
    # ns = 2
    # nsh = 2

    # # set cost
    # Q = np.diag([ 1e-1, 1e-8, 1e-8, 1e-8, 1e-3, 5e-3 ])

    # R = np.eye(nu)
    # R[0, 0] = 1e-3
    # R[1, 1] = 5e-3

    # Qe = np.diag([ 5e0, 1e1, 1e-8, 1e-8, 5e-3, 2e-3 ])

    # ocp.cost.cost_type = "LINEAR_LS"
    # ocp.cost.cost_type_e = "LINEAR_LS"
    # unscale = N / Tf

    # ocp.cost.W = unscale * scipy.linalg.block_diag(Q, R)
    # ocp.cost.W_e = Qe / unscale

    # Vx = np.zeros((ny, nx))
    # Vx[:nx, :nx] = np.eye(nx)
    # ocp.cost.Vx = Vx

    # Vu = np.zeros((ny, nu))
    # Vu[6, 0] = 1.0
    # Vu[7, 1] = 1.0
    # ocp.cost.Vu = Vu

    # Vx_e = np.zeros((ny_e, nx))
    # Vx_e[:nx, :nx] = np.eye(nx)
    # ocp.cost.Vx_e = Vx_e

    # ocp.cost.zl = 100 * np.ones((ns,))
    # ocp.cost.Zl = 0 * np.ones((ns,))
    # ocp.cost.zu = 100 * np.ones((ns,))
    # ocp.cost.Zu = 0 * np.ones((ns,))

    # # set intial references
    # ocp.cost.yref = np.array([1, 0, 0, 0, 0, 0, 0, 0])
    # ocp.cost.yref_e = np.array([0, 0, 0, 0, 0, 0])

    # # setting constraints
    # ocp.constraints.lbx = np.array([-12])
    # ocp.constraints.ubx = np.array([12])
    # ocp.constraints.idxbx = np.array([1])
    # ocp.constraints.lbu = np.array([model.dthrottle_min, model.ddelta_min])
    # ocp.constraints.ubu = np.array([model.dthrottle_max, model.ddelta_max])
    # ocp.constraints.idxbu = np.array([0, 1])
    # # ocp.constraints.lsbx=np.zero s([1])
    # # ocp.constraints.usbx=np.zeros([1])
    # # ocp.constraints.idxsbx=np.array([1])
    # ocp.constraints.lh = np.array(
    #     [
    #         constraint.along_min,
    #         constraint.alat_min,
    #         model.n_min,
    #         model.throttle_min,
    #         model.delta_min,
    #     ]
    # )
    # ocp.constraints.uh = np.array(
    #     [
    #         constraint.along_max,
    #         constraint.alat_max,
    #         model.n_max,
    #         model.throttle_max,
    #         model.delta_max,
    #     ]
    # )
    # ocp.constraints.lsh = np.zeros(nsh)
    # ocp.constraints.ush = np.zeros(nsh)
    # ocp.constraints.idxsh = np.array([0, 2])

    # # set intial condition
    # ocp.constraints.x0 = model.x0

    # # set QP solver and integration
    # ocp.solver_options.tf = Tf
    # # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    # ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    # ocp.solver_options.nlp_solver_type = "SQP_RTI"
    # ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    # ocp.solver_options.integrator_type = "ERK"
    # ocp.solver_options.sim_method_num_stages = 4
    # ocp.solver_options.sim_method_num_steps = 3

    # # ocp.solver_options.qp_solver_tol_stat = 1e-2
    # # ocp.solver_options.qp_solver_tol_eq = 1e-2
    # # ocp.solver_options.qp_solver_tol_ineq = 1e-2
    # # ocp.solver_options.qp_solver_tol_comp = 1e-2

    # # create solver
    # acados_solver = ac.AcadosOcpSolver(ocp, json_file="acados_ocp.json")

    # return constraint, model, acados_solver

def generate_straight():
    n = 1000
    L = 100
    r = 1
    width = 3
    w = width*np.ones(n)
    l = np.linspace(0, L, n)
    p = np.array([l, np.zeros(l.shape)]).T
    return p, w

def generate_circle():
    n = 100
    r = 7.5
    width = 3
    th = np.linspace(0, 2*np.pi, n)
    w = width*np.ones(n)
    p = r*np.array([np.cos(th), np.sin(th)]).T
    return p, w



# p,w = generate_circle()
p,w = generate_straight()
reftrack = np.column_stack(( p, w, w ))
cl = cl.Centerline(reftrack)
s, reftrack, kappa, normvec, psi = cl.discretize(0,cl.end(), 100)


if False:
    plt.subplot(221)
    plt.plot(s, label="s")
    plt.legend()
    plt.subplot(222)
    plt.scatter(reftrack[:,0], reftrack[:,1], label="reftrack")
    plt.legend()
    plt.subplot(223)
    plt.plot(kappa, label="kappa")
    plt.legend()
    plt.subplot(224)
    plt.plot(psi, label="psi")
    plt.legend()
    plt.show()

# # load model
acados_settings(s, kappa)

# # dimensions
# nx = model.x.size()[0]
# nu = model.u.size()[0]
# ny = nx + nu
# Nsim = int(T * N / Tf)

# # initialize data structs
# simX = np.ndarray((Nsim, nx))
# simU = np.ndarray((Nsim, nu))
# s0 = model.x0[0]
# tcomp_sum = 0
# tcomp_max = 0

# # simulate
# for i in range(Nsim):
#     # update reference
#     sref = s0 + sref_N
#     for j in range(N):
#         yref = np.array([s0 + (sref - s0) * j / N, 0, 0, 0, 0, 0, 0, 0])
#         # yref=np.array([1,0,0,1,0,0,0,0])
#         acados_solver.set(j, "yref", yref)
#     yref_N = np.array([sref, 0, 0, 0, 0, 0])
#     # yref_N=np.array([0,0,0,0,0,0])
#     acados_solver.set(N, "yref", yref_N)

#     # solve ocp
#     t = time.time()

#     status = acados_solver.solve()
#     if status != 0:
#         print("acados returned status {} in closed loop iteration {}/{}.".format(status, i, Nsim))

#     elapsed = time.time() - t

#     # manage timings
#     tcomp_sum += elapsed
#     if elapsed > tcomp_max:
#         tcomp_max = elapsed

#     # get solution
#     x0 = acados_solver.get(0, "x")
#     u0 = acados_solver.get(0, "u")
#     for j in range(nx):
#         simX[i, j] = x0[j]
#     for j in range(nu):
#         simU[i, j] = u0[j]

#     # update initial condition
#     x0 = acados_solver.get(1, "x")
#     acados_solver.set(0, "lbx", x0)
#     acados_solver.set(0, "ubx", x0)
#     s0 = x0[0]

#     # check if one lap is done and break and remove entries beyond
#     if x0[0] > Sref[-1] + 0.1:
#         # find where vehicle first crosses start line
#         N0 = np.where(np.diff(np.sign(simX[:, 0])))[0][0]
#         Nsim = i - N0  # correct to final number of simulation steps for plotting
#         simX = simX[N0:i, :]
#         simU = simU[N0:i, :]
#         break

