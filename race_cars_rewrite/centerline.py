import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

class BSpline:

    """
    Wrapper of scipy spline libary
    """
    def __init__(self, x, t_length, closed=True):
        t = np.linspace(0, t_length, num=x.size, endpoint=True)
        that, c, k = interpolate.splrep(t, x, s=0, k=4)
        self.spline = interpolate.BSpline(that, c, k, extrapolate=True)
        self.tmin = t.min()
        self.tmax= t.max()
        self.closed = closed

    def value(self, t):
        if isinstance(t, (list, tuple, np.ndarray)):
            return self.value_array(t)

        # if self.closed:
        #     # This normalization can cause wierd numerical issues since Bspline is
        #     # not periodic.
        #     while t > self.tmax:
        #         t -= self.tmax
        #     while t < 10**-10:
        #         t += self.tmax
        # else:
        #     if t > self.tmax:
        #         t = self.tmax
        #     if t < self.tmin:
        #         t =  self.tin

        return self.spline(t)

    def value_array(self, t):
        spline_values = np.zeros(t.shape)
        for i in range(t.size):
            spline_values[i] = self.value(t[i])
        return spline_values

class BSpline2D:
    def __init__(self, X, t_length, closed=True):
        self.x_spline = BSpline(X[:,0], t_length, closed=True )
        self.y_spline = BSpline(X[:,1], t_length, closed=True )
        return

    def value(self, t):
        x = self.x_spline.value(t)
        y = self.y_spline.value(t)
        p = np.array([x,y]).T
        return p


class Centerline:

    def __init__(self, reftrack):

        # Parameterize everything on t
        s = self.compute_length(reftrack[:,:2], exact=True)
        self.s_length       = s[-1]
        self.s_spline       = BSpline(s, self.s_length, closed=True)
        self.p_spline       = BSpline2D(reftrack[:,:2], self.s_length, closed=True )
        self.w_spline       = BSpline2D(reftrack[:, 2:], self.s_length, closed=True )

        def compute_kappa(p):
            def det(v,q):
                # z axis of cross product
                return v[0]*q[1] - q[0]*v[1]
            def curvature(p1,p2,p3):
                # menger curvature
                return 2*det(p2 - p1, p3 - p1)/(np.linalg.norm(p2-p3)*np.linalg.norm(p2-p3)*np.linalg.norm(p1-p3))
            # Evalulating kappa, inclusive start and end point
            kappa = np.zeros(p.shape[0])
            kappa[0] = curvature(p[-1], p[0], p[1])
            kappa[-1] = curvature(p[-1], p[0], p[1])
            for i in range(1, p.shape[0] -1):
                kappa[i] = curvature(p[i-1], p[i], p[i+1])
            return kappa

        # discretize centerline
        n_points = 100
        t = np.linspace(0, self.s_length, n_points, endpoint=True)
        p = self.p_spline.value(t)

        # Make spline of curvature
        kappa = compute_kappa(p)
        self.kappa_spline = BSpline( kappa, self.s_length, closed=True )

        # Make spline of normal vector spline
        R = np.array([[0, -1],
                      [1, 0]])
        T = p[1:] - p[:-1]
        N = T@R.T
        N /= np.linalg.norm(N, ord=1, axis=1, keepdims=True)
        self.N_spline       = BSpline2D(N, self.s_length, closed=True )



    def compute_length(self, wpts, exact=False, N_exact=1000):

        if exact:
            # inexact differ approx 0.16% from exact length
            x_fine_spline = BSpline(wpts[:,0], t_length=1, closed=True )
            y_fine_spline = BSpline(wpts[:,1], t_length=1, closed=True )
            t = np.linspace(0,1, N_exact)
            x = x_fine_spline.value(t)
            y = y_fine_spline.value(t)
            wpts = np.array([x,y]).T

        dp = wpts[1:] - wpts[:-1]
        l = np.sqrt(dp[:,0]**2 + dp[:,1]**2)
        s = np.cumsum(l)
        s = np.append([0],s)
        return s

    def proj(self, p, psi, t0 = None, method="bruteforce"):
        # Quite shitty method
        tt = None
        N, dt = 1000, 10
        if method=="bruteforce":
            if t0 is not None:
                tt = np.linspace(t0 - dt, t0 + dt, N)
            else:
                tt = np.linspace(0, self.t_length, N)

            e = np.inf
            tproj = 0
            for z in tt:
                pbar = self.p_spline.value(z)
                d = np.linalg.norm(p - pbar)
                if d < e:
                    e = d
                    tproj = z
        return tproj

    def discretize(self, t0, t1, N):
        t = np.linspace(t0, t1, N)
        reftrack        = np.zeros((N,4))
        reftrack[:, :2] = self.p_spline.value(t)
        reftrack[:, 2:] = self.w_spline.value(t)
        kappa           = self.kappa_spline.value(t)
        normvec         = self.N_spline.value(t)
        s               = self.s_spline.value(t)
        return s, reftrack, kappa, normvec


    def end(self):
        return self.s_length


if False:
    n = 10
    a = np.zeros((n,2))
    b = np.linspace(0,2*np.pi, n)
    a = 10*np.array([np.cos(b), np.sin(b)]).T
    plt.plot(a[:,0], a[:,1], label="blue")

    b = np.linspace(0,2*np.pi, 30)
    a_spline = BSpline2D(a, b[-1])
    p = a_spline.value(b)
    plt.plot(p[:,0], p[:,1], label="red")
    plt.legend()
    plt.show()
