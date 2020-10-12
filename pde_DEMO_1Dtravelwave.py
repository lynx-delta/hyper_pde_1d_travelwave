
import numpy as np
import pde_hyper_1Dtravelwave


def func_uI(x):
    # initial displacement distribution
    y = 1.3 + 2*np.exp(-((x - 15)**2)/15)
    return y


def func_hx(x):
    y = 1*np.exp(-((x - 60)**2)/10)
    return y


def func_hx_1(x):
    y_0 = np.zeros(x.size)
    y_0[((60-8) < x) & (x < (60+8))] = 1
    y = 1*np.cos(np.pi*(x - 60)/(2*8))
    return y*y_0



if __name__ == "__main__":

    import importlib
    importlib.reload(pde_hyper_1Dtravelwave)
    u_I = func_uI
    v_I = 0
    h_sb = func_hx
    L = 100        # length of spatial domain
    T = 52         # total time to solve for
    
    mx = 300       # number of gridpoints in space
    mt = 1000      # number of gridpoints in time
    
    pde_hyper_1Dtravelwave.travelwave_solve(mx, mt, L, T, h_sb, u_I, v_I,
                                            output='plot')
