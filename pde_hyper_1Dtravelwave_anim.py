
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import pde_hyper_1Dtravelwave


def func_uI(x):
    # initial displacement distribution
    y = 1.0 + 2*np.exp(-((x - 15)**2)/15)
    return y


def func_hx(x):
    y = 0.7*np.exp(-((x - 60)**2)/10)
    return y


def func_hx_1(x):
    y_0 = np.zeros(x.size)
    y_0[((60-8) < x) & (x < (60+8))] = 1
    y = 0.7*np.cos(np.pi*(x - 60)/(2*8))
    return y*y_0


def animate_update(i, line, mx, mt, L, h_sb, u_I, v_I):
    x, h_sb, u_jpls = pde_hyper_1Dtravelwave.travelwave_solve(mx, mt, L, i,
                                                              h_sb, u_I, v_I,
                                                              'data')
    line.set_ydata(u_jpls)  # update the data
    return line,


u_I = func_uI
v_I = 0
h_sb = func_hx
L = 100            # length of spatial domain
mx = 400           # number of gridpoints in space
mt = 1000          # number of gridpoints in time

# Initialise data for plot
x, h_sb, u_jpls = pde_hyper_1Dtravelwave.travelwave_solve(mx, mt, L, 0,
                                                              h_sb,u_I, v_I,
                                                              'data')
# Initialise plot
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(x, h_sb, 'bo', markersize=2)
line = ax1.plot(x, u_jpls, color='red', linewidth=1.5)[0]
    
# Additional parameters to pass to function
param = (line, mx, mt, L, h_sb, u_I, v_I)
    
# Animate figure
anm = anim.FuncAnimation(fig, animate_update, frames=100,
                   fargs=param, interval=150, blit=True)
plt.show()




