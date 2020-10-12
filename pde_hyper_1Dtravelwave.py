
import numpy as np
import matplotlib.pyplot as plt


def travelwave_solve(mx, mt, L, T, h_sb, u_I, v_I, output='plot'):
    '''
    Solver for 1D wave equation with variable wave speed 
    using finite difference method (especially to solve the motion
    of a tsunami in the open ocean)
    
    param mx: number of gridpoints in space
              integer object
    param mt: number of gridpoints in time
              integer object
    param L: length of spatial domain
             integer or float object
    param T: total time to solve for
             integer or float object
    param h_sb: seabed shape / height along x
               integer object, float object or function
               e.g. function
               def u_I(x):
                   y = 1*np.exp(-((x - 50)**2)/10)
                   return y
    param u_I: initial displacement distribution
               integer object, float object or function
               e.g. function
               def u_I(x):
                   y = 1.5 + 2*np.exp(-((x - 8)**2)/12)
                   return y
    param v_I: initial velocity distribution
               integer object, float object or function
               e.g. function
               def u_I(x):
                   y = np.sin(pi*x)
                   return y
    param output: output='data',
                  x: array with values of x (gridpoints in space)
                  h_sb: array with seabed shape along x
                  u_jpls: array with solution of u(x, T)
                  output='plot',
                  plot u_jpls against x
                  default = 'plot'
    return: depends on the option --> see param output
    '''
    
    # Set up the numerical environment variables
    x = np.linspace(0, L, mx+1)    # mesh points in space
    t = np.linspace(0, T, mt+1)    # mesh points in time
    dx = x[1] - x[0]               # gridspacing in x
    dt = t[1] - t[0]               # gridspacing in t
    
    # Initialise initial condition
    if isinstance(h_sb, (int, float)):
        h_sb = np.full((1, mx+1), h_sb)[0]  # constant for all x
    elif callable(h_sb):
        h_sb = h_sb(x)
    if isinstance(u_I, (int, float)):
        u_I = np.full((1, mx+1), u_I)[0]  # constant for all x
    elif callable(u_I):
        u_I = u_I(x)
    if isinstance(v_I, (int, float)):
        v_I = np.full((1, mx+1), v_I)[0]  # constant for all x
    elif callable(v_I):
        v_I = v_I(x)
    
    # Define h_x using seabed shape h_sb and h_0
    h_0 = min(u_I)
    h_x = h_0 - h_sb
    h_0max = max(u_I)
    # Make sure that there is enough water in the sea
    if max(h_sb) >= h_0:
        print('There is not enough water in the sea for a 1D simulation!')
        return
    # Courant number (stability check)
    hmax = max(h_x)
    if hmax*dt/dx > 1:                 
        mt = np.ceil(T/(dx/hmax))      # Change dt so that lambda is
        t = np.linspace(0, T, mt+1)    # smaller or equal to 1.0
        dt = t[1] - t[0]
        
    # Set mesh constant
    Csqr = (dt/dx)**2
    # Set up the solution variables
    u_j = np.zeros(x.size)             # u at current time step
    u_jp1s = np.zeros(x.size)          # u at next time step
    # Set initial condition
    u_jmns = u_I                # u at previous time step / initial condition
    
    # First timestep
    u_j[1:-1] = u_jmns[1:-1] + dt*v_I[1:-1] + 0.5*Csqr*(0.5*(h_x[1:-1] \
                + h_x[2:])*(u_jmns[2:] - u_jmns[1:-1]) - 0.5*(h_x[1:-1] \
                + h_x[:-2])*(u_jmns[1:-1] - u_jmns[:-2]))
    # Open boundary at x = 0
    u_j[0] = u_jmns[0] + h_x[0]*(dt/dx)*(u_jmns[1] - u_jmns[0])
    # Open boundary at x = L
    u_j[-1] = u_jmns[-1] - h_x[-1]*(dt/dx)*(u_jmns[-1] - u_jmns[-2])
    
    # Loop for regular timesteps
    for n in range(2, mt+1):
        u_jp1s[1:-1] = -u_jmns[1:-1] + 2*u_j[1:-1] + Csqr*(0.5*(h_x[1:-1] \
                      + h_x[2:])*(u_j[2:] - u_j[1:-1]) - 0.5*(h_x[1:-1] \
                      + h_x[:-2])*(u_j[1:-1] - u_j[:-2]))
        # Open boundary at x = 0
        u_jp1s[0] = u_j[0] + h_x[0]*(dt/dx)*(u_j[1] - u_j[0])
        # Open boundary at x = L
        u_jp1s[-1] = u_j[-1] - h_x[-1]*(dt/dx)*(u_j[-1] - u_j[-2])
        # Update arrays for next timestep
        # All 3 u_j to reallocate memory 'pointers'
        u_jmns, u_j, u_jp1s = u_j, u_jp1s, u_jmns
    
    # Output
    if output in ('PLOT', 'Plot', 'plot'):
        plot_1Dtravelwave(x, u_jp1s, h_sb, h_0max, T)
    elif output in ('DATA', 'Data', 'data'):
        return x, h_sb, u_jp1s
    else:
        print("Check input argument 'output'!")
        return

      
def plot_1Dtravelwave(x, u_jpls, h_sb, h_0max, T):   
    '''
    Plotting function: plots u_jpls (solution) against x
                       including h_sb (seabed shape)
    '''
    
    # Plot
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(1, 1, 1)
    label_1 = '$seabed$'
    label_2 = '$u(x, ' + str(T) + ')$'
    ax1.plot(x, h_sb, 'bo', markersize=2, label=label_1)
    ax1.plot(x, u_jpls, color='red', linewidth=1.5, label=label_2)
    # Limits for x-axis
    xrange = abs((max(x)-min(x)))
    xmin = min(x) - 0.04*xrange
    xmax = max(x) + 0.04*xrange
    # Limits for y-axis
    yrange = abs(h_0max - min(h_sb))
    ymin = min(h_sb) - 0.08*yrange
    ymax = h_0max + 0.08*yrange
    # Set axes limits
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    ax1.grid(True)
    # Add legend
    legend = ax1.legend(loc='best', fontsize='medium')
    legend.get_frame().set_alpha(1)
    # Use tight layout
    plt.tight_layout()
    plt.show()   



if __name__ == "__main__":
    import sys
    travelwave_solve(*sys.argv[1:])

