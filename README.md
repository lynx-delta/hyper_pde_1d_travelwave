## hyper_pde_1d_travelwave

Hyperbolic partial differential equation solver: visualization of 1-dimensional travelling wave (for teaching purposes)



### Wave Equation with Non-constant Wavespeed

The scripts in this repository are used to solve and visualize the hyperbolic Partial Differential Equation (PDE) governing the motion of a tsunami (or similar) in the open
ocean, in the case of a variable height seabed:  

![Wave equation with non-constant wavespeed](data/formula_1.png)  

with *u(x, t)* the non-dimensional sea surface height and *h(x)* the non-dimensional still-water height (see figure below).  

![Sea surface and still-water height](data/figure_1.png)  

The initiation of the tsunami (e.g. a subsea earthquake) can be modelled by the initial displacement  

![Initiation of tsunami](data/formula_2.png)  

centered at *x*<sub>1</sub>, with (additional) sea surface height <i>&alpha;</i><sub>1</sub> and spread <i>&sigma;</i><sub>1</sub>. The subsea hill (or depression) is modelled as  

![Subsea hill](data/formula_3.png)  

with the center at *x*<sub>B</sub>, the elevation &alpha<sub>B</sub> and the spread &sigma<sub>B</sub>. To model infinite spatial domains, an open boundary at *x* = 0 respectively at *x* = L given by the condition  

![Open boundary](data/formula_4.png)  

is implemented.  



### Solution and Visualization

The hyperbolic PDE is solved employing finite difference methods. To model, solve and visualize the hyperbolic PDE for time *t*, use the files `pde_hyper_1Dtravelwave.py` and `pde_DEMO_1Dtravelwave.py`.  

![Solution for time t](data/figure_2.png)  

For an animated solution, use file `pde_hyper_1Dtravelwave_anim.py` (`pde_hyper_1Dtravelwave.py` is imported and used for solving the PDE for a certain time frame).


