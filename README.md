# Mass Spring Damper IEKF

This is an example of a, pretty general, **MATLAB** implememntation of the Information Extended Kalman Filter.
Maximum Likelihood Estimation of filter parameters is considered.

Real world data is generated algorithmicaly. To that end an "actual system" state is integrated in paralled with the IEKF. From this "actual" state, "actual" meassurements (noise free)  are obtained. These, in turn, are contaminated with noise. Running in parallel with the IEKF allow to use the filter output as a feedback for control of the system. 

The code is generic enough to be applied to general nonlinear proccess and sensor equations

## Globals
Extensive use of global variables is used on purpose to keep the code as simple as posible.


function handle `u_actual_func = @(t) (...)` should be defined so that it uses mu_x global variable to generate the output.

## The example (mass spring damper)
The standard mass spring damper:

m ddx + c dx + k (x-rho0)=f_ext

<img src="https://render.githubusercontent.com/render/math?math=m%20%5Cddot%7Bx%7D%20%2Bk(x-%5Crho_0)%20%2B%20c%20%5Cdot%7Bx%7D%20%3D%20f_%7Bext%7D">

Sensor is an accelerometer:

z=ddx(x,dx)

<img src="https://render.githubusercontent.com/render/math?math=z%3D%20%5Cddot%7Bx%7D(x%2C%5Cdot%7Bx%7D)">

![Problem description](https://github.com/jabierros/Mass-Spring-Damper-IEKF/blob/main/mass_spring_damper.png)


