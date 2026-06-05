## Mass Spring Damper
The standard mass spring damper:

m ddx + c dx + k (x-rho0)=f_ext

<img src="https://render.githubusercontent.com/render/math?math=m%20%5Cddot%7Bx%7D%20%2Bk(x-%5Crho_0)%20%2B%20c%20%5Cdot%7Bx%7D%20%3D%20f_%7Bext%7D">

Sensor is an accelerometer:

z=ddx(x,dx)

<img src="https://render.githubusercontent.com/render/math?math=z%3D%20%5Cddot%7Bx%7D(x%2C%5Cdot%7Bx%7D)">

![Problem description](https://github.com/jabierros/IEKF-Examples/blob/main/Mass_Spring_Damper/mass_spring_damper.png)

This example utilizes decoupled configuration and state structures (`KF`, `TrueSystem`, `SimOpts`) to simulate the system and run the filter.

### Code Files (LibIEKF_0.1)
- [main_symbolic_EKF.m](https://github.com/jabierros/IEKF-Examples/blob/LibIEKF_0.1/Mass_Spring_Damper/main_symbolic_EKF.m) (Symbolic equations generation)
- [main_numeric_Information_EKF.m](https://github.com/jabierros/IEKF-Examples/blob/LibIEKF_0.1/Mass_Spring_Damper/main_numeric_Information_EKF.m) (Numerical IEKF simulation & optimization loop)

Refs: help to MathML rendering in github https://jsfiddle.net/8ndx694g/

---
*The links below correspond to the older version (`LibIEKF_0`):*

See OLD SYMBOLIC mlx file: 

https://htmlpreview.github.io/?https://github.com/jabierros/IEKF-Examples/blob/LibIEKF_0/Mass_Spring_Damper/main_symbolic_EKF_mlx.html

See OLD NUMERIC mlx file: 

https://htmlpreview.github.io/?https://github.com/jabierros/IEKF-Examples/blob/LibIEKF_0/Mass_Spring_Damper/main_numeric_Information_EKF_mlx.html
