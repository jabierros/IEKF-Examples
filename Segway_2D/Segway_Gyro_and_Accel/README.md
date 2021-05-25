## Segway 2D with gyroscope and acelerometer


Figure depicts the system. `main_symbolic_EKF.m` defines the dynamic equations, and exports required functions
`main_numeric_Information_EKF.m` runs the filter and estimates filter parameters

x_=[q,dq]
dx_/dt=f(x_)=[dq;ddq(q,dq)]

Sensors are a gyro @ C aligned with direction y(xyz)||2(Chassis) and 2 accelerometers @ C aligned along directions 1(Chassis) and 2(Chassis):

z=[dtheta
   ddx(x,dx) cos(theta)-(+g sin(theta))
   ddx(x,dx) sin(theta)-(-g cos(theta))]


![Problem description](https://github.com/jabierros/IEKF-Examples/blob/main/Segway_2D/Segway_Gyro_and_Accel/Segway_2D_figure.png)

