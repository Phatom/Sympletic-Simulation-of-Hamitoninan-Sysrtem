%计算每次检验误差时用到的delta，即步长减半后计算的函数值与减半前计算的函数值的差
function [delta,q,p]=solve_delta(dqdt,dpdt,q0,p0,steps)
tspan=0:steps:steps;
[q1,p1]=symp4_Neri(dqdt,dpdt,tspan,q0,p0);
tspan=0:steps/2:steps;
[q2,p2]=symp4_Neri(dqdt,dpdt,tspan,q0,p0);
q=q2(3);p=p2(3);
delta=max(abs(q2(3)-q1(2)),abs(p2(3)-p1(2)));
end