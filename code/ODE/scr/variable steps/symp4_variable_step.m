%采用变步长的四阶辛算法解ode的示例，在四阶辛算法的基础上加入了每一步计算前判断步长的部分，
%输出setq，setp为方程求解的结果，setsteps为每一步的步长，在方程求解精度要求较高的情况下，我们可以看到步长变化很明显。
function [setq,setp,setsteps] = symp4_variable_step(dqdt,dpdt,steps,t_max,q0,p0)%steps为初始步长,t_max为终止计算点
t=0; i=1;
tic
while t < t_max
    [steps,q,p] = solve_steps(dqdt,dpdt,p0,q0,steps);
    p0=p;q0=q;
    setq(i)=q; setp(i)=p;
    setsteps(i) = steps;
    t=t+steps;
    i = i+1;
end
toc

% figure
% plot(setq(1:2),setp(1:2))
% hold on
% for i=2:10000
%     plot(setq(i:(i+1)),setp(i:(i+1)));
% end