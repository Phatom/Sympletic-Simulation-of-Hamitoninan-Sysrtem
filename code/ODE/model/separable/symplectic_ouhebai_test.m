% Use the symp4_Neri to solve ouhebai
% Only curve in phase space are plotted



cc
n = 2;
t_max = 10000;
time_steps=100000;
tspan = linspace(0,t_max,time_steps);


% I find it's hard to use symbolic computation to give something like dqdt = @(q)[q(1),q(2)];
% m = 1;
% k = 1;
% l = 1;
% g = 1;
% 
% disp('calculating diravitaves')
% tic
% syms p1 q1 p2 q2
% H = @(p1,q1,p2,q2)((p1.^2+p2.^2)/(2*m*l^2)+k*l^2/2*(q2-q1).^2+m*g*l/2*(q1.^2+q2.^2));
% diff_p = [diff(H,p1),diff(H,p2)];
% dqdt = matlabFunction(diff_p);
% diff_q = [diff(H,q1),diff(H,q2)];
% dpdt = matlabFunction(diff_q);
% toc

dqdt = @(p)[p(1),p(2)];
dpdt = @(q)[2*q(1)-q(2),-q(1)+2*q(2)];

q10 = 1;
q20 = 1;
p10 = 1;
p20 = 1;

fprintf('symp4_Neri:\n');
tic
[se4.q,se4.p] = symp4_Neri(dqdt,dpdt,tspan,[q10,q20],[p10,p20]);
toc

figure
subplot(1,2,1)
plot(se4.q(:,1),se4.p(:,1))
subplot(1,2,2)
plot(se4.q(:,2),se4.p(:,2))
