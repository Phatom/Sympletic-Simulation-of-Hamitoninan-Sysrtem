% Compare ode45 with different order symplectic integrator for harmonic oscillator
% curve in phase space are plotted
% Yan Zhang & Yaoqi Zhang, 10-July-2018

cc
n=1;
t_max = 10000;
time_steps=100000;
tspan = linspace(0,t_max,time_steps);

disp('calculating diravitaves')
tic
syms  p q
H = @(p,q)(p.^2/2-cos(q));
diff_p = diff(H,p);
dqdt = str2func(['@(p)',vectorize(diff_p)]);
diff_q = diff(H,q);
dpdt = str2func(['@(q)',vectorize(diff_q)]);
toc

q0=0;
p0=2;

fprintf('symp4_Neri:\n');
tic
[se4.q,se4.p] = symp4_Neri(dqdt,dpdt,tspan,q0,p0);
toc

% ode45
fprintf('ode45:\n');
% Setup
odefun=@(t,y)([y(2);-sin(y(1))]);
y0=[q0;p0];
% Solver options
opt=odeset('Stats','on');
% Solve
tic
[~,y]=ode45(odefun,tspan,y0,opt); 
toc 
ode45ns.q=y(:,1);
ode45ns.p=y(:,2);

[SN,~,~] = ellipj(tspan,1);
q_theory = SN;

figure

subplot(1,2,1)
plot(ode45ns.q,ode45ns.p,se4.q,se4.p); % plot
xlabel('q')
ylabel('p')
title('Pendulum:phase space(long time)')
legend('ode45','symp4')
subplot(1,2,2)
plot(se4.q,se4.p); % plot
title('Pendulum:phase space(long time)')
xlabel('q')
ylabel('p')
axis equal