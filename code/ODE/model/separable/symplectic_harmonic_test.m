% Compare ode45 with different order symplectic integrator for harmonic oscillator
% curve in phase space, q-t curve and energy error are plotted.

% Yan Zhang & Yaoqi Zhang, 10-July-2018
clear all
close all
clc


n=1;
t_max = 10000;
time_steps=1000000;

% Dynamic Setup ?????
% ?????????????????windows?matlab2017a???????????
% windows????2018a????mac????2018a???
% ????????????????????????????
disp('calculating diravitaves')
tic
syms H(p,q) p q
H = @(p,q)(p.^2+q.^2)/2;
diff_p = diff(H,p);
dqdt = matlabFunction(diff_p);
diff_q = diff(H,q);
dpdt = matlabFuntcion(diff_q);
toc

% % Dynamic Setup
% disp('calculating diravitaves')
% tic
% syms H(p,q) p q
% H = @(p,q)(p.^2+q.^2)/2;
% diff_p = diff(H,p);
% dqdt = str2func(['@(p)',vectorize(diff_p)]);
% diff_q = diff(H,q);
% dpdt = str2func(['@(q)',vectorize(diff_q)]);
% toc

% Initial SETUP
q0=0.5;
p0=0.5;

% Time interval
tspan = linspace(0,t_max,time_steps);


% SYMPLECTIC SOLVER
% order1: Symmplectic Euler
fprintf('symp1_Euler:\n');
tic
[se.q,se.p]=symp1_Euler(dqdt,dpdt,tspan,q0,p0); % symp1_Euler is in the package [1]
toc
% order2: St?rmer–Verlet
fprintf('symp2_Verlet:\n');
tic
[se2.q,se2.p] = symp2_Verlet(dqdt,dpdt,tspan,q0,p0);
toc
% order4: Neri BCH
fprintf('symp4_Neri:\n');
tic
[se4.q,se4.p] = symp4_Neri(dqdt,dpdt,tspan,q0,p0);
toc


% ode45
fprintf('ode45:\n');
% Setup
odefun=@(t,y)([y(2);-y(1)]);
y0=[q0;p0];
% Solver options
opt=odeset('Stats','on');
% Solve
tic
[~,y]=ode45(odefun,tspan,y0,opt); 
toc 
ode45ns.q=y(:,1);
ode45ns.p=y(:,2);

 
% Test and PLOT

% phase space
figure
set(gca,'nextplot','replacechildren'); % axis settings
psf=plot(ode45ns.q,ode45ns.p,se.q,se.p,se2.q,se2.p,se4.q,se4.p); % plot
psf(1).Color = 'r';
psf(2).Color = 'b';
psf(3).Color = 'g';
psf(4).Color = 'k';
axis square;
axis([-2 2 -2 2]);
title('Harmonic Oscillator: Phase Space(1)')
xlabel(gca,'q (= y)');
ylabel('p (= y'')');
hold on
% Theory
theta = 0:0.05*pi:2*pi;
xs = sqrt(0.5) * cos(theta);
ys = sqrt(0.5) * sin(theta);
plot(xs,ys)
legend('Symp-1','Symp-2','Symp-4','ode45','Theory')
text(1.5,-1.9,sprintf('t:%0.1f',tspan(end)));

figure
plot(se.q,se.p,se2.q,se2.p,se4.q,se4.p,xs,ys);
legend('Symp-1','Symp-2','Symp-4','Theory');
title('Harmonic Oscillator: Phase Space(2)')
xlabel(gca,'q (= y)');
ylabel('p (= y'')');

% q-t、p-t图（）
q_theory = q0 * cos(tspan) + p0 * sin(tspan);
% 短时间 ode45 1000步的bug都那么明显吗（默认精度下）？
csf1 = plot(tspan(1:1000),q_theory(1:1000),tspan(1:1000),ode45ns.q(1:1000),...
    tspan(1:1000),se.q(1:1000),tspan(1:1000),se2.q(1:1000),tspan(1:1000),se4.q(1:1000));
csf1(1).Color = 'k';
csf1(2).Color = 'g';
csf1(3).Color = 'y';
csf1(4).Color = 'r';
csf1(5).Color = 'b';
legend('Theory','ode45','Symp-1','Symp-2','Symp-4')
title('q-t曲线（短时间）')
xlabel('t')
ylabel('q')
% 长时间（只画ode45和symp4BCH）
figure
subplot(2,1,1)
plot(tspan,ode45ns.q);
legend('ode45')
title('q-t曲线（长时间）')
xlabel('t')
ylabel('q')
subplot(2,1,2)
plot(tspan,se4.q)
legend('ode45')
title('q-t曲线（长时间）')
xlabel('t')
ylabel('q')






% Energy fluctuation
E0 = H(q0,p0);
[E1,static1] = Energy_test(H,se.q,se.p);
[E2,static2] = Energy_test(H,se2.q,se2.p);
[E4,static4] = Energy_test(H,se4.q,se4.p);
[E45,static45] = Energy_test(H,ode45ns.q,ode45ns.p);
figure
plot(E45-E0)
hold on
plot(E1-E0)
plot(E2-E0)
plot(E4-E0)
xlabel('t')
ylabel('energy error')
legend('ode45','Symp-1','Symp-2','Symp-4')
title('Harmonic Oscillator:Energy error(1)')
figure
plot(E1-E0,'g')
hold on
plot(E2-E0,'b')
plot(E4-E0,'k')
xlabel('t')
ylabel('energy error')
legend('Symp-1','Symp-2','Symp-4')
title('Harmonic Oscillator:Energy error(2)')
disp('Energy Error for each method:')
disp('ode45')
disp(static45);
disp('symp1_Euler')
disp(static1);
disp('symp2_Verlet')
disp(static2);
disp('symp4_Neri')
disp(static4);


