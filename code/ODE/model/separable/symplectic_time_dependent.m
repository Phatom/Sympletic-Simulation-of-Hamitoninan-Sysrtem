% Symplectic integrator for time depedent hamiltonian.
% Mathiu Equation as an example
clear all
close all
clc
n=1;
t_max = 50;
time_steps=500;
tspan = linspace(0,t_max,time_steps);
omega = 1;
h=1;
epsilon=0.2;
disp('calculating diravitaves')


% Dynamic Setup ?????
% ?????????????????windows?matlab2017a???????????
% windows????2018a????mac????2018a???
% ????????????????????????????
tic
syms p q t
H = @(p,q,t)(p.^2-omega^2*(1+h*cos(2*omega+epsilon)*t).*q.^2)/2;
diff_p = diff(H,p);
dqdt = matlabFunction(diff_p);
diff_q = diff(H,q);
dpdt = matlabFunction(diff_q);
toc



% disp('calculating diravitaves')
% tic
% syms  p q
% H =  @(p,q,t)(p.^2-omega^2*(1+h*cos(2*omega+epsilon)*t).*q.^2)/2;
% diff_p = diff(H,p);
% dqdt = str2func(['@(p)',vectorize(diff_p)]);
% diff_q = diff(H,q);
% dpdt = str2func(['@(q)',vectorize(diff_q)]);
% toc

q0=0;
p0=0.5;


disp('symp4_Neri')
tic
dim = length(q0);  
Nt = length(tspan); 
dts = diff(tspan);
q =  zeros(dim,Nt); 
q(:,1)=q0;
p = zeros(dim,Nt); 
p(:,1)=p0;

c = [1,1-2^(1/3),1-2^(1/3),1]./(2*(2-2^(1/3)));
d = [1,-2^(1/3),1,0]./(2-2^(1/3)); 

for nt = 2:Nt
   dt = dts(nt-1);
   q(:,nt) = q(:,nt-1) + dt * c(1) * feval(dpdt, p(:,nt-1),tspan(nt))';
   p(:,nt) = p(:,nt-1) - dt * d(1) * feval(dqdt, q(:,nt))'; 
   q(:,nt) = q(:,nt) + dt * c(2) * feval(dpdt, p(:,nt),tspan(nt))';
   p(:,nt) = p(:,nt) - dt * d(2) * feval(dqdt, q(:,nt))';
   q(:,nt) = q(:,nt) + dt * c(3) * feval(dpdt, p(:,nt),tspan(nt))';
   p(:,nt) = p(:,nt) - dt * d(3) * feval(dqdt, q(:,nt))';
   q(:,nt) = q(:,nt) + dt * c(4) * feval(dpdt, p(:,nt),tspan(nt))';
end

q = q.';
p = p.';
toc


figure
plot(tspan,q)
title('Parametric Resonance')
xlabel('t')
ylabel('q')
text(5.8,-3.5,'\epsilon=')
text(7,-3.5,num2str(epsilon))