% ¿Ì¬€ phi_4 phi_2 ode45 
cc
n = 1;
dt = 0.1;
time_steps=2000; % When time_steps is bigger than 10000 , ode45 failed (using the default RelTol and AbsTol ).
tspan = 0:dt:dt*time_steps;

disp('calculating diravitaves')
tic
syms p q
H = @(p,q)(p.^2+1).*(q.^2+1)/2;
diff_p = diff(H,p);
dqdt = matlabFunction(diff_p);
diff_q = diff(H,q);
dpdt = matlabFunction(diff_q);
disp(dqdt)
disp(dpdt)
toc

q0 = -3;
p0 = 0;

disp('theory')
tic
[~,cn,~] = ellipj(sqrt(1+q0^2)*tspan,q0^2/(1+q0^2));
q_theory = q0 * cn;
toc

disp('symp2')
tic
[se2.q,se2.p,se2.x,se2.y] = nons(dqdt,dpdt,tspan,q0,p0,2);
toc

disp('symp4')
tic
[se4.q,se4.p,se4.x,se4.y] = nons(dqdt,dpdt,tspan,q0,p0,4);
toc

disp('ode45')
odefun=@(t,y)([y(2)*(y(1)^2+1);-y(1)*(y(2)^2+1)]);
y0=[q0,p0];
% Solve
tic
[~,y]=ode45(odefun,tspan,y0); 
toc 
ode45ns.q=y(:,1);
ode45ns.p=y(:,2);



% Test and Plot
% Phase Space
figure(1)
subplot(2,2,1)
plot(tspan,q_theory)
xlabel('q')
ylabel('t')
legend('Theory')
title('q-t')
subplot(2,2,2)
plot(tspan,ode45ns.q)
xlabel('q')
ylabel('t')
legend('ode45')
title('q-t')
subplot(2,2,3)
plot(tspan,se2.q)
xlabel('q')
ylabel('t')
legend('symp2')
title('q-t')
subplot(2,2,4)
plot(tspan,se4.q)
xlabel('q')
ylabel('t')
legend('symp4')
title('q-t')

% q-t
figure(2)
subplot(3,1,1)
plot(ode45ns.q,ode45ns.p)
title('Phase Space')
xlabel('q')
ylabel('p')
legend('ode45')
subplot(3,1,2)
plot(se2.q,se2.p)
xlabel('q')
ylabel('p')
legend('symp2')
subplot(3,1,3)
plot(se4.q,se4.p)
xlabel('q')
ylabel('p')
legend('symp4')

% Energy
[E45, static45] = Energy_test( H,ode45ns.q,ode45ns.p );
[E2,static2] = Energy_test(H,se2.q,se2.p);
[E4,static4] = Energy_test(H,se4.q,se4.p);
E0 = feval(H,p0,q0);
figure(3)
subplot(3,1,1)
plot(E45-E0)
legend('ode45')
title('Energy Error')
subplot(3,1,2)
plot(E2-E0)
legend('symp2')
subplot(3,1,3)
plot(E4-E0)
legend('symp4')
disp('ode45')
disp(static45)
disp('symp2')
disp(static2)
disp('symp4')
disp(static4)