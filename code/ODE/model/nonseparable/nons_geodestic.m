%  Schwarzschild geodesics 
cc
n = 3;
dt = 0.2;
time_steps=5000;
tspan = 0:dt:dt*time_steps;



q0 = [0,5,0];
p0 = [0.982,0,-4.472];

dqdt = @(p,q)[ -p(1)./(1-2./q(2)), p(2)*(1-2./q(2)), -p(3)./q(2).^2];
dpdt = @(p,q)[ 0, -(p(1)./q(2)).^2./(1-2./q(2))^2-(p(2)./q(2)).^2+p(3).^2/q(2).^3 ,0];


disp('symp4')
tic
[se4.q,se4.p,se4.x,se4.y] = nons(dqdt,dpdt,tspan,q0,p0,4);
toc

% polar
figure
polar(se4.q(3,:),se4.q(2,:))
title('Precesion of Merury')

% 动画 (comet太快了，啥都看不清)
[mx,my] = pol2cart(se4.q(3,:),se4.q(2,:));
figure
axis([-10,10,-10,10]);
title('Precesion of Merury')
h = animatedline(mx(1),my(1));
 for ix = 2:length(mx)
    h.addpoints(mx(ix),my(ix));
    drawnow;
end
