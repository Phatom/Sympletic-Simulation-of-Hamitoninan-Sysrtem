%���ñ䲽�����Ľ����㷨��ode��ʾ�������Ľ����㷨�Ļ����ϼ�����ÿһ������ǰ�жϲ����Ĳ��֣�
%���setq��setpΪ�������Ľ����setstepsΪÿһ���Ĳ������ڷ�����⾫��Ҫ��ϸߵ�����£����ǿ��Կ��������仯�����ԡ�
function [setq,setp,setsteps] = symp4_variable_step(dqdt,dpdt,steps,t_max,q0,p0)%stepsΪ��ʼ����,t_maxΪ��ֹ�����
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