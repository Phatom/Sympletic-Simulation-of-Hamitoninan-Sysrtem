%�жϲ�����α仯�����Ƚ�epsilon��delta�Ĵ�С��ϵ���жϲ���Ӧ�����ӻ��Ǽ�С
function [steps,q,p]=solve_steps(dqdt,dpdt,p0,q0,steps)
epsilon=5e-16;
[delta,~,~]=solve_delta(dqdt,dpdt,q0,p0,steps);
s = delta - epsilon;
if s > 0
        while delta > epsilon
        steps = steps/2;
        [delta,q,p]=solve_delta(dqdt,dpdt,q0,p0,steps);
        end
elseif s < 0
        while delta < epsilon
        steps = steps*2;
        [delta,~,~]=solve_delta(dqdt,dpdt,q0,p0,steps);
        end
        steps = steps/2;%����ѭ����delta>epsilon��������Ҫ������Ҫ��Сһ�β���ʹdelta<epsilon
        [~,q,p]=solve_delta(dqdt,dpdt,q0,p0,steps);
end
end