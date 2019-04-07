%判断步长如何变化，即比较epsilon与delta的大小关系来判断步长应该增加还是减小
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
        steps = steps/2;%由于循环后delta>epsilon，不符合要求，所以要减小一次步长使delta<epsilon
        [~,q,p]=solve_delta(dqdt,dpdt,q0,p0,steps);
end
end