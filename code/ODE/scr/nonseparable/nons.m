function [q,p,x,y] = nons(dqdt,dpdt,tspan,q0,p0,order,varargin)
%NONS 此处显示有关此函数的摘要
%
    dt = diff(tspan);
    dim = length(q0);
    Nt = length(tspan);
    omega = 2;
    q =  zeros(dim,Nt);
    p = zeros(dim,Nt);
    x = zeros(dim,Nt);
    y = zeros(dim,Nt);
    q(:,1) = q0';
    p(:,1) = p0';
    x(:,1) = q0';
    y(:,1) = p0';

    for it = 1:Nt-1
        [q(:,it+1),p(:,it+1),x(:,it+1),y(:,it+1)] = update(dt(it),order,q(:,it),p(:,it),x(:,it),y(:,it));
    end



    function [omega] = heuristic_estimate_for_omega()
        c = 10;
        omega = (c*max(dt))^(-order);
    end

    function [q,p,x,y] = A_update(delta,q,p,x,y)
        x = x + delta * feval(dqdt,y,q,varargin{:})';
        p = p - delta * feval(dpdt,y,q,varargin{:})';
    end

    function [q,p,x,y] = B_update(delta,q,p,x,y)
        q = q + delta * feval(dqdt,p,x,varargin{:})';
        y = y - delta * feval(dpdt,p,x,varargin{:})';
    end

    function [qn,pn,xn,yn] = C_update(delta,q,p,x,y)
        c = cos(2*omega*delta);
        s = sin(2*omega*delta);
        qn = 0.5*((q+x) + c*(q-x) + s*(p-y));
        pn = 0.5*((p+y) - s*(q-x) + c*(p-y));
        xn = 0.5*((q+x) - c*(q-x) - s*(p-y));
        yn = 0.5*((p+y) + s*(q-x) - c*(p-y));
    end

    function [q,p,x,y] = update(dt,order,q,p,x,y)
        if order == 2
            [q,p,x,y] = A_update(0.5*dt,q,p,x,y);
            [q,p,x,y] = B_update(0.5*dt,q,p,x,y);
            [q,p,x,y] = C_update(dt,q,p,x,y);
            [q,p,x,y] = B_update(0.5*dt,q,p,x,y);
            [q,p,x,y] = A_update(0.5*dt,q,p,x,y);
        else
            gamma = 1 / (2 - 2^(1/(order+1)));
            [q,p,x,y] = update(gamma*dt,order-2,q,p,x,y);
            [q,p,x,y] = update((1-2*gamma)*dt,order-2,q,p,x,y);
            [q,p,x,y] = update(gamma*dt,order-2,q,p,x,y);
        end
    end

end

