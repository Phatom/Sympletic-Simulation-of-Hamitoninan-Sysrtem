function [q,p] = symp2_Verlet( dqdt,dpdt,tspan,q0,p0,varargin )
%symp2_Verlet A second order symplectic solver, implicit method.
%  [Q P] = symp2_Verlet(dqdt,dpdt,tspan,q0,p0) solves Hamilton's equations 
%  when the Hamiltonian is separebel, i.e. H(p,q)=T(p)+Q(q) Functions DQDT(P) 
%  and DPDT(Q) must return T'(p) and - V'(q) in the form 
%  of a N-dimensional column vectors. Vectors Q0 and P0 are the initial 
%  conditions at T0. Each row in the solution arrays Q and P corresponds to a 
%  time specified in TSPAN.  
%
%  [Q P] = symp2_Verlet(dqdt,dpdt,tspan,q0,p0,varargin) passes the additional
%  parameters R1, R2, ... to functions DQDT(P,R1,R2,...) and DPDT(Q,R1,
%  R2,...).
%
%
%  See also symp1_Euler, symp4_Ruth, symp4_Neri.

%  Yaoqi Zhang, 10-July-2018
N = length(q0);
Nt = length(tspan);
hs = diff(tspan);
q = zeros(N,Nt); q(:,1)=q0;
p = zeros(N,Nt); p(:,1)=p0;
% updating
for nt = 2:Nt
    h = hs(nt-1);
    q(:,nt) = next_q(q(:,nt-1),p(:,nt-1));
    p(:,nt) = next_p(q(:,nt),p(:,nt-1));
    q(:,nt) = next_q(q(:,nt),p(:,nt));
end
q = q.';
p = p.';
	%the coefficient from Yoshida,?1990
    function [q_next] = next_q(q_now,p_now)
        q_next = q_now + 0.5 * dpdt(p_now) * h;
    end

    function [p_next] = next_p(q_next,p_now)
        p_next = p_now - dqdt(q_next) * h;
    end

end

