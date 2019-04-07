function [q,p] = symp1_Euler(dqdt,dpdt,tspan,q0,p0,varargin)
%symp1_Euler a first order symplectic solver, implicit method.
%  [Q P] = symp1_Euler(dqdt,dpdt,tspan,q0,p0) solves Hamilton's equations 
%  when the Hamilton is separebel, i.e. H(p,q)=T(p)+Q(q)
%  Functions DQDT(P) and DPDT(Q) must return T'(p) and - V'(q) in the form 
%  of a N-dimensional column vectors. Vectors Q0 and P0
%  are the initial conditions at T0. Each row in the solution arrays Q and
%  P corresponds to a time specified in TSPAN.  
%
%  [Q P] = symp1_Euler(dqdt,dpdt,tspan,q0,p0,varargin) passes the additional
%  parameters R1, R2, ... to functions DQDT(P,R1,R2,...) and DPDT(Q,R1,
%  R2,...).
%
%
%  See also symp2_Verlet, symp4_Ruth, symp4_Neri.

%  Francisco J. Beron-Vera, 28-Mar-2005
%  $Revision: 1.0 $  $Date: 28-Mar-2005 14:58:31 $
%  Modified by Yaoqi Zhang, 10-July-2018

N = length(q0);  
Nt = length(tspan); 
hs = diff(tspan);
q = zeros(N,Nt); q(:,1)=q0;
p = zeros(N,Nt); p(:,1)=p0;
% updating
for nt = 2:Nt
   h = hs(nt-1);
   q(:,nt) = q(:,nt-1) + h * feval(dqdt, p(:,nt-1), varargin{:});
   p(:,nt) = p(:,nt-1) - h * feval(dpdt, q(:,nt)  , varargin{:});
end
q = q.';
p = p.';