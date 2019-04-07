function [q,p] = symp4_Ruth(dqdt,dpdt,tspan,q0,p0,varargin)
%symp4_Ruth A fourth order symplectic solver, implicit method.
%  [Q P] = symp4_Ruth(DQDT,DPDT,TSPAN,Q0,P0) solves Hamilton's equations 
%  when the Hamiltonian is separebel, i.e. H(p,q)=T(p)+Q(q) Functions DQDT(P) 
%  and DPDT(Q) must return T'(p) and - V'(q) in the form of a N-dimensional 
%  column vectors. Vectors Q0 and P0 are the initial 
%  conditions at T0. Each row in the solution arrays Q and P corresponds to a 
%  time specified in TSPAN.  
%
%  [Q P] = symp4_Ruth(dqdt,dpdt,tspan,q0,p0,varargin) passes the additional
%  parameters R1, R2, ... to functions DQDT(P,R1,R2,...) and DPDT(Q,R1,
%  R2,...).
%
%
%  See also symp1_Euler, symp2_Verlet, symp4_Neri.

%  Yaoqi Zhang, 10-July-2018

N = length(q0);  
Nt = length(tspan); 
hs = diff(tspan);
q =  zeros(N,Nt); q(:,1)=q0;
p = zeros(N,Nt); p(:,1)=p0;
%coefficients from E Forest&RD Ruth,1989 
c = [0,1/(2-2^(1/3)),1/(1-2^(1/3)),1/(2-2^(1/3))];
d = [2+2^(1/3)+2^(-1/3),1-2^(1/3)-2^(-1/3),1-2^(1/3)-2^(-1/3),2+2^(1/3)+2^(-1/3)]./6;
% updating
for nt = 2:Nt
   h = hs(nt-1);
   p(:,nt) = p(:,nt-1) - h * d(1) * feval(dqdt, q(:,nt-1),varargin{:});
   p(:,nt) = p(:,nt) - h * d(2) * feval(dqdt, q(:,nt-1),varargin{:});
   q(:,nt) = q(:,nt-1) + h * c(2) * feval(dpdt, p(:,nt),varargin{:});
   p(:,nt) = p(:,nt) - h * d(3) * feval(dqdt, q(:,nt),varargin{:});
   q(:,nt) = q(:,nt) + h * c(3) * feval(dpdt, p(:,nt),varargin{:});
   p(:,nt) = p(:,nt) - h * d(4) * feval(dqdt, q(:,nt),varargin{:});
   q(:,nt) = q(:,nt) + h * c(4) * feval(dpdt, p(:,nt),varargin{:});
end
q = q.';
p = p.';