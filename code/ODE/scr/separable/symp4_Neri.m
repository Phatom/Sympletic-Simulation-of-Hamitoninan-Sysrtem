function [q,p] = symp4_Neri(dqdt,dpdt,tspan,q0,p0,varargin)
%symp4_Neri A fourth order symplectic solver, implicit method.
%  [Q P] = symp4_Neri(DQDT,DPDT,TSPAN,Q0,P0) solves Hamilton's equations 
%  when the Hamiltonian is separebel, i.e. H(p,q)=T(p)+Q(q) Functions DQDT(P) 
%  and DPDT(Q) must return T'(p) and - V'(q) in the form of a N-dimensional 
%  column vectors. Vectors Q0 and P0 are the initial 
%  conditions at T0. Each row in the solution arrays Q and P corresponds to a 
%  time specified in TSPAN.  
%
%  [Q P] = symp4_Neri(dqdt,dpdt,tspan,q0,p0,varargin) passes the additional
%  parameters R1, R2, ... to functions DQDT(P,R1,R2,...) and DPDT(Q,R1,
%  R2,...).
%
%
%  See also symp1_Euler, symp2_Verlet, symp4_Neri.

%  Yaoqi Zhang, 10-July-2018

dim = length(q0);  
Nt = length(tspan); 
dts = diff(tspan);
q =  zeros(dim,Nt); 
q(:,1)=q0;
p = zeros(dim,Nt); 
p(:,1)=p0;

%coefficients from Yoshida,‎1990
c = [1,1-2^(1/3),1-2^(1/3),1]./(2*(2-2^(1/3)));
d = [1,-2^(1/3),1,0]./(2-2^(1/3)); 

for nt = 2:Nt
   dt = dts(nt-1);
   q(:,nt) = q(:,nt-1) + dt * c(1) * feval(dpdt, p(:,nt-1),varargin{:})';
   p(:,nt) = p(:,nt-1) - dt * d(1) * feval(dqdt, q(:,nt),varargin{:})'; 
   q(:,nt) = q(:,nt) + dt * c(2) * feval(dpdt, p(:,nt),varargin{:})';
   p(:,nt) = p(:,nt) - dt * d(2) * feval(dqdt, q(:,nt),varargin{:})';
   q(:,nt) = q(:,nt) + dt * c(3) * feval(dpdt, p(:,nt),varargin{:})';
   p(:,nt) = p(:,nt) - dt * d(3) * feval(dqdt, q(:,nt),varargin{:})';
   q(:,nt) = q(:,nt) + dt * c(4) * feval(dpdt, p(:,nt),varargin{:})';
end

q = q.';
p = p.';