function [ E, static ] = Energy_test( Hfun,q,p )
%ENERGY_TEST calculate the energy of a certain point in phase space for a given Hamilton H.
%   Input:
%	Hfun a function handle, describes the Hamilton of certain system.
%	q and p are state of the system, a vector wanted.
%   Output:
%   E is the energy for different time which should be a constant over time theoretically.
%	static:
%		static.fluc_E is \frac{\max{E}-\min{E}}{\max{E}} describes the energy fluctuation
%		static.avg_E and static.std_E gives the average and the standered error of E 
E = feval(Hfun,q,p);
static.fluc_E = (max(E) - min(E))/max(E);
static.avg_E = mean(E);
static.std_E = std(E);
end

