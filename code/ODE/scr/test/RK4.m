function y = RK4(F_xy,h,x,y0)
% an RK4 method to solve ODE
	y = zeros(1,length(x)); 
	y(1) = y0;                                          % initial condition

	for i=1:(length(x)-1)                              % calculation loop
		k_1 = F_xy(x(i),y(i));
		k_2 = F_xy(x(i)+0.5*h,y(i)+0.5*h*k_1);
		k_3 = F_xy((x(i)+0.5*h),(y(i)+0.5*h*k_2));
		k_4 = F_xy((x(i)+h),(y(i)+k_3*h));
		
		y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
	end
end