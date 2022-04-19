function output = gasViz(rho)
    sg = 0.6;
    T = 150;
    MW = sg*29;
	K = (9.379 + 0.01607*MW)*(T+460)^1.5/(209.2+19.26*MW+T);
	X = 3.448 + 986.4/(T+460) + 0.01009*MW;
	Y = 2.447 - 0.2224*X;
	output = 1e-4*K*exp(X*(rho/62.4)^Y);
end