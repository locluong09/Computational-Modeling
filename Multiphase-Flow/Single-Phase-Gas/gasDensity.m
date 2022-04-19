function output = gasDensity(p,z)
    sg = 0.6; 
    T = 150;
    MW = sg*29;
	output = MW*p/z/10.7316/(T+460);
end