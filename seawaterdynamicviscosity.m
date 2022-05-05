function [mu_sw] = seawaterdynamicviscosity(T,S)
%SEAWATERDYNAMICVISCOSITY as a function of temperature and salinity
%   J. D. Isdale, C. M. Spence, and J. S. Tudhope, Physical properties of 
% sea water solutions: viscosity, Desalination, 10(4), 319 - 328, 1972.
%   Inputs:
%   T : temperature in degrees Celsius
%   S : practical salinity in g/kg
%   Valid for T = 10 to 180 degrees C and S = 0 to 150 g/kg

A = 1.474*10^-3 + 1.5*10^-5*T - 3.927*10^-8*T^2;
B = 1.073*10^-5 - 8.5*10^-8*T + 2.23*10^-10*T^2;
mu_w = exp(-10.7019 + 604.129/(139.18 + T));
mu_sw = mu_w*(1 + A*S + B*S^2);

end