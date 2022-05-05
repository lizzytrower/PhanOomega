function [Omega,preciprate] = Oomegasolver(D,Rouse,intermittency,...
    rateconstant,rxnorder,rho_s,rho_f,nu,young)
%OOMEGASOLVER Function calculates paleo-Omega values based on ooid size
%   This function is an application of experimental work developed by
%   Trower et al. (2017) and tested in modern environments (marine - Trower
%   et al., 2018; lacustrine - Trower et al., 2020) to reconstruct calcium
%   carbonate saturation state in ancient seawater or lake water based on
%   measurements of ooid size (D) and estimates of transport mode (Rouse)
%   and intermittency (intermittency), with precipitation rate kinetics 
%   described by a rate constant (rateconstant) and reaction order
%   (rxnorder).

%   Required units for each input variable:
%   D - m (median grain diameter, assumed to be intermediate axis dimension
%   Rouse - dimensionless (Rouse number, recommended to use 2.5, reflecting
%           transport near the threshold of suspension)
%   intermittency - dimensionless (this parameter is called "f" in
%                   publications and must be in the range (0, 1])
%   rateconstant - umol/m^2/hr
%   rxnorder - dimensionless

%   IMPORTANT: This code as currently written is not designed to be run
%   with vectors as inputs, so it is recommended to run it with a for-loop
%   if you want to calculate multiple Omega estimates.

%   This code was developed by Lizzy Trower at the University of Colorado
%   Boulder in Matlab 2018b, last updated November 2020.

%The following parameters can be changed if needed to better reflect the
%system of interest.
%rho_s = 2800; %[kg/m^3] density of sediment, this value is set for aragonite
%rho_f = 1025; %[kg/m^3] density of fluid, this value is set for 25C seawater
R = (rho_s - rho_f)/rho_f;
g = 9.8; %[m/s^2]
%nu = 9.37*10^-7; %kinematic viscosity of fluid, this value is set for 25C seawater

%young = 56*10^9; %[Pa] young's modulus
strength = 1*10^6; %[Pa] tensile strength
kv = 9*10^5; %[dimensionless] see Trower et al. (2017) for explanation
tauc = 0.03; %Critical Shields number.  0.03 is good for sand.
Stc = 9; %critical Stokes number, see Trower et al. (2017) for explanation

CSF = 1;  %1 is for spheres, 0.8 is for natural
PS = 6;  %6 is for spheres, 3.5 is for natural
%values for spheres are usually the appropriate choices for ooids

H = 2; %[m] water depth

Dstar = (R.*g.*D.^3)./(nu.^2);
X = log10(Dstar);
R1 = -3.76715+1.92944.*X - 0.09815.*(X.^2) - 0.00575.*(X.^3) +...
    0.00056.*(X.^4);
R2 = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(X-4.6)) + ...
    0.3.*(0.5-CSF).*((1-CSF).^2).*(X-4.6);
R3 = (0.65-((CSF./2.83).*tanh(X-4.6))).^(1+((3.5-PS)./2.5));
Wstar = R3.*10.^(R2+R1);
ws = (R.*g.*nu.*Wstar).^(1./3);
cdrag = (4/3).*(R.*g.*D)./(ws.^2);

beta = 1;

ustar = ws./(Rouse.*.41.*beta);

tau = ustar^2/(R*g*D); %[dimensionless]
tstage = tau/tauc; %[dimensionless]
A1 = 0.36; %[dimensionless]

%this runs a separate script to calculate abrasion rate
susp_abrasion_calculations_abrcalc

%transform into volumetric abrasion rate
abrasionrate = Rabrasion*4*pi().*(D/2).^2; %[m^3/hr]

%calculate precipitation rate for an ooid at D_eq, incorporating the role
%of surface roughness (i.e., assuming the surface area available for
%precipitation is greater than the geometric surface area
SSA = 23.4; %ratio of actual surface area to geometric surface area, this
            %estimate is based on qualitative comparisons of surfaces in
            %Walter and Morse (1984) to SEM imagery of ooid surfaces; see
            %Trower et al. (2017)
SA_ssa = pi().*D.^2.*SSA; %[m^2] %adjusted surface area
preciprate = intermittency.*abrasionrate.*rho_s.*1000./100.0869./SA_ssa;

Omega = (preciprate./rateconstant).^(1./rxnorder) + 1;

end

