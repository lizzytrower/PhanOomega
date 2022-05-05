function [D50] = ooidsizecorrection(b50,k2)
%OOIDSIZECORRECTION Kellerhaus et al 1975 grain size correction algorithm
%   This function implements eqn. 4 from Kellerhaus et al. (1975) to apply
%   a correction to measured b50 (median minor axis in 2D) from thin
%   sections. This was written with ooids in mind, where B = C
%   (intermediate and minor axis dimensions are equal; i.e. radial symmetry
%   about major axis) is a reasonable assumption and therefore k2 ~ 1. It
%   may not be the best correction for other grain shapes.

%   Written by Lizzy Trower in December 2021 with Matlab R2021b

load("Kellerhaus_correction.mat",'Kellerhaus_correction');
Cb_ratio = interp1(Kellerhaus_correction(:,1),...
    Kellerhaus_correction(:,2),k2);

D50 = (1 - Cb_ratio/100).*b50./(2*k2).*(2.*(1 + k2.^2)).^0.5;


end