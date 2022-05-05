function [k,n] = calciteinterp_L(temp)
%uses data from Lopez et al., 2009
caln = [1.55, 1.84 ,2.3, 2.55];
callogk = [0.2, 0.33, 0.4, 0.51];
calT = [5, 25, 40, 55];
k = 10^interp1(calT, callogk, temp);
n = interp1(calT, caln, temp);

end