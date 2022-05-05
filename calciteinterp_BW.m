function [k,n] = calciteinterp_BW(temp)
%uses data from Burton & Walter 1989
caln = [0.6, 1.9, 2.3];
calk = [14, 3.9, 3.7];
calT = [5, 25, 37];
k = interp1(calT, calk, temp);
n = interp1(calT, caln, temp);

if temp > 37
    k = 3.7;
    n = 2.3;
end

end