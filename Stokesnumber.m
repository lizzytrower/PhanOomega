function [Stk] = Stokesnumber(D,Rouse,rho_s,rho_f,nu)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
R = (rho_s - rho_f)/rho_f;
g = 9.8; %[m/s^2]

CSF = 1;  %1 is for spheres, 0.8 is for natural
PS = 6;  %6 is for spheres, 3.5 is for natural
%values for spheres are usually the appropriate choices for ooids

tauc = 0.03; %Critical Shields number.  0.03 is good for sand.

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

z0 = 3*D/30; %[m] roughness coefficient set by grainsize
dz = (H-z0)/1000; %[m]
z = z0:dz:H; %[m]
Uf = sum((ustar/0.41)*log(z/z0))*dz/H; %[m/s] depth-averaged flow velocity

%compute bed load height and velocity
hb = D*1.44*(tstage-1)^0.5; %[m] height of the bed load layer 
Us = (R*g*D)^0.5*1.56*(tstage-1)^0.56; %[m/s] bed load velocity

if Us>Uf
Us=Uf;
end

%compute suspended load
if hb < H 
    hb(hb<D)=D; %sets minimum ht of bed load layer
    b = hb; %bottom of the suspended load layer - same as height of bedload
    betta = 2; %Based on Scheingross et al. (2014) best fit
    P = ws/(0.41*ustar*betta); %Rouse number
    res = 1000;
    di5 = (log(H)-log(b))/res;
    i5 = log(b):di5:log(H);
    z = exp(i5);
    z(length(z))=H;
    dz = diff(z);
    dz = [dz(1) dz];
    a1=sum((((1-(z(z>z0)./H))/(1-(b/H))).*(b./z(z>z0))).^P.*log(z(z>z0)./z0).*dz(z>z0))...
        /(Uf*H)*(ustar/0.41);
    cb = 1/(Us*hb +Uf*H*a1);

    %find concentration profile
    c=0;
    c(1) = cb;
    c(2:length(z)+1) = cb.*(((1-(z./H))./(1-(b./H))).*(b./z)).^P ;
    z=[0 z];
    c(z==H)=0;

%calculate the fall distance
    gradc(1:length(c))=0;
    gradc(2:length(c)) = -diff(c);    
    Hfall = (1./cb).*sum(z.*gradc);
else
    hb = H;
    cb = 1/(Us.*hb);
    Hfall = hb;   
    a1=0;
end

if cb == 0
    Hfall = 0;
end

sig = ustar;
dx = sig/100; %the number of bins to subdivide - can play with this for resolution
X = -6*sig:dx:6*sig; %spread distribution for six sigma = w'/ws
f = normpdf(X,0,sig); %centered at zero normal gausian distribution
X = X./ws;  % to normalize as w'/ws same as psi

Scos = 1;  %cosine of the angle of the bed for impacts.

wfall = Scos*((2*(2/3)*D*g/cdrag*R)*(1-exp(-cdrag*rho_f/rho_s*...
    (Hfall/Scos)/(2/3*D))))^0.5;

psifall = wfall./ws;
settlematrix = psifall + X;
settlematrix1=settlematrix;
settlematrix(settlematrix<0) = 0; %no negative impacts

%Stokes number correction
wi_st = settlematrix;
wi = mean(wi_st);

Stk = D.*wi.*ws.*rho_s./(9.*nu.*rho_f);

end