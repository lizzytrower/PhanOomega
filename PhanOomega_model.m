%%
%PhanOomega: reconstructing the carbonate chemistry of Phanerozoic seawater
%using ooid size data

%This code was written by Lizzy Trower in Matlab R2021b. This version was
%last edited in May 2022.

%%
%Load ooid size data

clear

%This data has already been pulled from a spreadsheet and organized into
%this .mat file
load('AKdata_formation.mat')
ages = plot_age;

reps = 10000;
ooidD_dist = zeros(length(ages),reps);
for Dcount = 1:length(ages)
    zz = rand(1,reps);
    Dpdf = makedist('Normal','mu',D50_corr(Dcount),'sigma',D_range(Dcount));
    Dpdf_trunc = truncate(Dpdf,100,inf);
    ooidD_dist(Dcount,:) = icdf(Dpdf_trunc,zz);
end

%%
%Temperature data

%Load Scotese et al 2021 tropical temperature data

troptemps = readmatrix('Part 4. Phanerozoic_Paleotemperature_Summaryv4.xlsx',...
    'Sheet','TropvsGAT 1my','Range','A2:C542');
troptemps(:,2) = [];
troptemps_rescale = interp1(troptemps(:,1),troptemps(:,2),ages);

troptemps_dist = zeros(length(ages),reps);
for tempcount = 1:length(ages)
    z0 = rand(1,reps);
    temppdf = makedist('Normal','mu',troptemps_rescale(tempcount),...
        'sigma',0.05*troptemps_rescale(tempcount));
    troptemps_dist(tempcount,:) = icdf(temppdf,z0);
end

%Optional code to plot Scotese et al tropical temperature data

% figure
% scatter(ages,prctile(troptemps_dist,50,2),'filled','m')
% hold on
% scatter(ages,prctile(troptemps_dist,10,2),'m')
% scatter(ages,prctile(troptemps_dist,90,2),'m')
% ylabel(['tropical T (' char(176) 'C)'])
% box on
% xlim([0 550])
% set(gca, 'xdir', 'reverse')
% xlabel('age (Ma)')
% ylim([20 40])

%%
%Load Phanerozoic seawater salinity curve from Hay et al 2006

salinity = readmatrix('Hay_etal_salinity.xlsx','Range','A2:B57');
sal_rescale = interp1(salinity(:,1),salinity(:,2),ages);
sal_rescale(1:4) = sal_rescale(5);
sal_rescale(end-7:end) = sal_rescale(end-8);

sal_dist = zeros(length(ages),reps);
for salcount = 1:length(ages)
    z0 = rand(1,reps);
    salpdf = makedist('Normal','mu',sal_rescale(salcount),...
        'sigma',0.05*sal_rescale(salcount));
    sal_dist(salcount,:) = icdf(salpdf,z0);
end


%%
%Load Demicco et al 2005 Ca and Mg data

[concCa_D] = readmatrix('Demicco_panelD_data.xlsx','Sheet','Ca',...
    'Range','A2:B45');
[concMg_D] = readmatrix('Demicco_panelD_data.xlsx','Sheet','Mg',...
    'Range','A2:B44');

%%
%Load other ion data from Demicco et al 2005

[concSO4] = readmatrix('Demicco_panelD_data.xlsx','Sheet','SO4',...
    'Range','A2:B43');
[concNa] = readmatrix('Demicco_panelD_data.xlsx','Sheet','Na',...
    'Range','A2:B45');
[concK] = readmatrix('Demicco_panelD_data.xlsx','Sheet','K',...
    'Range','A2:B45');

%%
%pCO2 data

%Load pCO2 data from Foster et al. 2017
pCO2_420to0 = readmatrix('41467_2017_BFncomms14845_MOESM2875_ESM.xlsx',...
    'Range','A3:E842');
pCO2_420to0_sigma = pCO2_420to0(:,5) - pCO2_420to0(:,2);

%Load model pCO2 from Berner 2006
pCO2_540to420 = readmatrix('GEOCARBSULF.xlsx','Range','A2:C19');
pCO2_540to420(:,2) = [];
pCO2_540to420(15:end,:) = [];

%combine
pCO2_stack = cat(1,pCO2_420to0(:,1:2),pCO2_540to420);

pCO2_rescale = interp1(pCO2_stack(:,1),pCO2_stack(:,2),ages);
pCO2_sigma_rescale = interp1(pCO2_420to0(:,1),pCO2_420to0_sigma,ages);

pCO2_dist = zeros(length(ages),reps);
for pCO2count = 1:find(ages>420,1,'last')
    zpCO2 = rand(1,reps);
    pCO2pdf = makedist('Normal','mu',pCO2_rescale(pCO2count),...
        'sigma',0.2*pCO2_rescale(pCO2count));
    pCO2_dist(pCO2count,:) = icdf(pCO2pdf,zpCO2);
end
for pCO2count = 1+find(ages>420,1,'last'):length(ages)
    zpCO2 = rand(1,reps);
    pCO2pdf = makedist('Normal','mu',pCO2_rescale(pCO2count),...
        'sigma',pCO2_sigma_rescale(pCO2count));
    pCO2_dist(pCO2count,:) = icdf(pCO2pdf,zpCO2);
end

%%
%Load fluid inclusion [Ca] data and linearly interpolate

%As written, this data is not used, but could be used in place of Demicco
%et al. [Ca] and [Mg] data

% %Lowenstein et al 2003 data
% [Ca_L] = readmatrix('Lowenstein_etal_2003.xlsx','Range','A2:B7');
% 
% %Horita et al 2003 data
% [Ca_H] = readmatrix('Horita_etal_2003.xlsx','Sheet','Mg, Ca, SO4',...
%     'Range','A2:D13');
% Ca_H = Ca_H(:,[1,3]);
% 
% Ca_alldata = cat(1,Ca_H,Ca_L);
% Ca_alldata = sortrows(Ca_alldata);
% Ca_interp = interp1(Ca_alldata(:,1),Ca_alldata(:,2),ages);
% Ca_interp(end-7:end) = 10.3;
% 
% concCa_dist = zeros(length(ages),reps);
% for Cacount = 1:length(ages)
%     zCa = rand(1,reps);
%     Capdf = makedist('Normal','mu',Ca_interp(Cacount),...
%         'sigma',5);
%     Capdf_trunc = truncate(Capdf,0.1,inf);
%     concCa_dist(Cacount,:) = icdf(Capdf_trunc,zCa);
% end

%%
%Rescale datasets so they're in the same age bins as the ooid size data

concCa_rescale = interp1(concCa_D(:,1),concCa_D(:,2),ages);
concCa_rescale(end-7:end) = 10.3;

concCa_dist = zeros(length(ages),reps);
for Cacount = 1:length(ages)
    zCa = rand(1,reps);
    Capdf = makedist('Normal','mu',concCa_rescale(Cacount),...
        'sigma',5);
    Capdf_trunc = truncate(Capdf,0.1,inf);
    concCa_dist(Cacount,:) = icdf(Capdf_trunc,zCa);
end

concMg_rescale = interp1(concMg_D(:,1),concMg_D(:,2),ages);
concMg_rescale(end-7:end) = 53;

concNa_rescale = interp1(concNa(:,1),concNa(:,2),ages);
concNa_rescale(end-7:end) = 469;

concK_rescale = interp1(concK(:,1),concK(:,2),ages);
concK_rescale(end-7:end) = 10.2;

concSO4_rescale = interp1(concSO4(:,1),concSO4(:,2),ages);
concSO4_rescale(end-9:end) = 28;

%Load and rescale [Cl] data from Hay et al. 2006
concCl = readmatrix('Hay_etal_chloride.xlsx','Range','A2:E32');
concCl_rescale = interp1(concCl(:,1),concCl(:,5),ages);
concCl_rescale(1:4) = concCl_rescale(5);

%%
%Run Omegasolver to estimate paleo-Omega from ooid size data

%This uses parfor to run for loops in parallel. "parfor" can be replaced by
%"for" if running in parallel isn't possible, but recommended to keep
%parallelized if possible - this code takes a while to run!

fdist = zeros(length(ages),reps);
for fcount = 1:length(ages)
    fpdf = makedist('Normal','mu',0.15,'sigma',0.1);
    fpdf_trunc = truncate(fpdf,0.01,0.3);
    z0 = rand(1,reps);
    fdist(fcount,:) = icdf(fpdf_trunc,z0);
end

Oomega = zeros(length(ages),reps);
Rouse = 2.5;
rho_calc = 2700; %(kg/m^3)
rho_arag = 2800; %(kg/m^3)
Stokes = zeros(length(ages),reps);

parfor agecount = 1:length(ages)
    try
        for MCcount = 1:reps
            try
                temp = troptemps_dist(agecount,MCcount);
                sal = sal_dist(agecount,MCcount);
                f = fdist(agecount,MCcount);
                D = ooidD_dist(agecount,MCcount);
                rho_sw = seawaterdensity(temp,sal);
                mu_sw = seawaterdynamicviscosity(temp,sal);
                nu_sw = mu_sw/rho_sw;
        
                    if mineral(agecount) == {'calcite'}
                        [k,n] = calciteinterp_BW(temp);
                        young = 56*10^9;
                        Oomega(agecount,MCcount) = Oomegasolver(D*10^-6,Rouse,...
                            f,k,n,rho_calc,rho_sw,nu_sw,young);
                        Stokes(agecount,MCcount) = Stokesnumber(D*10^-6,Rouse,...
                            rho_calc,rho_sw,nu_sw);
                    elseif mineral(agecount) == {'aragonite'}
                        [k,n] = aragoniteinterp(temp);
                        young = 87*10^9;
                        Oomega(agecount,MCcount) = Oomegasolver(D*10^-6,Rouse,...
                            f,k,n,rho_arag,rho_sw,nu_sw,young);
                        Stokes(agecount,MCcount) = Stokesnumber(D*10^-6,Rouse,...
                            rho_arag,rho_sw,nu_sw);
                    elseif mineral(agecount) == {'bimineralic'}
                        [k,n] = aragoniteinterp(temp);
                        young = 87*10^9;
                        Oomega(agecount,MCcount) = Oomegasolver(D*10^-6,Rouse,...
                            f,k,n,rho_arag,rho_sw,nu_sw,young);
                        Stokes(agecount,MCcount) = Stokesnumber(D*10^-6,Rouse,...
                            rho_arag,rho_sw,nu_sw);

%swap between previous codeblock and following one to change whether
%bimineralic ooids are treated as aragonite or calcite

%                     elseif mineral(agecount) == {'bimineralic'}
%                         [k,n] = calciteinterp_BW(temp);
%                         young = 56*10^9;
%                         Oomega(agecount,MCcount) = Oomegasolver(D*10^-6,Rouse,...
%                             f,k,n,rho_calc,rho_sw,nu_sw,young);
%                         Stokes(agecount,MCcount) = Stokesnumber(D*10^-6,Rouse,...
%                             rho_calc,rho_sw,nu_sw);
                    end
            catch
                disp(MCcount)
            end
        end
    catch
        disp(agecount)
    end
end

mustBeNonNan(Oomega)

%%
%Run PHREEQC to calculate [DIC], Alk, and pH for the reconstructed Omega
%values

Alk = zeros(length(ages),reps);
DIC = zeros(length(ages),reps);
pH = zeros(length(ages),reps);
Omega_ar = zeros(length(ages),reps);
Omega_ca = zeros(length(ages),reps);

%This section also uses "parfor" so that the loop can run faster. Can be
%replaced with "for" if needed, but this section is also slow!

%The "phreeqc_oomega" functions may need to be edited with the correct file
%location for wherever the phreeqc database is installed on your computer.

parfor count3 = 1:length(ages)
    for count4 = 1:reps
        try
            temp = troptemps_dist(count3,count4); 
            if mineral(count3) == {'calcite'}
                [Alk(count3,count4),DIC(count3,count4),pH(count3,count4),...
                    Omega_ar(count3,count4)] = ...
                    phreeqc_oomega_calcite(temp,...
                    concCa_dist(count3,count4),concMg_rescale(count3),...
                    concK_rescale(count3),concSO4_rescale(count3),...
                    concNa_rescale(count3),concCl_rescale(count3),...
                    pCO2_dist(count3,count4),Oomega(count3,count4));
                Omega_ar(count3,count4) = 10^Omega_ar(count3,count4);
                Omega_ca(count3,count4) = Oomega(count3,count4);
            elseif mineral(count3) == {'aragonite'}
                [Alk(count3,count4),DIC(count3,count4),pH(count3,count4),...
                    Omega_ca(count3,count4)] = ...
                    phreeqc_oomega_aragonite(temp,...
                    concCa_dist(count3,count4),concMg_rescale(count3),...
                    concK_rescale(count3),concSO4_rescale(count3),...
                    concNa_rescale(count3),concCl_rescale(count3),...
                    pCO2_dist(count3,count4),Oomega(count3,count4));
                Omega_ca(count3,count4) = 10^Omega_ca(count3,count4);
                Omega_ar(count3,count4) = Oomega(count3,count4);
            elseif mineral(count3) == {'bimineralic'}
                [Alk(count3,count4),DIC(count3,count4),pH(count3,count4),...
                    Omega_ca(count3,count4)] = ...
                    phreeqc_oomega_aragonite(temp,...
                    concCa_dist(count3,count4),concMg_rescale(count3),...
                    concK_rescale(count3),concSO4_rescale(count3),...
                    concNa_rescale(count3),concCl_rescale(count3),...
                    pCO2_dist(count3,count4),Oomega(count3,count4));
                Omega_ca(count3,count4) = 10^Omega_ca(count3,count4);
                Omega_ar(count3,count4) = Oomega(count3,count4);

%swap between previous codeblock and following one to change whether
%bimineralic ooids are treated as aragonite or calcite

%             elseif mineral(count3) == {'bimineralic'}
%                 [Alk(count3,count4),DIC(count3,count4),pH(count3,count4),...
%                     Omega_ar(count3,count4)] = ...
%                     phreeqc_oomega_calcite(temp,...
%                     concCa_dist(count3,count4),concMg_rescale(count3),...
%                     concK_rescale(count3),concSO4_rescale(count3),...
%                     concNa_rescale(count3),concCl_rescale(count3),...
%                     pCO2_dist(count3,count4),Oomega(count3,count4));
%                 Omega_ar(count3,count4) = 10^Omega_ar(count3,count4);
%                 Omega_ca(count3,count4) = Oomega(count3,count4);
            end
        catch
            disp([count3, count4])
        end
    end
end

%adjust Alk and DIC units
Alk = Alk*10^3;
DIC = DIC*10^3;

%%
%Export data!

%Recommend editing the filename to include the date and other relevant
%details about model configuration.

save('PhanOomega_data_DATE.mat','ages','DIC','Alk','pH','Omega_ar',...
    'Omega_ca','Oomega','troptemps_dist','mineral','depenv','ooidD_dist',...
    'Stokes')
