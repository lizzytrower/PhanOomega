function [Alk_out,DIC_out,pH_out,Omega_ar_out,c,vargout] = ...
    phreeqc_oomega_calcite(T_in,Ca_in,Mg_in,K_in,SO4_in,Na_in,...
Cl_in,pCO2_in,Omega_in)
% PHREEQCO2 Use PHREEQC to compute equilibrium carbonate system from any
% pair of DIC, ALK, pCO2, or pH
%
%   Ted Present, 2020
%   Requires installation of <a href="matlab:
%   web('https://www.usgs.gov/software/phreeqc-version-3')">IPhreeqcCOM</a>.
% 
%   
%
%   For more information, see the <a href="matlab:
%   web('https://www.usgs.gov/software/phreeqc-version-3')">USGS PHREEQC website</a>.
%
%   See also CO2SYS.
try
%% Defaults
% Default database location:
default_db = 'C:\Program Files\USGS\IPhreeqcCOM 3.7.1-15876\database\phreeqc.dat';

%% Convert numeric compositions to text to feed to COM server
s_temp = num2str(T_in);
s_Ca = num2str(Ca_in);
s_Mg = num2str(Mg_in);
s_K = num2str(K_in);
s_SO4 = num2str(SO4_in);
s_Na = num2str(Na_in);
s_Cl = num2str(Cl_in);
s_logPCO2 = num2str(log10(pCO2_in/1e6));

%% Script and run PHREEQC
% initialize IPhreeqcCOM server and database
ipc = actxserver('IPhreeqcCOM.Object');
ipc.LoadDatabase(default_db);
ipc.ClearAccumulatedLines;

% define solution
ipc.AccumulateLine ('SOLUTION 1 ');
ipc.AccumulateLine ('units mmol/kgw');
ipc.AccumulateLine([   'temp',     sprintf('\t'),      s_temp   ]);
ipc.AccumulateLine([   'Ca',       sprintf('\t'),      s_Ca     ]);
ipc.AccumulateLine([   'Mg',       sprintf('\t'),      s_Mg     ]);
ipc.AccumulateLine([   'K',        sprintf('\t'),      s_K      ]);
ipc.AccumulateLine([   'S(6)',     sprintf('\t'),      s_SO4,     sprintf('\t'), 'as SO4'     ]);
ipc.AccumulateLine([   'Na',       sprintf('\t'),      s_Na     ]);
ipc.AccumulateLine([   'Cl',       sprintf('\t'),      s_Cl     ]);
ipc.AccumulateLine(['pH',sprintf('\t'),'8',sprintf('\t'),'calcite',sprintf('\t'),num2str(log10(Omega_in))]);
ipc.AccumulateLine(['C(4)',sprintf('\t'),'2',sprintf('\t'),'CO2(g) ', s_logPCO2  ]);

% define results ouput
ipc.AccumulateLine ('SELECTED_OUTPUT');
ipc.AccumulateLine ('-reset false'); % default to not reporting outputs unless specified
ipc.AccumulateLine ('-high_precision true');
ipc.AccumulateLine ('-activities CO2 HCO3- CO3-2');
ipc.AccumulateLine ('-totals C(4)');
ipc.AccumulateLine ('si aragonite');
ipc.AccumulateLine ('si calcite');
ipc.AccumulateLine ('si dolomite');
ipc.AccumulateLine ('si CO2(g)');
%ipc.AccumulateLine ('si gypsum');
%ipc.AccumulateLine ('si halite');
ipc.AccumulateLine ('-ph  true');
ipc.AccumulateLine ('-temperature  true');
ipc.AccumulateLine ('-alkalinity  true');

% Run PHREEQC and get results
ipc.RunAccumulated;
output = ipc.GetSelectedOutputArray';

%% Return results
switch nargout
    case {0,1}
        c = output;
    otherwise
        c = output;
        vargout = []; % find a way to order output manually
%         idx = find(strcmp(output,'si_halite'));
%         vargout(idx) = output{idx,2};
end
Alk_out = cell2mat(c(3,2));
DIC_out = cell2mat(c(4,2));
pH_out = cell2mat(c(1,2));
Omega_ar_out = cell2mat(c(8,2));
catch
    disp('PHREEQC error')
    Alk_out = NaN;
    DIC_out = NaN;
    pH_out = NaN;
    Omega_ar_out = NaN;
end
end