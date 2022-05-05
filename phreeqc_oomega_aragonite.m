function [Alk_out,DIC_out,pH_out,Omega_ca_out,pCO2_out,c,vargout] = ...
    phreeqc_oomega_aragonite(T_in,Ca_in,Mg_in,K_in,SO4_in,Na_in,...
    Cl_in,pCO2_in,Omega_in)
% PHREEQC_OOMEGA_ARAGONITE Use PHREEQC to compute Alk, DIC, pH, and
% Omega_calcite from Omega_aragonite estimated via ooid size data

% This code was developed by Lizzy Trower with Matlab R2021b, based on a
% Matlab function developed by Ted Present (PHREEQCO2).

try

% Default database location:
default_db = ...
    'C:\Program Files\USGS\IPhreeqcCOM 3.7.1-15876\database\phreeqc.dat';
%depending on your installation, this file location or name may need to be
%updated!!!

% Convert numeric compositions to text to feed to COM server
s_temp = num2str(T_in);
s_Ca = num2str(Ca_in);
s_Mg = num2str(Mg_in);
s_K = num2str(K_in);
s_SO4 = num2str(SO4_in);
s_Na = num2str(Na_in);
s_Cl = num2str(Cl_in);
s_logPCO2 = num2str(log10(pCO2_in/1e6));

% Script and run PHREEQC
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
ipc.AccumulateLine(['pH',sprintf('\t'),'8',sprintf('\t'),'aragonite',sprintf('\t'),num2str(log10(Omega_in))]);
ipc.AccumulateLine(['C(4)',sprintf('\t'),'2',sprintf('\t'),'CO2(g) ', s_logPCO2  ]);

% define results ouput
ipc.AccumulateLine ('SELECTED_OUTPUT');
ipc.AccumulateLine ('-reset false'); % default to not reporting outputs unless specified
ipc.AccumulateLine ('-high_precision true');
ipc.AccumulateLine ('-molalities CO2 HCO3- CO3-2');
ipc.AccumulateLine ('-totals C(4)');
ipc.AccumulateLine ('si aragonite');
ipc.AccumulateLine ('si calcite');
ipc.AccumulateLine ('si dolomite');
ipc.AccumulateLine ('si CO2(g)');
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
        vargout = [];
end
Alk_out = cell2mat(c(3,2));
DIC_out = cell2mat(c(4,2));
pH_out = cell2mat(c(1,2));
Omega_ca_out = cell2mat(c(9,2));
pCO2_out = cell2mat(c(11,2));

catch
    disp('PHREEQC error')
    Alk_out = NaN;
    DIC_out = NaN;
    pH_out = NaN;
    Omega_ca_out = NaN;
end
end