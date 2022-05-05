%%
%PhanOomega: reconstructing the carbonate chemistry of Phanerozoic seawater
%using ooid size data

%This code loads and organizes data from a spreadsheet and generates a
%running plot to help make sure data is being read correctly.

%This code was written by Lizzy Trower in Matlab R2021b. This version was
%last edited in May 2022.

%%
clear

%%
%load in and organize Cambrian data

%load in the numerical data
num_Cambrian = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Cambrian','Range','A2:H998');
num_Cambrian(:,[3,7]) = []; %delete text column
%columns:
    %1: formation group
    %2: paleolatitude
    %3: grain number
    %4: major axis (um)
    %5: minor axis (um) **USE THIS ONE**
    %6: plot age (Ma) - midpoint of the given age bin
%end
D_Cambrian = num_Cambrian(:,5);
age_Cambrian = num_Cambrian(:,6);
formgp_Cambrian = num_Cambrian(:,1);
paleolat_Cambrian = num_Cambrian(:,2);

%load in the categorical data:
cat_Cambrian = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Cambrian','Range','K2:N998');
cat_Cambrian(:,2:3) = []; %delete the non-categorical columns
dep_Cambrian = categorical(cat_Cambrian(:,1));
min_Cambrian = categorical(cat_Cambrian(:,2));

%%
%plot data to make sure it looks ok
f_ooidsize = figure;
boxplot(D_Cambrian,formgp_Cambrian,'Positions',age_Cambrian,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Cambrian)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Ordovician data

%load in the numerical data
num_Ord = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Ordovician','Range','A2:H781');
num_Ord(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Ord = num_Ord(:,5);
age_Ord = num_Ord(:,6);
formgp_Ord = num_Ord(:,1) + max(formgp_Cambrian);
paleolat_Ord = num_Ord(:,2);

%load in the categorical data:
cat_Ord = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Ordovician','Range','K2:N781');
cat_Ord(:,2:3) = []; %delete the non-categorical columns
dep_Ord = categorical(cat_Ord(:,1));
min_Ord = categorical(cat_Ord(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Ord,formgp_Ord,'Positions',age_Ord,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Ord)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Silurian data

%load in the numerical data
num_Sil = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Silurian','Range','A2:H282');
num_Sil(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Sil = num_Sil(:,5);
age_Sil = num_Sil(:,6);
formgp_Sil = num_Sil(:,1) + max(formgp_Ord);
paleolat_Sil = num_Sil(:,2);

%load in the categorical data:
cat_Sil = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Silurian','Range','K2:N282');
cat_Sil(:,2:3) = []; %delete the non-categorical columns
dep_Sil = categorical(cat_Sil(:,1));
min_Sil = categorical(cat_Sil(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Sil,formgp_Sil,'Positions',age_Sil,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Sil)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Devonian data

%load in the numerical data
num_Dev = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Devonian','Range','A2:H601');
num_Dev(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Dev = num_Dev(:,5);
age_Dev = num_Dev(:,6);
formgp_Dev = num_Dev(:,1) + max(formgp_Sil);
paleolat_Dev = num_Dev(:,2);

%load in the categorical data:
cat_Dev = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Devonian','Range','K2:N601');
cat_Dev(:,2:3) = []; %delete the non-categorical columns
dep_Dev = categorical(cat_Dev(:,1));
min_Dev = categorical(cat_Dev(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Dev,formgp_Dev,'Positions',age_Dev,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Dev)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Carboniferous data

%load in the numerical data
num_Carb = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Carboniferous','Range','A2:H514');
num_Carb(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Carb = num_Carb(:,5);
age_Carb = num_Carb(:,6);
formgp_Carb = num_Carb(:,1) + max(formgp_Dev);
paleolat_Carb = num_Carb(:,2);

%load in the categorical data:
cat_Carb = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Carboniferous','Range','K2:N514');
cat_Carb(:,2:3) = []; %delete the non-categorical columns
dep_Carb = categorical(cat_Carb(:,1));
min_Carb = categorical(cat_Carb(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Carb,formgp_Carb,'Positions',age_Carb,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Carb)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Permian data

%load in the numerical data
num_Perm = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Permian','Range','A2:H1012');
num_Perm(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Perm = num_Perm(:,5);
age_Perm = num_Perm(:,6);
formgp_Perm = num_Perm(:,1) + max(formgp_Carb);
paleolat_Perm = num_Perm(:,2);

%load in the categorical data:
cat_Perm = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Permian','Range','K2:N1012');
cat_Perm(:,2:3) = []; %delete the non-categorical columns
dep_Perm = categorical(cat_Perm(:,1));
min_Perm = categorical(cat_Perm(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Perm,formgp_Perm,'Positions',age_Perm,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Perm)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Triassic data

%load in the numerical data
num_Trias = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Triassic','Range','A2:H1922');
num_Trias(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Trias = num_Trias(:,5);
age_Trias = num_Trias(:,6);
formgp_Trias = num_Trias(:,1) + max(formgp_Perm);
paleolat_Trias = num_Trias(:,2);

%load in the categorical data:
cat_Trias = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Triassic','Range','K2:N1922');
cat_Trias(:,2:3) = []; %delete the non-categorical columns
dep_Trias = categorical(cat_Trias(:,1));
min_Trias = categorical(cat_Trias(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Trias,formgp_Trias,'Positions',age_Trias,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Trias)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Jurassic data

%load in the numerical data
num_Juras = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Jurassic','Range','A2:H865');
num_Juras(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Juras = num_Juras(:,5);
age_Juras = num_Juras(:,6);
formgp_Juras = num_Juras(:,1) + max(formgp_Trias);
paleolat_Juras = num_Juras(:,2);

%load in the categorical data:
cat_Juras = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Jurassic','Range','K2:N865');
cat_Juras(:,2:3) = []; %delete the non-categorical columns
dep_Juras = categorical(cat_Juras(:,1));
min_Juras = categorical(cat_Juras(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Juras,formgp_Juras,'Positions',age_Juras,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Juras)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Cretaceous data

%load in the numerical data
num_Cret = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Cretaceous','Range','A2:H498');
num_Cret(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Cret = num_Cret(:,5);
age_Cret = num_Cret(:,6);
formgp_Cret = num_Cret(:,1) + max(formgp_Juras);
paleolat_Cret = num_Cret(:,2);

%load in the categorical data:
cat_Cret = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Cretaceous','Range','K2:N498');
cat_Cret(:,2:3) = []; %delete the non-categorical columns
dep_Cret = categorical(cat_Cret(:,1));
min_Cret = categorical(cat_Cret(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Cret,formgp_Cret,'Positions',age_Cret,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Cret)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Paleogene data

%load in the numerical data
num_Pg = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Paleogene','Range','A2:H292');
num_Pg(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Pg = num_Pg(:,5);
age_Pg = num_Pg(:,6);
formgp_Pg = num_Pg(:,1) + max(formgp_Cret);
paleolat_Pg = num_Pg(:,2);

%load in the categorical data:
cat_Pg = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Paleogene','Range','K2:N292');
cat_Pg(:,2:3) = []; %delete the non-categorical columns
dep_Pg = categorical(cat_Pg(:,1));
min_Pg = categorical(cat_Pg(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Pg,formgp_Pg,'Positions',age_Pg,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Pg)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Neogene data

%load in the numerical data
num_Ng = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Neogene','Range','A2:H301');
num_Ng(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Ng = num_Ng(:,5);
age_Ng = num_Ng(:,6);
formgp_Ng = num_Ng(:,1) + max(formgp_Pg);
paleolat_Ng = num_Ng(:,2);

%load in the categorical data:
cat_Ng = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Neogene','Range','K2:N301');
cat_Ng(:,2:3) = []; %delete the non-categorical columns
dep_Ng = categorical(cat_Ng(:,1));
min_Ng = categorical(cat_Ng(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Ng,formgp_Ng,'Positions',age_Ng,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Ng)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%load in and organize Quaternary data

%load in the numerical data
num_Quat = readmatrix('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Quaternary','Range','A2:H749');
num_Quat(:,[3,7]) = []; %delete the text columns
%columns:
    %1: formation group
    %2: paleolatitude
    %2: grain number
    %3: major axis (um)
    %4: minor axis (um) **USE THIS ONE**
    %5: plot age (Ma) - midpoint of the given age bin
%end
D_Quat = num_Quat(:,5);
age_Quat = num_Quat(:,6);
formgp_Quat = num_Quat(:,1) + max(formgp_Ng);
paleolat_Quat = num_Quat(:,2);

%load in the categorical data:
cat_Quat = readcell('Ooid Size_Supplementary Material.xlsx','Sheet',...
    'Quaternary','Range','K2:N749');
cat_Quat(:,2:3) = []; %delete the non-categorical columns
dep_Quat = categorical(cat_Quat(:,1));
min_Quat = categorical(cat_Quat(:,2));

%%
%plot data to make sure it looks ok
hold on
boxplot(D_Quat,formgp_Quat,'Positions',age_Quat,'Widths',5,...
    'PlotStyle','compact','Symbol','.','ColorGroup',min_Quat)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')

%%
%concatenate all data
D_all = cat(1,D_Cambrian,D_Ord,D_Sil,D_Dev,D_Carb,D_Perm,D_Trias,D_Juras,...
    D_Cret,D_Pg,D_Ng,D_Quat);
age_all = cat(1,age_Cambrian,age_Ord,age_Sil,age_Dev,age_Carb,age_Perm,age_Trias,age_Juras,...
    age_Cret,age_Pg,age_Ng,age_Quat);
formgp_all = cat(1,formgp_Cambrian,formgp_Ord,formgp_Sil,formgp_Dev,formgp_Carb,formgp_Perm,formgp_Trias,formgp_Juras,...
    formgp_Cret,formgp_Pg,formgp_Ng,formgp_Quat);
dep_all = cat(1,dep_Cambrian,dep_Ord,dep_Sil,dep_Dev,dep_Carb,dep_Perm,dep_Trias,dep_Juras,...
    dep_Cret,dep_Pg,dep_Ng,dep_Quat);
min_all = cat(1,min_Cambrian,min_Ord,min_Sil,min_Dev,min_Carb,min_Perm,min_Trias,min_Juras,...
    min_Cret,min_Pg,min_Ng,min_Quat);
paleolat_all = cat(1,paleolat_Cambrian,paleolat_Ord,paleolat_Sil,paleolat_Dev,...
    paleolat_Carb,paleolat_Perm,paleolat_Trias,paleolat_Juras,...
    paleolat_Cret,paleolat_Pg,paleolat_Ng,paleolat_Quat);

%%
%plot overview figure to check out data

%note that for some reason some of the boxes don't end up at quite the
%right x values, so this plot is only good for generally checking that
%eerything looks fine, use the 2nd output figure for correct x axis
%locations
f_output = figure;
subplot(2,1,1)
boxplot(D_all,formgp_all,'Positions',age_all,'Widths',5,'PlotStyle',...
    'compact','Symbol','.','ColorGroup',min_all)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
ylim([0 1500])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')
title('grouped by mineralogy')

subplot(2,1,2)
boxplot(D_all,formgp_all,'Positions',age_all,'Widths',5,'PlotStyle',...
    'compact','Symbol','.','ColorGroup',dep_all)
set(gca,'XTickLabel',{' '})
xticks(0:100:500)
xticklabels('auto')
xlim([0 550])
ylim([0 1500])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')
title('grouped by dep env')

%%
%simplify data for each formation group
%preallocate space in vectors
n_fmgroups = max(formgp_all);
D50_data = zeros(n_fmgroups,1);
D60_data = zeros(n_fmgroups,1);
D75_data = zeros(n_fmgroups,1);
D_stdev = zeros(n_fmgroups,1);
mineral = categorical(zeros(n_fmgroups,1));
depenv = categorical(zeros(n_fmgroups,1));
plot_age = zeros(n_fmgroups,1);
paleolat = zeros(n_fmgroups,1);

for n = 1:n_fmgroups
    mintemp = min_all(formgp_all == n);
    mineral(n) = mintemp(1);
    depenvtemp = dep_all(formgp_all == n);
    depenv(n) = depenvtemp(1);
    plot_agetemp = age_all(formgp_all == n);
    plot_age(n) = plot_agetemp(1);
    paleolat_temp = paleolat_all(formgp_all == n);
    paleolat(n) = paleolat_temp(1);
    D50_data(n) = prctile(D_all(formgp_all == n),50);
    D60_data(n) = prctile(D_all(formgp_all == n),60);
    D75_data(n) = prctile(D_all(formgp_all == n),75);
    D_stdev(n) = std(D_all(formgp_all == n));
end

%apply correction to measured b50
k2 = 1; %assume that ooids are radially symmetric about major axis
D50_corr = ooidsizecorrection(D50_data,k2);
D_range = D60_data - D50_data;

%%
f_outputgrouped = figure;

colors_comb = magma(4);

errorbar(plot_age(mineral == 'aragonite' & depenv == 'shoal'),...
    D50_corr(mineral == 'aragonite' & depenv == 'shoal'),...
    D_range(mineral == 'aragonite' & depenv == 'shoal'),...
    'o','Color',colors_comb(1,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(1,:))
hold on
errorbar(plot_age(mineral == 'aragonite' & depenv == 'tidal'),...
    D50_corr(mineral == 'aragonite' & depenv == 'tidal'),...
    D_range(mineral == 'aragonite' & depenv == 'tidal'),...
    'd','Color',colors_comb(1,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(1,:))
errorbar(plot_age(mineral == 'aragonite' & depenv == 'high_energy'),...
    D50_corr(mineral == 'aragonite' & depenv == 'high_energy'),...
    D_range(mineral == 'aragonite' & depenv == 'high_energy'),...
    's','Color',colors_comb(1,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(1,:))
errorbar(plot_age(mineral == 'aragonite' & depenv == 'low_energy'),...
    D50_corr(mineral == 'aragonite' & depenv == 'low_energy'),...
    D_range(mineral == 'aragonite' & depenv == 'low_energy'),...
    'v','Color',colors_comb(1,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(1,:))
errorbar(plot_age(mineral == 'aragonite' & depenv == 'other'),...
    D50_corr(mineral == 'aragonite' & depenv == 'other'),...
    D_range(mineral == 'aragonite' & depenv == 'other'),...
    'p','Color',colors_comb(1,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(1,:))

errorbar(plot_age(mineral == 'calcite' & depenv == 'shoal'),...
    D50_corr(mineral == 'calcite' & depenv == 'shoal'),...
    D_range(mineral == 'calcite' & depenv == 'shoal'),...
    'o','Color',colors_comb(2,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(2,:))
errorbar(plot_age(mineral == 'calcite' & depenv == 'tidal'),...
    D50_corr(mineral == 'calcite' & depenv == 'tidal'),...
    D_range(mineral == 'calcite' & depenv == 'tidal'),...
    'd','Color',colors_comb(2,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(2,:))
errorbar(plot_age(mineral == 'calcite' & depenv == 'high_energy'),...
    D50_corr(mineral == 'calcite' & depenv == 'high_energy'),...
    D_range(mineral == 'calcite' & depenv == 'high_energy'),...
    's','Color',colors_comb(2,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(2,:))
errorbar(plot_age(mineral == 'calcite' & depenv == 'low_energy'),...
    D50_corr(mineral == 'calcite' & depenv == 'low_energy'),...
    D_range(mineral == 'calcite' & depenv == 'low_energy'),...
    'v','Color',colors_comb(2,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(2,:))
errorbar(plot_age(mineral == 'calcite' & depenv == 'other'),...
    D50_corr(mineral == 'calcite' & depenv == 'other'),...
    D_range(mineral == 'calcite' & depenv == 'other'),...
    'p','Color',colors_comb(2,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(2,:))

errorbar(plot_age(mineral == 'bimineralic' & depenv == 'shoal'),...
    D50_corr(mineral == 'bimineralic' & depenv == 'shoal'),...
    D_range(mineral == 'bimineralic' & depenv == 'shoal'),...
    'o','Color',colors_comb(3,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(3,:))
errorbar(plot_age(mineral == 'bimineralic' & depenv == 'tidal'),...
    D50_corr(mineral == 'bimineralic' & depenv == 'tidal'),...
    D_range(mineral == 'bimineralic' & depenv == 'tidal'),...
    'd','Color',colors_comb(3,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(3,:))
errorbar(plot_age(mineral == 'bimineralic' & depenv == 'high_energy'),...
    D50_corr(mineral == 'bimineralic' & depenv == 'high_energy'),...
    D_range(mineral == 'bimineralic' & depenv == 'high_energy'),...
    's','Color',colors_comb(3,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(3,:))
errorbar(plot_age(mineral == 'bimineralic' & depenv == 'low_energy'),...
    D50_corr(mineral == 'bimineralic' & depenv == 'low_energy'),...
    D_range(mineral == 'bimineralic' & depenv == 'low_energy'),...
    'v','Color',colors_comb(3,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(3,:))
errorbar(plot_age(mineral == 'bimineralic' & depenv == 'other'),...
    D50_corr(mineral == 'bimineralic' & depenv == 'other'),...
    D_range(mineral == 'bimineralic' & depenv == 'other'),...
    'p','Color',colors_comb(3,:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',colors_comb(3,:))

xlim([0 550])
ylim([0 2000])
set(gca, 'xdir', 'reverse')
box on
xlabel('age (Ma)')
ylabel('D_5_0 (\mum)')

%%
%save variables that can be fed into the PhanOomega code
save('AKdata_formation.mat','D50_corr','D_stdev','D_range','D50_data',...
    'D75_data',"plot_age",'depenv','mineral','paleolat');