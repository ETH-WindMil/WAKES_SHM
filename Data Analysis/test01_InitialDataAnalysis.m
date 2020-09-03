%%-----------------------------------------------------------------------%%
% This script provides an initial exploratory analysis on the data obtained
% from the FAST-DWM simulations.
% Created by: David Avendano-Valencia - January 2020
%%-----------------------------------------------------------------------%%

clear
close all
clc

%-- Defining data folder
data_folder = 'DWM Datasets\';

%-- Loading wake data
load([data_folder,'DWMwakeInputData'],'wakeData');

%-- Loading Damage Equivalent Loads
load([data_folder,'BladeFatigueDataTable'],'BladeFatigueDataTable');
FatigueData = BladeFatigueDataTable;
clear BladeFatigueDataTable

w_exp = unique(FatigueData.m);

%% Extracting input and output data
close all
clc

ind_w_exp = 2;
indices = FatigueData.m == w_exp(ind_w_exp);

EOPnames = wakeData.Properties.VariableNames(4:8);
X = [wakeData.meanU wakeData.Ti wakeData.alpha wakeData.HFlowAng wakeData.distBetweenWTG];
Y = [FatigueData.RootMxb1(indices) FatigueData.RootMyb1(indices)];

%%
close all
clc

ind_turbine = wakeData.Turbine == 1;

for j=1:2
    figure('Position',[100+(j-1)*900 100 900 600])
    for i=1:4
        subplot(2,2,i)
        plot(X(ind_turbine,i),Y(ind_turbine,j),'o')
        hold on
        plot(X(~ind_turbine,i),Y(~ind_turbine,j),'s')
        xlabel(EOPnames(i))
        ylabel('DEL')
        grid on
        legend({'Up-wind','Waked'},'Location','northoutside','Orientation','horizontal')
    end
end
