%%-----------------------------------------------------------------------%%
% This script provides an exploratory exceedence probabilities of the DEL 
% data obtained from the FAST-DWM simulations. Comparing the DEL of the
% upwind and wake-affected wind turbines conditional on input environmental
% inflow conditions.
% Created by: Imad Abdallah - Sept 2020
%%-----------------------------------------------------------------------%%

clear
close all
clc

%%
% filter

% Mean wind speed range for upwind wind turbine:
Vrange = [3 10];
% Vrange = [11 14];
% Vrange = [15 25];

% Turbulence intensity range for upwind wind turbine:
Tirange = [0 50];

% Shear exponent range for upwind wind turbine:
Shearrange = [-2 2];

% Inflow skewnewss range for upwind wind turbine:
skewrange = [ -10 -1 ];
% skewrange = [ -0.9999 0.9999 ];
% skewrange = [ 1 10 ];

% Distace range between turbines
dist_betw_wtg = [3 5];
% dist_betw_wtg = [8 11];

sensor_component = 'TOWER'; % 'BLADE', 'TOWER', 'YAW'
ChooseSensorNumber = 1; % 1,2

% List of sensors


%%
%-- Defining data folder
data_folder = '.\DWM_Datasets\';

%-- Loading freestream and inflow wake data
load([data_folder,'DWMwakeInputData'],'wakeData');

if strcmp(sensor_component,'BLADE')
    vars = {'RootMxb1','RootMyb1'}; 
    %-- BLADE: Loading Damage Equivalent Loads
    load([data_folder,'BladeFatigueDataTable'],'BladeFatigueDataTable');
    FatigueData = BladeFatigueDataTable;
    clear BladeFatigueDataTable
end

if strcmp(sensor_component,'TOWER')
    vars = {'TwrBsMxt','TwrBsMyt'};     
    %-- TOWER: Loading Damage Equivalent Loads
    load([data_folder,'TowerFatigueDataTable'],'TowerFatigueDataTable');
    FatigueData = TowerFatigueDataTable;
    clear TowerFatigueDataTable
end

if strcmp(sensor_component,'YAW')
    vars = {'YawBrMyp','YawBrMzp'};     
    %-- YAW: Loading Damage Equivalent Loads
    load([data_folder,'YawFatigueDataTable'],'YawFatigueDataTable');
    FatigueData = YawFatigueDataTable;
    clear YawFatigueDataTable
end


%%
% filter according to inflow conditions at turbine 1 (up-wind wind turbine)
Vidx_upwind = find(wakeData.Turbine==1 & wakeData.meanU>=min(Vrange) & wakeData.meanU<=max(Vrange)... 
& wakeData.Ti>=min(Tirange) & wakeData.Ti<=max(Tirange)...
& wakeData.alpha>=min(Shearrange) & wakeData.alpha<=max(Shearrange)...
& wakeData.HFlowAng>=min(skewrange) & wakeData.HFlowAng<=max(skewrange)...
& wakeData.distBetweenWTG>=min(dist_betw_wtg) & wakeData.distBetweenWTG<=max(dist_betw_wtg));

runNumbers = wakeData.Run(Vidx_upwind);

Vidx_downwind = Vidx_upwind + 1;

figure
histogram(wakeData.meanU(Vidx_upwind),10)
hold on 
histogram(wakeData.meanU(Vidx_downwind),10)
xlabel('U [m/s]','FontSize',16)
set(gcf,'color','w');
grid minor
set(gca,'LineWidth',1,'FontSize',16)
legend({'Up-wind','Wake-affected'},'Location','NorthEast','FontSize',16)

figure
histogram(wakeData.Ti(Vidx_upwind),10)
hold on 
histogram(wakeData.Ti(Vidx_downwind),10)
xlabel('Ti [-]','FontSize',16)
set(gcf,'color','w');
grid minor
set(gca,'LineWidth',1,'FontSize',16)
legend({'Up-wind','Wake-affected'},'Location','NorthEast','FontSize',16)



%%
% Filter fatigue data
rowsTurbine1 = find(ismember(FatigueData.RunNumber,runNumbers)==1 & FatigueData.Turbine==1);% 
rowsTurbine2 = find(ismember(FatigueData.RunNumber,runNumbers)==1 & FatigueData.Turbine==2);% 



%%
% data for turbine 1
TTurbine1 = FatigueData(rowsTurbine1,vars);
lldsTurbine1 = TTurbine1.(genvarname(vars{ChooseSensorNumber}));
lldsTurbine1(isnan(lldsTurbine1))=[];          % remove NaNs
ind=find(lldsTurbine1<=0);lldsTurbine1(ind)=[];
[xfTurbine1, indf] = sort(lldsTurbine1);
n = length(indf);
ffTurbine1 = ((1:n)-0.3)/(n+0.4); %     quantile = (ii-0.3)/(n+0.4); % formula for the median rank
[fTurbine1,xiTurbine1] = ksdensity(xfTurbine1,'support','positive','function','pdf','bandwidth',0.35);

% data for turbine 2
TTurbine2 = FatigueData(rowsTurbine2,vars);
lldsTurbine2 = TTurbine2.(genvarname(vars{ChooseSensorNumber}));
lldsTurbine2(isnan(lldsTurbine2))=[];          % remove NaNs
ind=find(lldsTurbine2<=0);lldsTurbine2(ind)=[];
[xfTurbine2, indf] = sort(lldsTurbine2);
n = length(indf);
ffTurbine2 = ((1:n)-0.3)/(n+0.4); %     quantile = (ii-0.3)/(n+0.4); % formula for the median rank
[fTurbine2,xiTurbine2] = ksdensity(xfTurbine2,'support','positive','function','pdf','bandwidth',0.35);

%%
% plot
reps = length(lldsTurbine1)/length(Vidx_upwind);
figure
plot(repelem(wakeData.HFlowAng(Vidx_upwind),reps),lldsTurbine1,'r+','LineWidth',1.2)
hold on
plot(repelem(wakeData.HFlowAng(Vidx_upwind),reps),lldsTurbine2,'ks');
legend({'Up-wind','Wake-affected'},'Location','NorthWest','FontSize',16)
ylabel(strcat(vars{ChooseSensorNumber},'[kNm]'),'FontSize',16)
xlabel('Horizontal inflow skewness','FontSize',16)
set(gcf,'color','w');
grid minor
set(gca,'LineWidth',1,'FontSize',16)

%%
% plot
figure
semilogy(xfTurbine1,1-ffTurbine1,'r+')
hold on
semilogy(xfTurbine2,1-ffTurbine2,'ks')
legend({'Up-wind','Wake-affected'},'Location','SouthWest','FontSize',16)
xlabel(strcat(vars{ChooseSensorNumber},'[kNm]'),'FontSize',16)
ylabel('Pe','FontSize',16)
set(gcf,'color','w');
grid minor
set(gca,'LineWidth',1,'FontSize',16)


%%



