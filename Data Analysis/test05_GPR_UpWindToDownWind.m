%%-----------------------------------------------------------------------%%
% This script constructs univariate GPRs from the  
% Created by: David Avendano-Valencia - January 2020
%%-----------------------------------------------------------------------%%

clear
close all
clc

addpath('Core')

%-- Properties of the DWM-FAST simulation
component = 'yaw';                                                        % Select the type of component to build GPR on DELs
                                                                            % 'tower' - Tower base bending moment in the fore-aft and lateral directions
                                                                            % 'blade' - Blade root bending moment in the flapwise and edgewise directions
                                                                            % 'lss'   - Low speed shaft bending moment in the two directions perpendicular to rotation axis
                                                                            % 'yaw'   - Yaw mechanism bending moment in the fore-aft and lateral directions

% Wohler exponents for each wind turbine component type
switch component
    case 'blade'
        N_wexp = 5;
        w_exp = 9:13;
    case 'tower'
        N_wexp = 2;
        w_exp = 3:4;
    case 'lss'
        N_wexp = 2;
        w_exp = 6:7;
    case 'yaw'
        N_wexp = 3;
        w_exp = 8:10;
end

% Distance between wind turbines
distances = [3 5 8 11];
N_dist = 4;

%-- Properties of the Metropolis Hastings sampling algorithm for the GPR
nsamples = 1e4;                                                             % Number of samples of MHS

%% Computation loop

for ind_d = 1:N_dist
    for ind_wexp = 1:N_wexp
        disp(['Calculating GPR for ',component,...
            ' WExp:',num2str(w_exp(ind_wexp)),...
            ' Dist:',num2str(distances(ind_d))])
        [theta,smpl] = ObtainGPRforData(component,nsamples,ind_wexp,ind_d);
        save(['Results Up2Dn/UniGPR_Up2Dn',component,...
            '_WExp_',num2str(w_exp(ind_wexp)),...
            '_Dist_',num2str(distances(ind_d))],'theta','smpl')
    end
end

%% Plotting results
close all
clc

component = 'blade';
ind_wexp = 2;
ind_d = 1;
% load(['Results Sparse/UniGPR_',component,...
%             '_WExp_',num2str(w_exp(ind_wexp)),...
%             '_Dist_',num2str(distances(ind_d))],'theta','smpl')

HyperParNames = {'\sigma_u^2';
                 '\sigma_f^2';
                 '\lambda_1';
                 '\lambda_2';
                 '\lambda_3';
                 '\lambda_4';};

figure('Position',[100 100 600 900])
for i=1:6
    subplot(6,1,i)
    plot(log10(squeeze(smpl(:,i,:))))
    grid on
    xlabel('Samples')
    ylabel(HyperParNames{i})
end

figure('Position',[700 100 900 250])
for i=1:2
    subplot(1,2,i)
    boxplot(log10(smpl(:,:,i)))
end

figure('Position',[700 600 900 450])
for i=1:6
    subplot(2,3,i)
    for j=1:2
        histogram(log10(squeeze(smpl(:,i,j))))
        hold on
    end
    xlabel(HyperParNames{i})
    grid on
end

%% ------------------------------------------------------------------------
function [theta,smpl] = ObtainGPRforData(component,nsamples,ind_w_exp,ind_d)

%-- Defining data folder
data_folder = 'DWM Datasets/';

%-- Loading wake data
load([data_folder,'DWMwakeInputData'],'wakeData');

%-- Loading Damage Equivalent Loads of the specified component
switch component
    case 'blade'
        % ---------- Blade root edgewise and flapwise moments -------------
        load([data_folder,'BladeFatigueDataTable'],'BladeFatigueDataTable');
        FatigueData = BladeFatigueDataTable;
        clear BladeFatigueDataTable
        Y = [FatigueData.RootMxb1 FatigueData.RootMyb1]/1e4;
        
    case 'tower'
        % --------- Tower base side-to-side and fore-aft moments ----------
        load([data_folder,'TowerFatigueDataTable'],'TowerFatigueDataTable');
        FatigueData = TowerFatigueDataTable;
        clear TowerFatigueDataTable
        Y = [FatigueData.TwrBsMxt FatigueData.TwrBsMyt]/1e4;
        
    case 'lss'
        % -------------- Rotating LSS bending moments ---------------------
        load([data_folder,'LSSFatigueDataTable'],'LSSFatigueDataTable');
        FatigueData = LSSFatigueDataTable;
        clear LSSFatigueDataTable
        Y = [FatigueData.LSSGagMya FatigueData.LSSGagMza]/1e4;
        
    case 'yaw'
        % --------- Non-rotating yaw-bearing pitch and yaw moment ---------
        load([data_folder,'YawFatigueDataTable'],'YawFatigueDataTable');
        FatigueData = YawFatigueDataTable;
        clear YawFatigueDataTable
        Y = [FatigueData.YawBrMyp FatigueData.YawBrMzp]/1e4;
        
end

%-- Determining available Wohler exponents for the selected dataset
w_exp = unique(FatigueData.m);

%-- Indices to select the Wohler exponent
indices = FatigueData.m == w_exp(ind_w_exp);

%-- Indices of simulations using the same distance between wind turbines
dist = unique(wakeData.distBetweenWTG);
ind_dist = wakeData.distBetweenWTG == dist(ind_d);

%-- Input variables
X = [wakeData.meanU wakeData.Ti wakeData.alpha wakeData.HFlowAng];
X = X - repmat(min(X),size(X,1),1);
X = X ./ repmat(max(X),size(X,1),1);

%-- Slicing data
x = cell(2,1);
y = cell(2,1);
for j=1:2
    %-- Inputs from upwind turbine
    ind_turbine = wakeData.Turbine == 1;
    x{j} = X(ind_turbine&ind_dist,:)';
    
    %-- DELs from downwind turbine
    ind_turbine = wakeData.Turbine == 2;
    y{j} = Y(indices,:);
    y{j} = y{j}(ind_turbine&ind_dist,j)';
end

%-- Markov Chain Monte Carlo optimization of the GPR via Metropolis Hastings sampling
theta0 = ones(1,size(X,2)+2);                                               % Prior parameters

smpl = zeros(nsamples,6,2);
theta = zeros(6,2);

parfor i = 1:2
    
    %-- Obtaining a sample for Subset of Regressors
    ind = UniformSpaceSampling( x{i}, 180 );
    indx = false(size(x{i},2),1);
    indx(ind) = true;
    
    %-- Sampling from the GPR
    smp = optimize_gpr_bayes( x{i}, y{i}, theta0, nsamples, 'SoD', indx );
    smpl(:,:,i) = smp;
    theta(:,i) = median(smp);
end

end
