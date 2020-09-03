%%-----------------------------------------------------------------------%%
% This script constructs univariate GPRs from the  
% Created by: David Avendano-Valencia - January 2020
%%-----------------------------------------------------------------------%%

clear
close all
clc

addpath('C:\Users\ldav\Documents\GitHub\Nonlinear Regression AVE\Core\')

%-- Properties of the DWM-FAST simulation
component = 'blade';                                                        % Select the type of component to build GPR on DELs
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

%% Computation loop

ind_d = 3;
ind_wexp = 1;

disp(['Calculating GPR for ',component,...
    ' WExp:',num2str(w_exp(ind_wexp)),...
    ' Dist:',num2str(distances(ind_d))])

m = 10*(1:20);
theta = zeros(6,10,numel(m));
lnL = zeros(10,numel(m));
t = zeros(1,numel(m));


for n=1:numel(m)
    tic
    [theta(:,:,n),lnL(:,n)] = ObtainGPRforData(component,ind_wexp,ind_d,m(n));
    t(n) = toc;
end

save(['Results Sparse\NoVarVSPerformance_',component],'theta','lnL','t')

%%
close all
clc

% load(['Results Sparse\NoVarVSPerformance_',component],'theta','lnL','t')

plot(m,mean(lnL))


figure
for i=1:6
    subplot(2,3,i)
    semilogy(m,squeeze(median(theta(i,:,:),2)))
    grid on
end

%% ------------------------------------------------------------------------
function [theta,lnL] = ObtainGPRforData(component,ind_w_exp,ind_d,N)

%-- Defining data folder
data_folder = 'DWM Datasets\';

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
x = cell(4,1);
y = cell(4,1);
for i=1:2
    for j=1:2
        ind_turbine = wakeData.Turbine == i;
        x{i+2*(j-1)} = X(ind_turbine&ind_dist,:)';
        y{i+2*(j-1)} = Y(indices,:);
        y{i+2*(j-1)} = y{i+2*(j-1)}(ind_turbine&ind_dist,j)';
    end
end

x = x{1};
y = y{1};

%-- Markov Chain Monte Carlo optimization of the GPR via Metropolis Hastings sampling
theta0 = ones(1,size(X,2)+2);                                               % Prior parameters

theta = zeros(6,10);
lnL = zeros(10,1);


parfor i=1:10
    
    %-- Obtaining a sample for Subset of Regressors
    ind = UniformSpaceSampling( x, N );
    indx = false(1,size(x,2));
    indx(ind) = true;
    
    %-- Sampling from the GPR
    [theta(:,i),lnL(i)] = optimize_gpr( x, y, theta0, 'SoR', indx );
    
end

end
