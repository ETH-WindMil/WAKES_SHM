clear
close all
clc

%-- Add GPR functions to matlab's path
addpath('Core')

%-- Properties of the DWM-FAST simulation
component = 'tower';                                                        % Select the type of component to build GPR on DELs
                                                                            % 'tower' - Tower base bending moment in the fore-aft and lateral directions
                                                                            % 'blade' - Blade root bending moment in the flapwise and edgewise directions
                                                                            % 'lss'   - Low speed shaft bending moment in the two directions perpendicular to rotation axis
                                                                            % 'yaw'   - Yaw mechanism bending moment in the fore-aft and lateral directions

%-- Wohler exponents for each wind turbine component type
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

%-- Distance between wind turbines
distances = [3 5 8 11];
N_dist = 4;

%-- Hyperparameter names
HyperParNames = {'\sigma_u^2';                                              % Names of GPR hyperparameters
                 '\sigma_f^2';
                 '\lambda_1';
                 '\lambda_2';
                 '\lambda_3';
                 '\lambda_4';};
             
%-- Output names
switch component
    case 'blade'
        OutputNames = {'Edgewise bending moment - Upwind';
                       'Edgewise bending moment - Downwind';
                       'Flapwise bending moment - Upwind';
                       'Flapwise bending moment - Downwind'};
    case 'tower'
        OutputNames = {'Side-to-side bending moment - Upwind';
                       'Side-to-side bending moment - Downwind';
                       'Fore-aft bending moment - Upwind';
                       'Fore-aft bending moment - Downwind'};
                   
	case 'yaw'
        OutputNames = {'Yaw bearing pitch moment - Upwind';
                       'Yaw bearing yaw moment - Downwind';
                       'Yaw bearing pitch moment - Upwind';
                       'Yaw bearing yaw moment - Downwind'};
    case 'lss'
        OutputNames = {'LSS bending moment Y axis - Upwind';
                       'LSS bending moment Z axis - Downwind';
                       'LSS bending moment Y axis - Upwind';
                       'LSS bending moment Z axis - Downwind'};
end
             
%% Calculating the prediction error based on the optimized GPR

%-- Load optimization results
ind_wexp = 2;                                                               % Wohler exponent index
ind_d = 1;                                                                  % Distance index
load(['Results Sparse/UniGPR_',component,...
            '_WExp_',num2str(w_exp(ind_wexp)),...
            '_Dist_',num2str(distances(ind_d))],'theta')

%-- Loading DEL dataset
[x,y,~] = LoadDWMData(component,ind_wexp,ind_d);

rss_sss = zeros(4,2);

for i=1:2
    X = x{2*i};
    Y = y{2*i};
    N = size(X,2);
    
    %-- Obtaining a sample for Subset of Regressors
    ind = UniformSpaceSampling( X, 180 );
    indx = false(size(X,2),1);
    indx(ind) = true;
    
    %-- Predicting DELs in the input grid
    f_hat = gpr_predict(X,X,Y,theta(:,i),'SoD',indx);
    
    rss_sss(i,1) = sum( ( f_hat(indx) - Y(indx) ).^2 )./ sum( Y(indx).^2 );
    rss_sss(i,2) = sum( ( f_hat(~indx) - Y(~indx) ).^2 )./ sum( Y(~indx).^2 );
        
end

%-- Load optimization results
load(['Results Up2Dn/UniGPR_Up2Dn',component,...
            '_WExp_',num2str(w_exp(ind_wexp)),...
            '_Dist_',num2str(distances(ind_d))],'theta')

%-- Loading DEL dataset
[x,y] = LoadDWMData(component,ind_wexp,ind_d);

for i=1:2
    X = x{2*i-1};
    Y = y{2*i};
    N = size(X,2);
    
    %-- Obtaining a sample for Subset of Regressors
    ind = UniformSpaceSampling( X, 180 );
    indx = false(size(X,2),1);
    indx(ind) = true;
    
    %-- Predicting DELs in the input grid
    f_hat = gpr_predict(X,X,Y,theta(:,i),'SoD',indx);
    
    rss_sss(i+2,1) = sum( ( f_hat(indx) - Y(indx) ).^2 )./ sum( Y(indx).^2 );
    rss_sss(i+2,2) = sum( ( f_hat(~indx) - Y(~indx) ).^2 )./ sum( Y(~indx).^2 );
        
end

%%
close all
clc

figure('position',[100 400 900 280])

barh(log10(rss_sss))
xlabel('$\log_{10} RSS/SSS$','Interpreter','latex')
legend({'Training';'Validation'},'Location','northoutside','Orientation','horizontal')
set(gca,'YTickLabel',{OutputNames{2} OutputNames{4}})
grid on
xlim([-4.5 0])