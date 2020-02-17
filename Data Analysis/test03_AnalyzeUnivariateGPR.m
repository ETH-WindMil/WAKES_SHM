clear
close all
clc

%-- Add GPR functions to matlab's path
addpath('Core\')

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

%-- Load results
ind_wexp = 1;                                                               % Wohler exponent index
ind_d = 2;                                                                  % Distance index
load(['Results\UniGPR_',component,...
            '_WExp_',num2str(w_exp(ind_wexp)),...
            '_Dist_',num2str(distances(ind_d))],'theta','smpl')
nsamples = size(smpl,1);

HyperParNames = {'\sigma_u^2';                                              % Names of GPR hyperparameters
                 '\sigma_f^2';
                 '\lambda_1';
                 '\lambda_2';
                 '\lambda_3';
                 '\lambda_4';};
             
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
             
%% Calculating the response surface based on the optimized GPR
close all
clc

%-- Loading training dataset
[x,y,rangeX] = LoadDWMData(component,ind_wexp,ind_d);

n = 100;
f_hat = zeros(n^2,4);
sigmaY2 = zeros(n^2,4);

for i=1:4
    
    X = x{i};
    Y = y{i};
    N = size(X,2);
    
    %-- Calculating a grid of wind speeds and turbulence intensity values
    x_grid = [linspace(0,1,n)' linspace(0,1,n)'];
    [X1_grid,X2_grid] = meshgrid(x_grid(:,1),x_grid(:,2));
    Xast = [X1_grid(:) X2_grid(:) median(X(:,3))*ones(n^2,1) median(X(:,4))*ones(n^2,1)];
    
    %-- Predicting DELs in the input grid
    [f_hat(:,i),sigmaY2(:,i)] = gpr_predict(X,Y,Xast',theta(:,1));
        
end

%% Plotting results
close all
clc

X1g = rangeX(2,1)*x_grid(:,1)+rangeX(1,1);
X2g = rangeX(2,2)*x_grid(:,2)+rangeX(1,2);

clr = lines(4);


figure('Position',[100 100 900 800])
for i=1:4
    
    %-- Plotting
    subplot(2,2,i)
    
    mask = ones(size(f_hat(:,i)));
    mask(log(sigmaY2(:,i))>-10) = NaN;
    Z = reshape(f_hat(:,i).*mask,n,n);
    cmap = [ones(1,3); parula(256)];
    
    imagesc( X1g, X2g, Z )
    hold on
    plot( rangeX(2,1)*X(1,:)+rangeX(1,1), rangeX(2,2)*X(2,:)+rangeX(1,2), '.','Color',clr(2,:))
    xlabel('Wind speed [m/s]')
    ylabel('Turbulence intensity [%]')
    xlim([0 25])
    ylim([10 50])
    
    colormap(cmap)
    grid on
    title(OutputNames{i})
        
end

figure('Position',[100 100 900 450])
figure('Position',[1000 100 900 450])
figure('Position',[1000 550 900 450])
for i=1:2
    
    err = (f_hat(:,2*i)-f_hat(:,2*i-1))./f_hat(:,2*i-1);
    mask = ones(size(f_hat(:,i)));
    mask( log( sigmaY2(:,2*i) ) > -10 ) = 0;
    err(~logical(mask)) = max(err);
    Z = reshape(100*err,n,n);
    cmap = [parula(256); ones(1,3)];
    
    figure(2)
    subplot(1,2,i)
    imagesc( X1g, X2g, Z )
    colormap(cmap)
    hold on
    D1 = plot(rangeX(2,1)*x{2*i}(1,:)+rangeX(1,1),rangeX(2,2)*x{2*i}(2,:)+rangeX(1,2),...
        '.','Color',clr(2,:));
    D2 = plot(rangeX(2,1)*x{2*i-1}(1,:)+rangeX(1,1),rangeX(2,2)*x{2*i-1}(2,:)+rangeX(1,2),...
        '.','Color',clr(3,:));
    grid on
    xlabel('Wind speed [m/s]')
    ylabel('Turbulence intensity [%]')
    cbar = colorbar('Location','northoutside');
    cbar.Label.String = 'Relative difference [%]';
    xlim([0 25])
    ylim([10 50])
    legend([D1,D2],{'Data points upwind';'Data points downwind'},'Location','southeast')
    
    err = err.^2./(sigmaY2(:,2*i)+sigmaY2(:,2*i-1));
    mask = ones(size(f_hat(:,i)));
    mask( log( sigmaY2(:,2*i) ) > -10 ) = inf;
    Z = reshape(err.*mask,n,n);
    
    figure(3)
    subplot(1,2,i)
    imagesc( X1g, X2g, Z )
    colormap(cmap)
    hold on
    D1 = plot(rangeX(2,1)*x{2*i}(1,:)+rangeX(1,1),rangeX(2,2)*x{2*i}(2,:)+rangeX(1,2),...
        '.','Color',clr(2,:));
    D2 = plot(rangeX(2,1)*x{2*i-1}(1,:)+rangeX(1,1),rangeX(2,2)*x{2*i-1}(2,:)+rangeX(1,2),...
        '.','Color',clr(3,:));
    grid on
    xlabel('Wind speed [m/s]')
    ylabel('Turbulence intensity [%]')
    cbar = colorbar('Location','northoutside');
    cbar.Label.String = 'Normalized mean difference';
    xlim([0 25])
    ylim([10 50])
    legend([D1,D2],{'Data points upwind';'Data points downwind'},'Location','southeast')
    
    err = f_hat(:,2*i)-f_hat(:,2*i-1);
    err = err.^2./sigmaY2(:,2*i-1);
    dKL = ( sigmaY2(:,2*i).\sigmaY2(:,2*i-1) ) + err + log( sigmaY2(:,2*i)./sigmaY2(:,2*i-1) );
    Z = reshape(dKL.*mask,n,n);
    
    figure(4)
    subplot(1,2,i)    
    imagesc( X1g, X2g, Z )
    colormap(cmap)
    hold on
    D1 = plot(rangeX(2,1)*x{2*i}(1,:)+rangeX(1,1),rangeX(2,2)*x{2*i}(2,:)+rangeX(1,2),...
        '.','Color',clr(2,:));
    D2 = plot(rangeX(2,1)*x{2*i-1}(1,:)+rangeX(1,1),rangeX(2,2)*x{2*i-1}(2,:)+rangeX(1,2),...
        '.','Color',clr(3,:));
    grid on
    xlabel('Wind speed [m/s]')
    ylabel('Turbulence intensity [%]')
    cbar = colorbar('Location','northoutside');
    cbar.Label.String = 'Kullback-Leibler divergence';
    xlim([0 25])
    ylim([10 50])
    legend([D1,D2],{'Data points upwind';'Data points downwind'},'Location','southeast')
    
end


%% ------------------------------------------------------------------------
function [x,y,rangeX] = LoadDWMData(component,ind_w_exp,ind_d)

%-- Defining data folder
data_folder = 'D:\Databases\SHM\Simulation\Wake Simulations - Imad DWM\WakeSims\DWM_FAST_PostProc\';

%-- Loading wake data
load([data_folder,'DWMwakeInputData'],'wakeData');

%-- Loading Damage Equivalent Loads
switch component
    case 'blade'
        load([data_folder,'BladeFatigueDataTable'],'BladeFatigueDataTable');
        FatigueData = BladeFatigueDataTable;
        clear BladeFatigueDataTable
        Y = [FatigueData.RootMxb1 FatigueData.RootMyb1]/1e4;
        
    case 'tower'
        load([data_folder,'TowerFatigueDataTable'],'TowerFatigueDataTable');
        FatigueData = TowerFatigueDataTable;
        clear TowerFatigueDataTable
        Y = [FatigueData.TwrBsMxt FatigueData.TwrBsMyt]/1e4;
        
    case 'lss'
        load([data_folder,'LSSFatigueDataTable'],'LSSFatigueDataTable');
        FatigueData = LSSFatigueDataTable;
        clear LSSFatigueDataTable
        Y = [FatigueData.LSSGagMya FatigueData.LSSGagMza]/1e4;
        
    case 'yaw'
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
rangeX = [min(X);
          max(X)-min(X)];
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

end