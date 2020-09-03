clear
close all
clc

%-- Add GPR functions to matlab's path
addpath('C:\Users\ldav\Documents\GitHub\Nonlinear Regression AVE\Core')

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
load(['Results Sparse\UniGPR_',component,...
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
                   
        OutNames = {'Edgewise bending moment';
                    'Flapwise bending moment'};
    case 'tower'
        OutputNames = {'Side-to-side bending moment - Upwind';
                       'Side-to-side bending moment - Downwind';
                       'Fore-aft bending moment - Upwind';
                       'Fore-aft bending moment - Downwind'};
                   
        OutNames = {'Side-to-side bending moment';
                    'Fore-aft bending moment'};
                   
	case 'yaw'
        OutputNames = {'Yaw bearing pitch moment - Upwind';
                       'Yaw bearing pitch moment - Downwind';
                       'Yaw bearing yaw moment - Upwind';
                       'Yaw bearing yaw moment - Downwind'};
                   
        OutNames = {'Yaw bearing pitch moment';
                    'Yaw bearing yaw moment'};
    case 'lss'
        OutputNames = {'LSS bending moment Y axis - Upwind';
                       'LSS bending moment Y axis - Downwind';
                       'LSS bending moment Z axis - Upwind';
                       'LSS bending moment Z axis - Downwind'};
                   
        OutNames = {'LSS bending moment Y axis';
                    'LSS bending moment Z axis'};
                   
end
             
%% Calculating the response surface based on the optimized GPR
close all
clc

%-- Loading training dataset
[x,y,rangeX] = LoadDWMData(component,ind_wexp,ind_d);

n = 30;
f_hat = zeros(n^4,4);
sigmaY2 = zeros(n^4,4);

for i=1:4
    
    X = x{i};
    Y = y{i};
    N = size(X,2);
    
    %-- Calculating a grid of wind speeds and turbulence intensity values
    x_grid = linspace(0,1,n)';
    [X1_grid,X2_grid,X3_grid,X4_grid] = ndgrid( x_grid );
    Xast = [X1_grid(:) X2_grid(:) X3_grid(:) X4_grid(:)];
    
    %-- Obtaining a sample for Subset of Data
    ind = UniformSpaceSampling( x{i}, 180 );
    indx = false(size(x{i},2),1);
    indx(ind) = true;
    
    %-- Predicting DELs in the input grid
    [f_hat(:,i),sigmaY2(:,i)] = gpr_predict(Xast',X,Y,theta(:,i),'SoD',indx);
        
end

save(['Results\ResponseSimulation_',component],'f_hat','sigmaY2','x_grid','n','rangeX')

%% Calculating the integrated KL divergence for each input
close all
clc

load(['Results\ResponseSimulation_',component],'f_hat','sigmaY2','x_grid','n','rangeX')

%-- Names of the inputs
InputNames = {'Wind speed [m/s]';
              'Turbulence intensity [%]';
              'Shear exponent';
              'Inflow horizontal skewness (deg)'};

%-- Calculating the KL divergence for the complete simulation
err = f_hat(:,[2 4]) - f_hat(:,[1 3]);
err = err.^2./sigmaY2(:,[1 3]);
dKL = ( ( sigmaY2(:,[2 4]).\sigmaY2(:,[1 3]) ) + err + log( sigmaY2(:,[2 4])./sigmaY2(:,[1 3]) ) - 1 )/2;

DKL{1} = reshape(dKL(:,1),n,n,n,n);                                         % KL divergence on sensor 1
DKL{2} = reshape(dKL(:,2),n,n,n,n);                                         % KL divergence on sensor 2
                                                                            % ( Sensor depends on the selected component )
input_ind = [1 2 3 4];

figure('Position',[100 100 900 600])

for j=1:4
    
    Xg = rangeX(2,j)*x_grid+rangeX(1,j);                                    % Input vector
    indx_selector = true(1,4);                                              % Selector of dimensions to integrate
    indx_selector(j) = false;
    
    subplot(2,2,j)
    
    for i=1:2
        
        semilogy( Xg, squeeze( mean( DKL{i}, input_ind(indx_selector) ) ) )
        grid on
        hold on
        
    end
    
    ylabel('$D_{KL}( \mathcal{M}_1 || \mathcal{M}_2 )$','Interpreter','latex')
    xlabel(InputNames{j})
    ylim([1e-2 1e1])
    legend(OutNames,'Location','northoutside')
    
end

%% Plotting results
close all
clc

X1g = rangeX(2,1)*x_grid(:,1)+rangeX(1,1);
X2g = rangeX(2,2)*x_grid(:,2)+rangeX(1,2);

clr = lines(4);

%-- Plotting the estimated response surfaces

figure('Position',[100 100 900 800])
for i=1:4
    
    %-- Plotting
    subplot(2,2,i)
    
    mask = ones(size(f_hat(:,i)));
    mask(log(sigmaY2(:,i)) > 4+min(log(sigmaY2(:,i)))) = NaN;
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

%% Plotting the error surfaces
close all
clc

figure('Position',[100 100 900 450])
figure('Position',[1000 100 900 450])
figure('Position',[1000 550 900 450])
for i=1:2
    
    mask = ones(size(f_hat(:,i)));
    mask( log10( sigmaY2(:,2*i) ) > 2+min(log10( sigmaY2(:,2*i) )) ) = 0;
    mask( log10( sigmaY2(:,2*i-1) ) > 2+min(log10( sigmaY2(:,2*i-1) )) ) = 0;
    
    err = f_hat(:,2*i) - f_hat(:,2*i-1);
    err(~logical(mask)) = max(err);
    Z = reshape(err,n,n);
    cmap = [parula(256); ones(1,3)];
    
    figure(1)
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
    cbar.Label.String = 'Pointwise difference';
    xlim([0 25])
    ylim([10 50])
    legend([D1,D2],{'Data points upwind';'Data points downwind'},'Location','southeast')
    
    
    mask = ones(size(f_hat(:,i)));
    mask( log10( sigmaY2(:,2*i) ) > 2+min(log10( sigmaY2(:,2*i) )) ) = inf;
    mask( log10( sigmaY2(:,2*i-1) ) > 2+min(log10( sigmaY2(:,2*i-1) )) ) = inf;
    
    err = f_hat(:,2*i)-f_hat(:,2*i-1);
    err = err.^2./(sigmaY2(:,2*i)+sigmaY2(:,2*i-1));
    Z = reshape(err.*mask,n,n);
    
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
    cbar.Label.String = 'Normalized mean difference';
    xlim([0 25])
    ylim([10 50])
    legend([D1,D2],{'Data points upwind';'Data points downwind'},'Location','southeast')
    
    err = f_hat(:,2*i)-f_hat(:,2*i-1);
    err = err.^2./sigmaY2(:,2*i-1);
    dKL = ( ( sigmaY2(:,2*i).\sigmaY2(:,2*i-1) ) + err + log( sigmaY2(:,2*i)./sigmaY2(:,2*i-1) ) - 1 )/2;
    Z = reshape(dKL.*mask,n,n);
    
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
    cbar.Label.String = 'Kullback-Leibler divergence';
    xlim([0 25])
    ylim([10 50])
    legend([D1,D2],{'Data points upwind';'Data points downwind'},'Location','southeast')
    
end
