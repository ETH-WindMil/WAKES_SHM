%%-----------------------------------------------------------------------%%
% This script performs a comparison of the performance figures obtained in
% the Monte-Carlo optimization of the GPR of the DELs in the up-wind and
% wake-affected wind turbines at different components of the wind turbine.
%
% Created by: David Avendano-Valencia - January 2020
%%-----------------------------------------------------------------------%%

clear
close all
clc

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

%% Plotting optimization results
close all
clc

HyperParNames = {'$\sigma_u^2$';
                 '$\sigma_f^2$';
                 '$\lambda_1$';
                 '$\lambda_2$';
                 '$\lambda_3$';
                 '$\lambda_4$';};

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
end

Nburn = 1e3;
ind_wexp = 4;
for ind_d = 1:N_dist
    
    load(['Results\UniGPR_',component,...
        '_WExp_',num2str(w_exp(ind_wexp)),...
        '_Dist_',num2str(distances(ind_d))],'theta','smpl')
    
    
    figure('Position',[700 100 900 520])
    for i=1:4
        subplot(2,2,i)
        boxplot(log10(smpl(Nburn+1:end,:,i)),'Notch','on')
        set(gca,'FontSize',11)
        set(gca,'XTickLabel',HyperParNames,'TickLabelInterpreter','latex')
        grid on
        title(OutputNames{i})
        ylim([-9 3]), set(gca,'YTick',-9:3:3)
        ylabel('Log-hyperparameter value')
    end
end

%% Hyperparameter confidence intervals vs. distance between wind turbines
close all
clc

Nburn = 1e3;
ind_wexp = 1;

theta_mn = zeros(N_dist,4,4);
theta_up = zeros(N_dist,4,4);
theta_dn = zeros(N_dist,4,4);

for ind_d = 1:N_dist
    
    load(['Results\UniGPR_',component,...
        '_WExp_',num2str(w_exp(ind_wexp)),...
        '_Dist_',num2str(distances(ind_d))],'theta','smpl')
    
    for i=1:4
        
        theta_mn(ind_d,:,i) = median( log10(smpl(Nburn+1:end,3:end,i)) );
        theta_up(ind_d,:,i) = prctile( log10(smpl(Nburn+1:end,3:end,i)), 97.5 ) - theta_mn(ind_d,:,i);
        theta_dn(ind_d,:,i) = theta_mn(ind_d,:,i) - prctile( log10(smpl(Nburn+1:end,3:end,i)),  2.5 );
        
    end
end

figure('Position',[100 100 900 600])
for i=1:4
    subplot(2,2,i)
    errorbar(repmat(distances,N_dist,1)',theta_mn(:,:,i),theta_dn(:,:,i),theta_up(:,:,i),...
        'LineWidth',2)
    grid on
    set(gca,'FontSize',11)
    title(OutputNames{i})
    ylim([-11 3]), set(gca,'YTick',-9:3:3)
    ylabel('Log-hyperparameter value')
    xlabel('Distance between wind turbines')
    legend(HyperParNames(3:end),'Interpreter','latex','Orientation','horizontal',...
        'Location','south')
end

%% Hyperparameter confidence intervals vs. Wohler exponent
close all
clc

Nburn = 1e3;
ind_d = 2;

theta_mn = zeros(N_wexp,4,4);
theta_up = zeros(N_wexp,4,4);
theta_dn = zeros(N_wexp,4,4);

for ind_wexp = 1:N_wexp
    
    load(['Results\UniGPR_',component,...
        '_WExp_',num2str(w_exp(ind_wexp)),...
        '_Dist_',num2str(distances(ind_d))],'theta','smpl')
    
    for i=1:4
        
        theta_mn(ind_wexp,:,i) = median( log10(smpl(Nburn+1:end,3:end,i)) );
        theta_up(ind_wexp,:,i) = prctile( log10(smpl(Nburn+1:end,3:end,i)), 97.5 ) - theta_mn(ind_wexp,:,i);
        theta_dn(ind_wexp,:,i) = theta_mn(ind_wexp,:,i) - prctile( log10(smpl(Nburn+1:end,3:end,i)),  2.5 );
        
    end
end

figure('Position',[100 100 900 600])
for i=1:4
    subplot(2,2,i)
    errorbar(repmat(w_exp,N_dist,1)',theta_mn(:,:,i),theta_dn(:,:,i),theta_up(:,:,i),...
        'LineWidth',2)
    grid on
    set(gca,'FontSize',11)
    title(OutputNames{i})
    ylim([-11 3]), set(gca,'YTick',-9:3:3)
    ylabel('Log-hyperparameter value')
    xlabel('Distance between wind turbines')
    legend(HyperParNames(3:end),'Interpreter','latex','Orientation','horizontal',...
        'Location','south')
end