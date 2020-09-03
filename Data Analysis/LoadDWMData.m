%% ------------------------------------------------------------------------
function [x,y,rangeX] = LoadDWMData(component,ind_w_exp,ind_d)

%-- Defining data folder
data_folder = 'DWM Datasets/';

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