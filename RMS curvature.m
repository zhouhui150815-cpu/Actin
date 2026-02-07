% =================================================
% Batch calculate RMS curvature for actin stress fibers (CSV)
% =================================================
clear; clc;
% ----------------- setting -----------------
dataFolder = "F:\ZX1 1TB\logistics\Papers\2025 HZ's paper GHZ\Fig.4 FA-actin\Actin\Control\cell1\analysis"; % 改成你的 CSV 文件夹路径
files = dir(fullfile(dataFolder, '*.csv'));  
N_interp = 200;   
numFibers = length(files);
meanKappa = zeros(numFibers,1);
rmsKappa  = zeros(numFibers,1);
fiberNames = cell(numFibers,1);
% ----------------- Main loop -----------------
for i = 1:numFibers
    % Read CSV file
    data = readmatrix(fullfile(dataFolder, files(i).name));
    x = data(:,1);
    y = data(:,2);
    
    % remove NaN
    valid = ~(isnan(x) | isnan(y));
    x = x(valid);
    y = y(valid);
    
    % save filename
    fiberNames{i} = files(i).name;
    
    % 1. calculate contour
    ds = sqrt(diff(x).^2 + diff(y).^2);
    s = [0; cumsum(ds)];
    
    % 2. uniform
    s_uniform = linspace(0, s(end), N_interp);
    x_u = interp1(s, x, s_uniform, 'spline');
    y_u = interp1(s, y, s_uniform, 'spline');
    
    % 3. derivate
    dx = gradient(x_u, s_uniform);
    dy = gradient(y_u, s_uniform);
    ddx = gradient(dx, s_uniform);
    ddy = gradient(dy, s_uniform);
    
    % 4. mean roo
    kappa = abs(dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^(3/2);
    
    % remove outlier
    kappa = kappa(3:end-2);
    
    % 5. stat curvature
    meanKappa(i) = mean(kappa);
    rmsKappa(i) = sqrt(mean(kappa.^2));
end
% -----------------output -----------------
results = table(fiberNames, meanKappa, rmsKappa, ...
    'VariableNames', {'FiberName','MeanCurvature','RMSCurvature'});
disp(results);
% Save result
writetable(results, fullfile(dataFolder, 'FiberCurvatureResults.csv'));
