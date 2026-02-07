clear; clc;
% ================= Folder =================
folderName = uigetdir(pwd,'选择包含 CSV 文件的文件夹');
if folderName == 0
    return;
end
% ================= parameters =================
px_per_um = 6.15; % pixel to um unit
Ninterp  = 500;
Nmode    = 30;
maxFiles = 20;
% ================= read CSV files =================
files = dir(fullfile(folderName,'*.csv'));
if isempty(files)
    error('文件夹中没有 CSV 文件。');
end
% --- 按文件名数字排序 ---
fileNames = {files.name};
fileNums = zeros(size(fileNames));
for k = 1:numel(fileNames)
    [~, name, ~] = fileparts(fileNames{k});
    tmp = str2double(name);
    if isnan(tmp); tmp = Inf; end
    fileNums(k) = tmp;
end
[~, idxSort] = sort(fileNums);
files = files(idxSort);

nFiles = min(maxFiles, numel(files));

% ================= Start =================
A_mat = NaN(Nmode,nFiles);
q_mat = NaN(Nmode,nFiles);

FiberLength = NaN(nFiles,1);
qc_all      = NaN(nFiles,1);   
TotalEnergy = NaN(nFiles,1);    
BRI         = NaN(nFiles,1);   
q_eff       = NaN(nFiles,1);

% ================= Main loop =================
for k = 1:nFiles

    data = readmatrix(fullfile(folderName,files(k).name));
    if size(data,2) < 2, continue; end

    % -------- cooridnate --------
    x = data(:,1)/px_per_um;
    y = data(:,2)/px_per_um;
    valid = ~(isnan(x)|isnan(y));
    x = x(valid); 
    y = y(valid);
    if numel(x) < 10, continue; end

    % -------- Fiber length --------
    ds = hypot(diff(x),diff(y));
    s  = [0; cumsum(ds)];
    L  = s(end);
    FiberLength(k) = L;
    if L <= 0, continue; end

    % -------- Interp --------
    s_u = linspace(0,L,Ninterp)';
    xu  = interp1(s,x,s_u,'pchip');
    yu  = interp1(s,y,s_u,'pchip');

    % -------- tangiential angle --------
    dx = gradient(xu,s_u);
    dy = gradient(yu,s_u);
    theta = unwrap(atan2(dy,dx));
    theta = theta - mean(theta);

    % -------- Fourier decomposition --------
    q  = (1:Nmode)' * pi / L;
    Aq = zeros(Nmode,1);
    for n = 1:Nmode
        Aq(n) = sqrt(2/L) * trapz(s_u, theta .* cos(q(n)*s_u));
    end

    A_mat(:,k) = Aq;
    q_mat(:,k) = q;

    % -------- qc--------
    qc_all(k) = 2*pi/L;

    % ======== bending energy ========
    TotalEnergy(k) = sum(abs(Aq).^2);
    BRI(k) = TotalEnergy(k) / L;

    % -------- characteristic q: q_eff--------
    E_q = q.^2 .* abs(Aq).^2;
    if sum(E_q) > 0
        q_eff(k) = sum(q .* E_q) / sum(E_q);
    end
end

% ================= Figure =================
figure('Color','w','Position',[100 100 1400 600]);

% -------- Amplitude vs wavenumber --------
subplot(1,2,1); hold on;

Z_all = abs(A_mat(:));
Zmin = 0;
Zmax = prctile(Z_all,95);

Yvals = [0 20 40 60 80 100 120];

for k = 1:nFiles
    q = q_mat(:,k);
    Z = abs(A_mat(:,k));

    X = [q q];
    Y = Yvals(k) * ones(Nmode,1);
    ZZ = [Z Z];
    C  = ZZ;

    surface(X,Y,ZZ,C, ...
        'FaceColor','none', ...
        'EdgeColor','interp', ...
        'LineWidth',2);
end

xlabel('q (\mum^{-1})');
ylabel('Fiber index');
zlabel('|A_q| (\mum^{1/2})');
title('Amplitude spectrum');

colormap(jet);
cb = colorbar;
ylabel(cb,'|A_q| (\mum^{1/2})');
clim([Zmin Zmax]);
zlim([Zmin Zmax]);

grid on;
view(45,30);

% -------- BRI --------
subplot(1,2,2);
bar(BRI);
xlabel('Fiber index');
ylabel('BRI = E / L (\mum^{-1})');
title('Total spectral bending index');
grid on;

% ================= save as Excel =================
ResultTable = table( ...
    {files(1:nFiles).name}', ...
    FiberLength, ...
    qc_all, ...
    q_eff, ...
    TotalEnergy, ...
    BRI, ...
    'VariableNames', ...
    {'FileName','FiberLength_um','qc_umInv','q_eff_umInv','TotalEnergy','BRI'} );

outFile = fullfile(folderName,'BendingResistanceIndex.xlsx');
writetable(ResultTable,outFile);

fprintf('结果已保存到: %s\n', outFile);

