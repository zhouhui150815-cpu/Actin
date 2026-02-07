% ============================================================
% SCI-grade Microtubule Fourier Bending Spectrum Analysis
% Ensemble-averaged <|A_q|^2> with log-log slope fitting
% Reference: Gittes 1993, Brangwynne 2007
% Unit: μm
% ============================================================
clear; clc; close all;
% Parameters
Nmode = 30;              % Number of Fourier modes
numFolders = 3;           % Number of datasets
px2um = 1/6.15;           % Pixel -> μm conversion
meanResults = cell(numFolders,1);
qResults = cell(numFolders,1);
%% Loop over folders
for fidx = 1:numFolders
    folderPath = uigetdir(pwd, sprintf('Select folder %d (CSV)', fidx));
    files = dir(fullfile(folderPath, '*.csv'));

    all_Aq2 = [];

    for k = 1:length(files)
        filename = fullfile(folderPath, files(k).name);
        data = readmatrix(filename);

        % --- Load coordinates ---
        x = data(:,1) * px2um;
        y = data(:,2) * px2um;

        % --- Arc length ---
        ds = sqrt(diff(x).^2 + diff(y).^2);
        s = [0; cumsum(ds)];
        L = s(end);

        % Skip too-short filaments
        if L < 2
            continue;
        end

        % --- Resample to uniform arc-length ---
        Npts = length(s);
        s_uniform = linspace(0, L, Npts)';

        x = interp1(s, x, s_uniform, 'linear');
        y = interp1(s, y, s_uniform, 'linear');
        s = s_uniform;

        % --- Tangent angle θ(s) ---
        dx = gradient(x, mean(diff(s)));
        dy = gradient(y, mean(diff(s)));

        theta = atan2(dy, dx);
        theta = unwrap(theta);
        theta = theta - mean(theta);

        % --- Fourier coefficients ---
        qvals = (1:Nmode)' * pi / L;
        Aq = zeros(Nmode,1);

        for n = 1:Nmode
            phi_n = sqrt(2/L) * cos(qvals(n) * s);
            Aq(n) = trapz(s, theta .* phi_n);
        end

        % --- Save |A_q|^2 ---
        Aq2 = abs(Aq).^2;
        all_Aq2 = [all_Aq2, Aq2];
    end

    % --- Ensemble average ---
    mean_Aq2 = mean(all_Aq2, 2, 'omitnan');

    meanResults{fidx} = mean_Aq2;
    qResults{fidx} = qvals;

    % Save CSV
    writematrix([qvals mean_Aq2], fullfile(folderPath, 'Ensemble_Aq2_um.csv'));
end
% === Log-log plot with aligned q-range ===
figure; hold on;
colors = lines(numFolders);

% --- Compute global q range ---
q_min_common = max(cellfun(@(q) q(3), qResults));  % skip first 2 modes
q_max_common = min(cellfun(@(q) max(q), qResults));

for fidx = 1:numFolders
    q = qResults{fidx};
    Aq2 = meanResults{fidx};

    % Align q range
    mask = (q >= q_min_common) & (q <= q_max_common);
    q_plot = q(mask);
    Aq2_plot = Aq2(mask);

    % --- Plot experimental spectrum ---
    loglog(q_plot, Aq2_plot, 'o-', ...
        'Color', colors(fidx,:), ...
        'LineWidth', 1, ...
        'MarkerSize', 5, ...
        'DisplayName', sprintf('Folder %d', fidx));
end
% === Theory line q^-4===
q_theory = linspace(q_min_common, 10, 300);
C = meanResults{1}(5) * qResults{1}(5)^1;
theory = C ./ (q_theory.^1);

loglog(q_theory, theory, '--k', ...
    'LineWidth', 2, ...
    'DisplayName', 'q^{-1} theory');
% === Axis formatting ===
xlabel('Wavenumber q (\mum^{-1})');
ylabel('<|A_q|^2> (\mum)');
title('Microtubule Bending Spectrum (log-log)');
legend('show');

set(gca, 'XScale', 'log', 'YScale', 'log'); % force log scale
grid on;
box on;
set(gcf, 'Color', 'w');
disp('=== Analysis Complete ===');