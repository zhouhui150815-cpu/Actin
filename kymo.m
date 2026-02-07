%% ===============================
% Kymograph 
% 1. Read CSV file of all curvatures
data = readmatrix('control.csv'); 
% CSV 结构：
% First column: cell pair
% First row: time scale
% Others: RMS curvature
% 2. Read data：
order = data(2:end,1);    % cell pair
time  = data(1,2:end);    % time
Iraw  = data(2:end,2:end);% RMS curvature
% 3. Settings
rowHeight    = 1;      % row height
colWidth     = 40;      % column width
jetFraction  = 0.90;   % jet colormap 90%
numColors    = 256;    % colormap resolution
smoothWindow = 3;      % Smoothing along time
useInterp    = false;  % true=Interp smoothing, false=moving mean
cLow         = 0.02;   % color value min
cHigh        = 0.20;   % color value max

% 4. Time direction：smoothing/interp
if useInterp
    tFine = linspace(min(time), max(time), numel(time)*2); 
    I = interp1(time, Iraw.', tFine, 'pchip').'; 
    timePlot = tFine;
else
    I = movmean(Iraw, smoothWindow, 2); % move mean
    timePlot = time;
end
% 5. column/row dimensions
I_exp = repelem(I, rowHeight, colWidth);
order_exp = repelem(order, rowHeight);
time_exp  = repelem(timePlot, colWidth);
% 6. jet colormap (90%)
fullJet = jet(numColors);
maxIdx = round(jetFraction*numColors);
truncatedJet = fullJet(1:maxIdx,:);
% 7. plot kymograph
figure;
imagesc(time_exp, order_exp, I_exp);
% top-down cell pair
set(gca,'YDir','reverse');
xlabel('Time (s)');
ylabel('Order');
title('Kymograph');
colormap(truncatedJet);
% 8. color limits
clim([cLow cHigh]); 
% 9. Colorbar setting
cb = colorbar;
cb.Label.String = 'Intensity';
cb.FontSize = 11;
cb.Ticks = 0.02:0.02:0.18; 
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);
% 10. trick setting
xticks(timeplot);
yticks(order);     
% 11. display
axis tight;
set(gca, 'XTick', [], 'XTickLabel', timePlot, 'FontSize', 12, 'LineWidth', 1);