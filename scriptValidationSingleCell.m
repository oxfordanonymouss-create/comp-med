% Single-cell validation analysis. Loads singlecell_sweep.mat (produced by
% scriptSweep_SingleCell_Parallel.m) and produces the qualitative validation
% figures for the "Single-cell validation" section of the paper.
%
% Figures produced (shown on screen, not saved):
%   - dAPD90 vs [ATO], Male vs Female, one figure per cell type.
%   - AP traces across the concentration range, Male / Female tiles, one
%     figure per cell type.

clear

%% Load sweep outputs
sweepFile = 'singlecell_sweep.mat';
if ~exist(sweepFile, 'file')
    error(['%s not found. Run scriptSweep_SingleCell_Parallel.m first ' ...
        'to produce the sweep outputs.'], sweepFile);
end
load(sweepFile, 'apd90M', 'apd90F', 'concentrations', 'cellTypeNames', ...
    'traceM', 'traceF');

[nConc, nCT] = size(apd90M);
fprintf('Loaded single-cell sweep: %d concentrations x %d cell types\n', ...
    nConc, nCT);

% Beat window used when plotting traces (must match the sweep script)
bcl = 1000;
showLastN = 2;

%% Compute dAPD90 relative to baseline (first row = concentration 0)
dAPD90M = apd90M - apd90M(1, :);
dAPD90F = apd90F - apd90F(1, :);

%% Figure set 1: dAPD90 vs concentration — one figure per cell type
for ct = 1:nCT
    fig = figure(ct); clf
    semilogx(concentrations(2:end), dAPD90M(2:end, ct), 'b-o', 'LineWidth', 2); hold on
    semilogx(concentrations(2:end), dAPD90F(2:end, ct), 'r-o', 'LineWidth', 2);
    set(gca, 'FontSize', 16);
    xlabel('[ATO] (\mumol/L)', 'FontSize', 18);
    ylabel('\DeltaAPD_{90} (ms)', 'FontSize', 18);
    title(cellTypeNames{ct}, 'FontSize', 18);
    legend('Male', 'Female', 'Location', 'northwest', 'FontSize', 16);
end

%% Figure set 2: AP traces — one figure per cell type, two tiles (male / female)
colors = turbo(nConc);
for ct = 1:nCT
    fig = figure(nCT + ct); clf
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    % Male
    ax1 = nexttile;
    hold(ax1, 'on');
    for i = 1:nConc
        plot(ax1, traceM{i, ct}.time, traceM{i, ct}.V, 'Color', colors(i, :), 'LineWidth', 1);
    end
    set(ax1, 'FontSize', 16);
    xlabel(ax1, 'Time (ms)', 'FontSize', 18); ylabel(ax1, 'V (mV)', 'FontSize', 18);
    title(ax1, 'Male', 'FontSize', 18);
    xlim(ax1, [0 showLastN * bcl]); ylim(ax1, [-100 60]);
    % Female
    ax2 = nexttile;
    hold(ax2, 'on');
    for i = 1:nConc
        plot(ax2, traceF{i, ct}.time, traceF{i, ct}.V, 'Color', colors(i, :), 'LineWidth', 1);
    end
    set(ax2, 'FontSize', 16);
    xlabel(ax2, 'Time (ms)', 'FontSize', 18); ylabel(ax2, 'V (mV)', 'FontSize', 18);
    title(ax2, 'Female', 'FontSize', 18);
    xlim(ax2, [0 showLastN * bcl]); ylim(ax2, [-100 60]);
    % Shared colorbar outside the tiles so panel widths stay equal
    colormap(turbo);
    caxis(ax1, log10([concentrations(2), concentrations(end)]));
    caxis(ax2, log10([concentrations(2), concentrations(end)]));
    cb = colorbar(ax2);
    cb.Layout.Tile = 'east';
    cb.Label.String = '[ATO] (\mumol/L)';
    cb.Label.FontSize = 18;
    cb.FontSize = 16;
    cb.Ticks = log10([0.01, 0.1, 1, 10]);
    cb.TickLabels = {'0.01', '0.1', '1', '10'};
    title(t, cellTypeNames{ct}, 'FontSize', 18);
end
