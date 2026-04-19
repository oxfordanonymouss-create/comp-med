% Sensitivity analysis: re-run the single-cell concentration-response with
% perturbed Hill coefficient and IC50 values, and check whether the
% qualitative sex difference (female dAPD90 > male dAPD90) survives
% parameter uncertainty in the Drolet et al. (2004) calibration data.
%
% One-at-a-time perturbations around the baseline (h = 1, IC50_Kr = 0.14,
% IC50_Ks = 1.13):
%   - Hill coefficient h in {0.7, 1.0, 1.5}
%   - IC50_Kr in {0.13, 0.14, 0.15}   (+/- 1 SD; Drolet 0.14 +/- 0.01)
%   - IC50_Ks in {1.07, 1.13, 1.19}   (+/- 1 SD; Drolet 1.13 +/- 0.06)
%
% Parallel version: the (scenario, concentration, cell-type) grid is
% flattened and distributed across workers via parfor.
%% Setting parameters
clear

totalTic = tic;

beats = 100;
showLastN = 2;
ignoreFirst = beats - showLastN;
bcl = 1000;
options = [];

% Baseline IC50 values (umol/L) and Hill coefficient — Drolet et al. 2004
IC50_Kr_base = 0.14;
IC50_Ks_base = 1.13;
h_base = 1;

% IC50 standard deviations from Drolet et al. 2004
SD_Kr = 0.01;
SD_Ks = 0.06;

% Concentration range (driven by baseline IC50 values)
nConcPerChannel = 10;
concKr = logspace(log10(0.01 * IC50_Kr_base), log10(10 * IC50_Kr_base), nConcPerChannel);
concKs = logspace(log10(0.01 * IC50_Ks_base), log10(10 * IC50_Ks_base), nConcPerChannel);
concentrations = [0, unique([concKr, concKs])];
nConc = length(concentrations);

cellTypes = [0, 1, 2]; % endo, epi, mid
cellTypeNames = {'Endocardial', 'Epicardial', 'Midmyocardial'};
nCT = length(cellTypes);

%% Define one-at-a-time scenarios
% Each scenario is (label, h, IC50_Kr, IC50_Ks, group)
% group = 1 (Hill), 2 (IC50_Kr), 3 (IC50_Ks); baseline included in all groups
scenarios = struct('label', {}, 'h', {}, 'IC50_Kr', {}, 'IC50_Ks', {}, 'group', {});

scenarios(end+1) = struct('label', 'h = 0.7',          'h', 0.7, 'IC50_Kr', IC50_Kr_base,        'IC50_Ks', IC50_Ks_base,        'group', 1);
scenarios(end+1) = struct('label', '(base)',           'h', 1.0, 'IC50_Kr', IC50_Kr_base,        'IC50_Ks', IC50_Ks_base,        'group', 1);
scenarios(end+1) = struct('label', 'h = 1.5',          'h', 1.5, 'IC50_Kr', IC50_Kr_base,        'IC50_Ks', IC50_Ks_base,        'group', 1);
scenarios(end+1) = struct('label', 'IC50_{Kr} - 1SD',  'h', h_base, 'IC50_Kr', IC50_Kr_base-SD_Kr, 'IC50_Ks', IC50_Ks_base,        'group', 2);
scenarios(end+1) = struct('label', 'IC50_{Kr} + 1SD',  'h', h_base, 'IC50_Kr', IC50_Kr_base+SD_Kr, 'IC50_Ks', IC50_Ks_base,        'group', 2);
scenarios(end+1) = struct('label', 'IC50_{Ks} - 1SD',  'h', h_base, 'IC50_Kr', IC50_Kr_base,        'IC50_Ks', IC50_Ks_base-SD_Ks, 'group', 3);
scenarios(end+1) = struct('label', 'IC50_{Ks} + 1SD',  'h', h_base, 'IC50_Kr', IC50_Kr_base,        'IC50_Ks', IC50_Ks_base+SD_Ks, 'group', 3);

nScen = length(scenarios);

%% Start parallel pool
if isempty(gcp('nocreate'))
    parpool('local', 4);
end

%% Sweep (scenario, concentration, cell-type) — flattened grid, parallelized
tic
nTotal = nScen * nConc * nCT;
[sIdx, iIdx, ctIdx] = ind2sub([nScen, nConc, nCT], 1:nTotal);

apd90M_flat = zeros(nTotal, 1);
apd90F_flat = zeros(nTotal, 1);
traceM_flat = cell(nTotal, 1);
traceF_flat = cell(nTotal, 1);

% Pull scenario fields into plain arrays for parfor slicing
scen_h = [scenarios.h];
scen_IC50_Kr = [scenarios.IC50_Kr];
scen_IC50_Ks = [scenarios.IC50_Ks];

parfor idx = 1:nTotal
    s = sIdx(idx);
    i = iIdx(idx);
    ct = ctIdx(idx);
    X = concentrations(i);

    h_s        = scen_h(s);
    IC50_Kr_s  = scen_IC50_Kr(s);
    IC50_Ks_s  = scen_IC50_Ks(s);

    paramM = struct();
    paramM.bcl = bcl;
    paramM.model = @model_ToRORd_Land_Male;
    paramM.cellType = cellTypes(ct);
    paramM.IKr_Multiplier = 1 / (1 + (X / IC50_Kr_s)^h_s);
    paramM.IKs_Multiplier = 1 / (1 + (X / IC50_Ks_s)^h_s);

    paramF = struct();
    paramF.bcl = bcl;
    paramF.model = @model_ToRORd_Land_Female;
    paramF.cellType = cellTypes(ct);
    paramF.IKr_Multiplier = paramM.IKr_Multiplier;
    paramF.IKs_Multiplier = paramM.IKs_Multiplier;

    X0 = getStartingState('m_endo');

    [tM, XM] = modelRunner(X0, options, paramM, beats, ignoreFirst);
    cM = getCurrentsStructure(tM, XM, beats, paramM, 0);
    apd90M_flat(idx) = getAPD(cM.time, cM.V, 0.9, bcl);
    traceM_flat{idx} = struct('time', cM.time, 'V', cM.V);

    [tF, XF] = modelRunner(X0, options, paramF, beats, ignoreFirst);
    cF = getCurrentsStructure(tF, XF, beats, paramF, 0);
    apd90F_flat(idx) = getAPD(cF.time, cF.V, 0.9, bcl);
    traceF_flat{idx} = struct('time', cF.time, 'V', cF.V);

    fprintf('Scenario %d/%d | %s | [ATO] = %.4f uM | h = %.2f | IC50_Kr = %.3f | IC50_Ks = %.3f | M: %.1f ms | F: %.1f ms\n', ...
        s, nScen, cellTypeNames{ct}, X, h_s, IC50_Kr_s, IC50_Ks_s, apd90M_flat(idx), apd90F_flat(idx));
end

apd90M = reshape(apd90M_flat, [nScen, nConc, nCT]);
apd90F = reshape(apd90F_flat, [nScen, nConc, nCT]);
traceM = reshape(traceM_flat, [nScen, nConc, nCT]);
traceF = reshape(traceF_flat, [nScen, nConc, nCT]);

fprintf('Simulation time: %.1f seconds\n', toc);

%% dAPD90 relative to each scenario's own baseline (concentration 0)
dAPD90M = apd90M - apd90M(:, 1, :);
dAPD90F = apd90F - apd90F(:, 1, :);

%% Plot — 3x3 grid: rows = cell types, cols = parameter groups
groupTitles = {'Hill coefficient', 'IC_{50} (I_{Kr})', 'IC_{50} (I_{Ks})'};
lineStyles = {'--', '-', ':'}; % low, baseline, high (within each group)
x = concentrations(2:end);
baseIdx = find(strcmp({scenarios.label}, 'h = 1.0 (base)'));

figure(1); clf
tl = tiledlayout(nCT, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for ct = 1:nCT
    for g = 1:3
        nexttile;
        scenIdx = find([scenarios.group] == g);

        % Order within group: low, baseline, high
        if g == 1
            order = scenIdx; % already 0.7, 1.0, 1.5
        else
            order = [scenIdx(1), baseIdx, scenIdx(2)];
        end

        hold on
        legHandles = []; legLabels = {};
        for k = 1:length(order)
            s = order(k);
            ls = lineStyles{k};
            hM = semilogx(x, squeeze(dAPD90M(s, 2:end, ct)), ['b', ls, 'o'], 'LineWidth', 1.5, 'MarkerSize', 4);
            hF = semilogx(x, squeeze(dAPD90F(s, 2:end, ct)), ['r', ls, 'o'], 'LineWidth', 1.5, 'MarkerSize', 4);
            legHandles = [legHandles, hM, hF];
            legLabels = [legLabels, {sprintf('M %s', scenarios(s).label), sprintf('F %s', scenarios(s).label)}];
        end
        set(gca, 'XScale', 'log', 'FontSize', 11);
        if ct == nCT
            xlabel('[ATO] (\mumol/L)', 'FontSize', 12);
        end
        if g == 1
            ylabel({cellTypeNames{ct}, '\DeltaAPD_{90} (ms)'}, 'FontSize', 12);
        end
        if ct == 1
            title(groupTitles{g}, 'FontSize', 13);
        end
        if ct == 1 && g == 1
            legend(legHandles, legLabels, 'Location', 'northwest', 'FontSize', 8);
        end
        grid on
    end
end
title(tl, 'Sensitivity analysis — Hill coefficient and IC_{50} perturbations', 'FontSize', 15);

%% AP traces — one figure per (cell type, parameter group)
% 2 rows (Male / Female) x 3 cols (low / baseline / high parameter value)
% Each tile overlays AP traces at all non-zero concentrations, colored by [ATO]
colors = turbo(nConc);
sexNames = {'Male', 'Female'};
figNum = 1;
for ct = 1:nCT
    for g = 1:3
        figNum = figNum + 1;
        figure(figNum); clf
        scenIdx = find([scenarios.group] == g);
        if g == 1
            order = scenIdx; % 0.7, 1.0, 1.5
        else
            order = [scenIdx(1), baseIdx, scenIdx(2)];
        end

        tl_ap = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
        for sex = 1:2
            for k = 1:length(order)
                s = order(k);
                ax = nexttile;
                hold(ax, 'on');
                for i = 2:nConc % skip baseline (X = 0)
                    if sex == 1
                        plot(ax, traceM{s, i, ct}.time, traceM{s, i, ct}.V, ...
                            'Color', colors(i, :), 'LineWidth', 1);
                    else
                        plot(ax, traceF{s, i, ct}.time, traceF{s, i, ct}.V, ...
                            'Color', colors(i, :), 'LineWidth', 1);
                    end
                end
                set(ax, 'FontSize', 11);
                xlim(ax, [0 showLastN * bcl]); ylim(ax, [-100 60]);
                if sex == 1
                    title(ax, scenarios(s).label, 'FontSize', 12);
                end
                if k == 1
                    ylabel(ax, {sexNames{sex}, 'V (mV)'}, 'FontSize', 12);
                end
                if sex == 2
                    xlabel(ax, 'Time (ms)', 'FontSize', 12);
                end
                caxis(ax, log10([concentrations(2), concentrations(end)]));
            end
        end
        title(tl_ap, sprintf('%s — %s', cellTypeNames{ct}, groupTitles{g}), ...
            'FontSize', 14);
        colormap(turbo);
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.Label.String = '[ATO] (\mumol/L)';
        cb.Label.FontSize = 12;
        cb.Ticks = log10([0.01, 0.1, 1, 10]);
        cb.TickLabels = {'0.01', '0.1', '1', '10'};
    end
end

fprintf('Total elapsed time: %.1f seconds\n', toc(totalTic));

%% Helper functions
function apd = getAPD(time, V, level, bcl)
    % Measure APD on the last beat only
    tStart = time(end) - bcl;
    lastBeat = time >= tStart;
    time = time(lastBeat);
    V = V(lastBeat);

    baseline  = V(end);
    threshold = baseline + (1 - level) * (max(V) - baseline);
    above = V > threshold;
    startIdx = find(above, 1, 'first');
    endIdx = find(above & (1:length(above))' >= startIdx, 1, 'last');
    apd = time(endIdx) - time(startIdx);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
