% Virtual ablation: at three reference arsenic concentrations corresponding
% to the Mumford low / medium / high exposure groups (mapped to cell via
% scriptConcentrationMapping.m at T = 30 ms central threshold; 0.02, 0.03,
% 0.04 umol/L), measure how much
% of the female-minus-male dAPD90 gap remains when each sex-difference
% factor in the model is swapped from male to female, one at a time. The
% swap that drives the residual gap toward zero identifies the ion-channel
% sex difference responsible for the greater female response to acute
% IKr/IKs pore block.
%
% Swap factors are the cell-type-specific male/female ratios of the
% hard-coded conductance scalings in model_ToRORd_Land_Male.m and
% model_ToRORd_Land_Female.m.
%
% The residual gap is reported in milliseconds rather than as a
% percentage of the baseline gap: at low concentrations the baseline
% gap is small, so a percentage divides by a near-zero denominator and
% amplifies solver/detector noise into meaningless multi-hundred-percent
% swings. A millisecond-scale residual is directly interpretable and
% comparable across concentrations.
%
% Output: one figure with three bar-chart panels (one per concentration)
% and one console table.
clear

beats = 100;
showLastN = 2;
ignoreFirst = beats - showLastN;
bcl = 1000;
options = [];

% Reference concentrations: Mumford low / medium / high exposure groups
% (cell-side values from scriptConcentrationMapping.m)
refConcs = [0.02, 0.03, 0.04];  % umol/L
IC50_Kr = 0.14;
IC50_Ks = 1.13;
h = 1;
nRef = length(refConcs);

cellTypes     = [0, 1, 2];
cellTypeNames = {'Endo', 'Epi', 'Mid'};
nCT = length(cellTypes);

% Per-cell-type sex-difference swap factors = (male endo G) / (female endo G)
ablatedChannels = {'IKr_Multiplier', 'IKs_Multiplier', 'IK1_Multiplier', 'INaCa_Multiplier', 'IpCa_Multiplier'};
ablatedLabels   = {'I_{Kr}',         'I_{Ks}',         'I_{K1}',         'I_{NaCa}',         'I_{pCa}'};
swapFactor = [ ...
    1.00/0.79,  1.09/0.875, 0.80/0.632;  ... % IKr
    1.00/0.83,  1.04/0.87,  1.00/0.83;   ... % IKs
    1.00/0.86,  0.98/0.74,  1.30/1.18;   ... % IK1
    1.00/1.15,  1.00/1.15,  1.40/1.40;   ... % INaCa (mid: null)
    1.00/1.60,  0.88/1.60,  2.50/4.00    ... % IpCa
];
nAbl = length(ablatedLabels);

% Scenarios: male baseline, female baseline, one female-with-male-X per channel
scenarios = [{'Male', 'M', 0}; {'Female', 'F', 0}];
for k = 1:nAbl
    scenarios(end+1, :) = {sprintf('F + M %s', ablatedLabels{k}), 'F', k}; %#ok<SAGROW>
end
nScen = size(scenarios, 1);

% Concentrations: 0 (for each scenario's own baseline) + the three refs
concentrations = [0, refConcs];
nConc = length(concentrations);

%% Run simulations -- serial, 84 of them
tic
apd90 = zeros(nScen, nConc, nCT);
for ct = 1:nCT
    for s = 1:nScen
        for i = 1:nConc
            X = concentrations(i);

            param = struct();
            param.bcl = bcl;
            param.cellType = cellTypes(ct);
            if strcmp(scenarios{s, 2}, 'M')
                param.model = @model_ToRORd_Land_Male;
            else
                param.model = @model_ToRORd_Land_Female;
            end

            fracKr = 1 / (1 + (X / IC50_Kr)^h);
            fracKs = 1 / (1 + (X / IC50_Ks)^h);
            param.IKr_Multiplier = fracKr;
            param.IKs_Multiplier = fracKs;

            swapIdx = scenarios{s, 3};
            if swapIdx > 0
                field = ablatedChannels{swapIdx};
                if isfield(param, field)
                    current = param.(field);
                else
                    current = 1;
                end
                param.(field) = current * swapFactor(swapIdx, ct);
            end

            X0 = getStartingState('m_endo');
            [tm, Xm] = modelRunner(X0, options, param, beats, ignoreFirst);
            cm = getCurrentsStructure(tm, Xm, beats, param, 0);
            apd90(s, i, ct) = getAPD(cm.time, cm.V, 0.9, bcl);

            fprintf('%-4s | %-20s | [ATO] = %.4f uM | APD90 = %.1f ms\n', ...
                cellTypeNames{ct}, scenarios{s, 1}, X, apd90(s, i, ct));
        end
    end
end
fprintf('Simulation time: %.1f seconds\n', toc);

%% dAPD90 = drug - baseline, per scenario per ref concentration per cell type
% dAPD90 shape: [nScen, nRef, nCT]
dAPD90 = apd90(:, 2:end, :) - apd90(:, 1, :);

% Residual F-M dAPD90 gap (ms) after each single-channel swap, and the
% baseline F-M gap for reference. Shapes: [nAbl, nRef, nCT] and [nRef, nCT].
residualGap = zeros(nAbl, nRef, nCT);
baselineGap = zeros(nRef, nCT);
for ct = 1:nCT
    for r = 1:nRef
        dM = dAPD90(1, r, ct);
        baselineGap(r, ct) = dAPD90(2, r, ct) - dM;
        for k = 1:nAbl
            residualGap(k, r, ct) = dAPD90(2 + k, r, ct) - dM;
        end
    end
end

%% Figure -- three tiled bar panels, one per concentration
% In each panel: leftmost group = baseline F-M gap, then one group per
% single-channel swap showing the residual F-M gap after that swap.
%   bar close to 0            -> swap closes the gap
%   bar equal to baseline     -> swap has no effect
%   bar larger than baseline  -> swap widens the gap
%   bar opposite-signed       -> swap overshoots past the male
figure(1); clf
tl = tiledlayout(1, nRef, 'Padding', 'compact', 'TileSpacing', 'compact');

% Shared y-limits for direct visual comparison across panels
allBars = [baselineGap(:); residualGap(:)];
yMargin = 0.1 * max(abs(allBars));
yMin = min(allBars) - yMargin;
yMax = max(allBars) + yMargin;

barLabels = [{'Baseline F-M'}, ablatedLabels];
for r = 1:nRef
    nexttile;
    barData = [baselineGap(r, :); squeeze(residualGap(:, r, :))]; % (1 + nAbl) x nCT
    bar(barData);
    set(gca, 'XTickLabel', barLabels, 'XTickLabelRotation', 30, 'FontSize', 10);
    ylim([yMin, yMax]);
    yline(0, 'k-', 'LineWidth', 1);
    title(sprintf('[ATO] = %.2f \\mumol/L', refConcs(r)), 'FontSize', 12);
    if r == 1
        ylabel('Residual F - M \DeltaAPD_{90} (ms)', 'FontSize', 12);
        legend(cellTypeNames, 'Location', 'best', 'FontSize', 10);
    end
    grid on
end
title(tl, 'Virtual ablation -- residual F-M \DeltaAPD_{90} gap after each single-channel swap', ...
    'FontSize', 13);

%% Table
fprintf('\n=== VIRTUAL ABLATION -- residual F-M dAPD90 gap (ms) per channel swap ===\n');
for r = 1:nRef
    fprintf('\n[ATO] = %.2f umol/L\n', refConcs(r));
    fprintf('%-14s %10s %10s %10s\n', 'Row', 'Endo (ms)', 'Epi (ms)', 'Mid (ms)');
    fprintf('%-14s %+10.1f %+10.1f %+10.1f\n', 'Baseline F-M', ...
        baselineGap(r, 1), baselineGap(r, 2), baselineGap(r, 3));
    for k = 1:nAbl
        fprintf('%-14s %+10.1f %+10.1f %+10.1f\n', ...
            sprintf('F+M %s', ablatedLabels{k}), ...
            residualGap(k, r, 1), residualGap(k, r, 2), residualGap(k, r, 3));
    end
end
fprintf('\nRead as: bar close to 0            = swap closes the gap;\n');
fprintf('         bar equal to "Baseline F-M" = swap has no effect;\n');
fprintf('         bar larger than baseline    = swap widens the gap;\n');
fprintf('         bar opposite-signed         = swap overshoots past the male.\n');

%% Helper
function apd = getAPD(time, V, level, bcl)
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
