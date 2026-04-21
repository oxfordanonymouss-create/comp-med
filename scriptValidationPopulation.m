% Population validation analysis. Loads population_sweep.mat (produced by
% scriptSweep_Population_Parallel.m) and produces the qualitative validation
% figures for the "Population validation" section of the paper.
%
% Six figures are produced (shown on screen, not saved):
%   - Mean +/- SD dAPD90 vs cell concentration, sex-stratified, one panel
%     per cell type.
%   - Boxplots of dAPD90 distributions at low / medium / high cell
%     concentrations, M vs F, one row per cell type.
%   - Per-sex and aggregate prevalence of dAPD90 > T ms, one panel per
%     cell type, T = 30/40 overlaid (central T = 30, bracket {30, 40}).
%   - Last-beat AP traces for all population cells at low / medium / high
%     [ATO], M vs F, one row per cell type.
%   - Last-beat AP traces at Bangladesh-realistic calibrated [ATO]
%     (~0.024 / 0.030 / 0.052 umol/L grid points), M vs F, one row per
%     cell type.
%   - Standalone Mid population AP traces at [ATO] = 0.030 umol/L,
%     M vs F (same data as the Mid / 0.030 panel above).

clear

%% Load sweep outputs
sweepFile = 'population_sweep.mat';
if ~exist(sweepFile, 'file')
    error(['%s not found. Run scriptSweep_Population_Parallel.m first ' ...
        'to produce the sweep outputs.'], sweepFile);
end
load(sweepFile, 'apd90M', 'apd90F', 'concentrations', 'cellTypeNames', ...
    'tracesM', 'tracesF');

[nPop, nConc, nCT] = size(apd90M);
fprintf('Loaded population sweep: %d cells x %d concentrations x %d cell types\n', ...
    nPop, nConc, nCT);

%% Compute dAPD90 per cell, relative to each cell's own baseline (C = 0)
dAPD90M = apd90M - apd90M(:, 1, :);
dAPD90F = apd90F - apd90F(:, 1, :);

% Strip baseline column for plotting on log x-axis
xConc = concentrations(2:end);
nC = length(xConc);

%% Figure 1: Mean +/- SD dAPD90 vs concentration, per cell type
fig1 = figure('Position', [100 100 1200 350]); clf
for ct = 1:nCT
    mM = squeeze(mean(dAPD90M(:, 2:end, ct), 1));
    sM = squeeze(std (dAPD90M(:, 2:end, ct), 0, 1));
    mF = squeeze(mean(dAPD90F(:, 2:end, ct), 1));
    sF = squeeze(std (dAPD90F(:, 2:end, ct), 0, 1));

    subplot(1, nCT, ct);
    fill([xConc, fliplr(xConc)], [mM - sM, fliplr(mM + sM)], 'b', ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on
    fill([xConc, fliplr(xConc)], [mF - sF, fliplr(mF + sF)], 'r', ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(xConc, mM, 'b-o', 'LineWidth', 1.8);
    plot(xConc, mF, 'r-o', 'LineWidth', 1.8);
    set(gca, 'XScale', 'log', 'FontSize', 12);
    xlabel('[ATO] (\mumol/L)');
    ylabel('\DeltaAPD_{90} (ms)');
    title(cellTypeNames{ct});
    if ct == 1
        legend({'Male SD', 'Female SD', 'Male mean', 'Female mean'}, ...
            'Location', 'northwest', 'FontSize', 10);
    end
end
sgtitle('Population mean \DeltaAPD_{90} \pm SD across virtual cells');

%% Figure 2: Boxplots at three representative concentrations, per cell type
% Pick low / medium / high cell concentrations (well below IC50, near IC50,
% well above IC50 for Kr).
IC50_Kr = 0.14;
targetConcs = [0.01, IC50_Kr, 1.0];
boxIdx = zeros(1, length(targetConcs));
for k = 1:length(targetConcs)
    [~, boxIdx(k)] = min(abs(concentrations - targetConcs(k)));
end

fig2 = figure('Position', [100 100 1200 900]); clf
for ct = 1:nCT
    for k = 1:length(targetConcs)
        subplot(nCT, length(targetConcs), (ct-1)*length(targetConcs) + k);
        data  = [dAPD90M(:, boxIdx(k), ct); dAPD90F(:, boxIdx(k), ct)];
        group = [ones(nPop, 1); 2*ones(nPop, 1)];
        boxplot(data, group, 'Labels', {'Male', 'Female'}, ...
            'Colors', 'br', 'Symbol', 'k.');
        ylabel('\DeltaAPD_{90} (ms)');
        title(sprintf('%s, [ATO] = %.3f \\mumol/L', ...
            cellTypeNames{ct}, concentrations(boxIdx(k))), 'FontSize', 10);
        set(gca, 'FontSize', 11);
    end
end
sgtitle('Population \DeltaAPD_{90} distributions at low / medium / high [ATO]');

%% Figure 3: Prevalence of dAPD90 > T, per cell type
% Per-sex curves at central T = 30 ms; aggregate (M+F 50/50) at T in {30,40}
% (Mirams 2014 dQT/dAPD90 scaling 1.35 -> T=30; 1.0 -> T=40).
T_central = 30;
T_sens    = [30, 40];
sensColors = [0 0 0; 0.5 0.5 0.5];
sensStyles = {'-', '--'};

fig3 = figure('Position', [100 100 1200 350]); clf
for ct = 1:nCT
    subplot(1, nCT, ct);
    % Per-sex prevalence at central threshold
    prevM = mean(dAPD90M(:, 2:end, ct) > T_central, 1) * 100;
    prevF = mean(dAPD90F(:, 2:end, ct) > T_central, 1) * 100;
    semilogx(xConc, prevM, 'b-o', 'LineWidth', 1.8); hold on
    semilogx(xConc, prevF, 'r-o', 'LineWidth', 1.8);
    % Aggregate prevalence at three thresholds
    for s = 1:length(T_sens)
        agg = mean([dAPD90M(:, 2:end, ct); dAPD90F(:, 2:end, ct)] > T_sens(s), 1) * 100;
        semilogx(xConc, agg, sensStyles{s}, 'Color', sensColors(s, :), ...
            'LineWidth', 1.2);
    end
    set(gca, 'FontSize', 12);
    xlabel('[ATO] (\mumol/L)');
    ylabel('Prevalence (%)');
    title(cellTypeNames{ct});
    if ct == 1
        legend({sprintf('Male (T=%d)', T_central), ...
                sprintf('Female (T=%d)', T_central), ...
                sprintf('Aggregate (T=%d)', T_sens(1)), ...
                sprintf('Aggregate (T=%d)', T_sens(2))}, ...
            'Location', 'northwest', 'FontSize', 9);
    end
end
sgtitle('Population prevalence of \DeltaAPD_{90} > T (ms)');

%% Figure 4: AP traces at three representative concentrations, per cell type
% Overlay last-beat AP traces for all population cells, M vs F.
fig4 = figure('Position', [100 100 1200 900]); clf
for ct = 1:nCT
    for k = 1:length(targetConcs)
        subplot(nCT, length(targetConcs), (ct-1)*length(targetConcs) + k);
        hold on
        for j = 1:nPop
            trM = tracesM{j, boxIdx(k), ct};
            trF = tracesF{j, boxIdx(k), ct};
            tM = trM.time - trM.time(1);
            tF = trF.time - trF.time(1);
            plot(tM, trM.V, 'Color', [0 0 1 0.25], 'LineWidth', 0.8);
            plot(tF, trF.V, 'Color', [1 0 0 0.25], 'LineWidth', 0.8);
        end
        xlabel('Time (ms)');
        ylabel('V_m (mV)');
        xlim([0 bcl_guess(trM.time)]);
        title(sprintf('%s, [ATO] = %.3f \\mumol/L', ...
            cellTypeNames{ct}, concentrations(boxIdx(k))), 'FontSize', 10);
        set(gca, 'FontSize', 11);
        if ct == 1 && k == 1
            hM = plot(nan, nan, 'b-', 'LineWidth', 1.5);
            hF = plot(nan, nan, 'r-', 'LineWidth', 1.5);
            legend([hM, hF], {'Male', 'Female'}, 'Location', 'northeast', ...
                'FontSize', 9);
        end
    end
end
sgtitle('Population AP traces at low / medium / high [ATO]');

%% Figure 5: AP traces at Bangladesh-realistic calibrated concentrations
% Targets come from scriptConcentrationMapping.m (Mumford-anchor cell
% concentrations at T = 30 ms central threshold, clustered around
% 0.02-0.04 umol/L across cell types). The sweep grid has no entry close
% to 0.04 that is distinct from 0.0302, so the high column is picked at
% the next available grid point (0.052 umol/L), consistent with the
% T = 40 ms sensitivity envelope.
realisticTargets = [0.024, 0.030, 0.052];
realIdx = zeros(1, length(realisticTargets));
for k = 1:length(realisticTargets)
    [~, realIdx(k)] = min(abs(concentrations - realisticTargets(k)));
end

fig5 = figure('Position', [100 100 1200 900]); clf
for ct = 1:nCT
    for k = 1:length(realisticTargets)
        subplot(nCT, length(realisticTargets), (ct-1)*length(realisticTargets) + k);
        hold on
        for j = 1:nPop
            trM = tracesM{j, realIdx(k), ct};
            trF = tracesF{j, realIdx(k), ct};
            tM = trM.time - trM.time(1);
            tF = trF.time - trF.time(1);
            plot(tM, trM.V, 'Color', [0 0 1 0.25], 'LineWidth', 0.8);
            plot(tF, trF.V, 'Color', [1 0 0 0.25], 'LineWidth', 0.8);
        end
        xlabel('Time (ms)');
        ylabel('V_m (mV)');
        xlim([0 bcl_guess(trM.time)]);
        title(sprintf('%s, [ATO] = %.3f \\mumol/L', ...
            cellTypeNames{ct}, concentrations(realIdx(k))), 'FontSize', 10);
        set(gca, 'FontSize', 11);
        if ct == 1 && k == 1
            hM = plot(nan, nan, 'b-', 'LineWidth', 1.5);
            hF = plot(nan, nan, 'r-', 'LineWidth', 1.5);
            legend([hM, hF], {'Male', 'Female'}, 'Location', 'northeast', ...
                'FontSize', 9);
        end
    end
end
sgtitle('Population AP traces at Bangladesh-realistic [ATO] (Mumford-calibrated)');

%% Figure 6: Standalone Mid AP traces at [ATO] = 0.030 umol/L
% Same data as the Mid / 0.030 umol/L panel of Figure 5, shown on its own.
ctMid = find(strcmp(cellTypeNames, 'Mid'), 1);
[~, idx030] = min(abs(concentrations - 0.03));

fig6 = figure('Position', [100 100 700 500]); clf
hold on
for j = 1:nPop
    trM = tracesM{j, idx030, ctMid};
    trF = tracesF{j, idx030, ctMid};
    tM = trM.time - trM.time(1);
    tF = trF.time - trF.time(1);
    plot(tM, trM.V, 'Color', [0 0 1 0.25], 'LineWidth', 0.8);
    plot(tF, trF.V, 'Color', [1 0 0 0.25], 'LineWidth', 0.8);
end
xlabel('Time (ms)', 'FontSize', 14);
ylabel('V_m (mV)', 'FontSize', 14);
xlim([0 bcl_guess(trM.time)]);
set(gca, 'FontSize', 13);
hM = plot(nan, nan, 'b-', 'LineWidth', 1.5);
hF = plot(nan, nan, 'r-', 'LineWidth', 1.5);
legend([hM, hF], {'Male', 'Female'}, 'Location', 'northeast', 'FontSize', 12);
% title(sprintf('Mid population AP traces, [ATO] = %.3f \\mumol/L', ...
%     concentrations(idx030)), 'FontSize', 14);

%% Figure 7: Endo AP traces at the five highest [ATO] on the sweep grid
ctEndo = find(strcmp(cellTypeNames, 'Endo'), 1);
idxEndoHigh = (nConc-4):nConc;

fig7 = figure('Position', [100 100 1800 400]); clf
for k = 1:length(idxEndoHigh)
    subplot(1, length(idxEndoHigh), k);
    hold on
    for j = 1:nPop
        trM = tracesM{j, idxEndoHigh(k), ctEndo};
        trF = tracesF{j, idxEndoHigh(k), ctEndo};
        tM = trM.time - trM.time(1);
        tF = trF.time - trF.time(1);
        plot(tM, trM.V, 'Color', [0 0 1 0.25], 'LineWidth', 0.8);
        plot(tF, trF.V, 'Color', [1 0 0 0.25], 'LineWidth', 0.8);
    end
    xlabel('Time (ms)', 'FontSize', 14);
    ylabel('V_m (mV)', 'FontSize', 14);
    xlim([0 bcl_guess(trM.time)]);
    set(gca, 'FontSize', 13);
    title(sprintf('[ATO] = %.3f \\mumol/L', concentrations(idxEndoHigh(k))), ...
        'FontSize', 13);
    if k == 1
        hM = plot(nan, nan, 'b-', 'LineWidth', 1.5);
        hF = plot(nan, nan, 'r-', 'LineWidth', 1.5);
        legend([hM, hF], {'Male', 'Female'}, 'Location', 'northeast', ...
            'FontSize', 12);
    end
end

%% Console summary
fprintf('\n=== POPULATION VALIDATION SUMMARY ===\n');
for ct = 1:nCT
    fprintf('--- %s ---\n', cellTypeNames{ct});
    mM = squeeze(mean(dAPD90M(:, :, ct), 1));
    mF = squeeze(mean(dAPD90F(:, :, ct), 1));
    fprintf('  Mean dAPD90 at C = %.4f umol/L: Male %5.1f ms,  Female %5.1f ms\n', ...
        concentrations(boxIdx(1)), mM(boxIdx(1)), mF(boxIdx(1)));
    fprintf('  Mean dAPD90 at C = %.4f umol/L: Male %5.1f ms,  Female %5.1f ms\n', ...
        concentrations(boxIdx(2)), mM(boxIdx(2)), mF(boxIdx(2)));
    fprintf('  Mean dAPD90 at C = %.4f umol/L: Male %5.1f ms,  Female %5.1f ms\n', ...
        concentrations(boxIdx(3)), mM(boxIdx(3)), mF(boxIdx(3)));
    qualPass = all(mF(2:end) >= mM(2:end));
    fprintf('  Female mean dAPD90 >= Male at all C > 0: %s\n', ...
        ternary(qualPass, 'PASS', 'FAIL'));
end

%% Helpers
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function L = bcl_guess(t)
    % Guess display window from the last-beat trace duration.
    L = t(end) - t(1);
    if L <= 0 || ~isfinite(L), L = 1000; end
end
