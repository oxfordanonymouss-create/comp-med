% Concentration mapping: derive water -> cardiomyocyte free arsenite mapping
% by matching the model output against two clinical anchors:
%
%   - Chen et al. 2013: female:male slope ratio of dQTc/dWater = 4.3, at the
%     Bangladesh mean water arsenic of 95 ug/L. Matched against the model's
%     local slope ratio (dAPD_F/dC) / (dAPD_M/dC).
%
%   - Mumford et al. 2007: aggregate QTc-prolongation prevalence of
%     3.9 / 11.1 / 20.6 % across exposure groups (<=21, 100-300, 430-690 ug/L
%     water; midpoints 11, 200, 560). Matched against simulated aggregate
%     prevalence of dAPD90 > T ms. dQTc threshold is 450-410=40 ms; APD-space
%     threshold uses Mirams et al. 2014 scaling dQT ~= 1.35*dAPD90, giving
%     T = 30 central, with T in {30,40} as a sensitivity bracket spanning
%     scaling factors 1.35 to 1.0.
%
% For each cell type (endo, epi, mid) the script derives 4 (water, cell)
% mapping points (1 Chen + 3 Mumford) and fits a linear-through-origin bridge
% C_cell = k * C_water. The four points falling on a single line tests the
% linearity assumption empirically rather than assuming it.
%

clear

%% Constants (declared, not loaded -- see notes in script)
IC50_Kr           = 0.14;      % umol/L, Drolet et al. 2004
chenRatio         = 4.3;       % Chen 2013 female:male slope ratio
chenRatioCIlow    = 2;         % rough envelope from Chen's CIs (women 0.3-5.7,
chenRatioCIhigh   = 10;        % men ~0.7 non-significant)
chenWaterMean     = 95;        % ug/L, Bangladesh mean

mumfordWater      = [11, 200, 560]; % ug/L, midpoints of Mumford groups
mumfordPrev       = [3.9, 11.1, 20.6]; % %, aggregate prolongation prevalence

T_central         = 30;        % ms, central dAPD90 threshold: clinical dQTc=40 ms
                               % (Mumford QTc>=450 minus 410 baseline) mapped to
                               % APD space via Mirams 2014 scaling dQT~=1.35*dAPD
T_sens            = [30, 40];  % sensitivity bracket: scaling 1.35 (T=30) to 1 (T=40)

% PK ceiling per ug/L water (high-exposure 570 ug/L -> ~0.22 umol/L total
% blood; scale linearly to other water values for a per-point ceiling).
pkCeilingPerWater = 0.22 / 570;

%% Load sweep outputs
scFile  = 'singlecell_sweep.mat';
popFile = 'population_sweep.mat';
if ~exist(scFile, 'file')
    error('%s not found. Run scriptSweep_SingleCell_Parallel.m first.', scFile);
end
if ~exist(popFile, 'file')
    error('%s not found. Run scriptSweep_Population_Parallel.m first.', popFile);
end

sc  = load(scFile,  'apd90M', 'apd90F', 'concentrations', 'cellTypeNames');
pop = load(popFile, 'apd90M', 'apd90F', 'concentrations', 'cellTypeNames');

nCT       = length(sc.cellTypeNames);
ctNames   = sc.cellTypeNames;

% Pre-allocate result storage
xrefChen        = nan(nCT, 1);
xrefChenCIlow   = nan(nCT, 1);
xrefChenCIhigh  = nan(nCT, 1);
mumfordXref     = nan(nCT, length(mumfordPrev), length(T_sens));
kFit            = nan(nCT, 1);
ratioCurves     = cell(nCT, 1);  % stored for plotting
prevalenceCurves = cell(nCT, 1); % stored for plotting

%% ========================================================================
%% PART A -- Chen mapping point (single-cell slope-ratio match)
%% ========================================================================
fprintf('\n=== PART A: Chen calibration (single-cell slope-ratio match) ===\n');

% Single-cell concentration grid (drop baseline C = 0)
sc_C = sc.concentrations(2:end);
log10C = log10(sc_C);

% Fine evaluation grid for derivative computation
fineGrid_log10C = linspace(min(log10C), max(log10C), 400);
fineGrid_C      = 10.^fineGrid_log10C;

for ct = 1:nCT
    % Per-cell-type APD90 at each concentration
    apd_M = sc.apd90M(:, ct);
    apd_F = sc.apd90F(:, ct);

    % dAPD90 relative to baseline (drop baseline row)
    dM = apd_M(2:end) - apd_M(1);
    dF = apd_F(2:end) - apd_F(1);

    % PCHIP interpolation on log10(C)
    pp_dM = pchip(log10C, dM);
    pp_dF = pchip(log10C, dF);
    dM_fine = ppval(pp_dM, fineGrid_log10C);
    dF_fine = ppval(pp_dF, fineGrid_log10C);

    % Derivatives w.r.t. C (not log10(C)) via finite differences on linear C
    deriv_M = gradient(dM_fine, fineGrid_C);
    deriv_F = gradient(dF_fine, fineGrid_C);

    % R(C) = (dF/dC) / (dM/dC), with masking
    R = deriv_F ./ deriv_M;
    noiseFloor = 0.01 * max(abs(deriv_M));
    mask = (dF_fine > dM_fine) & (abs(deriv_M) > noiseFloor) & (deriv_M > 0);

    % Restrict crossing search to masked region
    R_masked = R; R_masked(~mask) = NaN;
    ratioCurves{ct} = struct('C', fineGrid_C, 'R', R, 'mask', mask, ...
        'dM', dM_fine, 'dF', dF_fine);

    % Find Xref where R is closest to chenRatio (within mask)
    if any(~isnan(R_masked))
        [~, idxBest] = min(abs(R_masked - chenRatio));
        xrefChen(ct) = fineGrid_C(idxBest);
        % CI band: where R is inside chenRatioCIlow..chenRatioCIhigh
        inBand = (R_masked >= chenRatioCIlow) & (R_masked <= chenRatioCIhigh);
        if any(inBand)
            xrefChenCIlow(ct)  = min(fineGrid_C(inBand));
            xrefChenCIhigh(ct) = max(fineGrid_C(inBand));
        end
    end

    fprintf('--- %s ---\n', ctNames{ct});
    if ~isnan(xrefChen(ct))
        ceilingHere = pkCeilingPerWater * chenWaterMean;
        fprintf('  Xref (R = 4.3): %.4f umol/L\n', xrefChen(ct));
        fprintf('  Xref CI band  : %.4f -- %.4f umol/L (R in [%.1f,%.1f])\n', ...
            xrefChenCIlow(ct), xrefChenCIhigh(ct), chenRatioCIlow, chenRatioCIhigh);
        fprintf('  Xref / IC50_Kr: %.3f  (low-dose regime if << 1)\n', ...
            xrefChen(ct) / IC50_Kr);
        fprintf('  PK ceiling at 95 ug/L water: %.4f umol/L  (Xref %s ceiling)\n', ...
            ceilingHere, ternary(xrefChen(ct) <= ceilingHere, 'PASSES', 'EXCEEDS'));
    else
        fprintf('  No valid Xref crossing (R(C) does not reach 4.3 in masked region).\n');
    end
end

%% ========================================================================
%% PART B -- Mumford mapping points (population prevalence match)
%% ========================================================================
fprintf('\n=== PART B: Mumford calibration (population prevalence match) ===\n');

pop_C = pop.concentrations;
pop_log10C = log10(max(pop_C, eps));   % protect log of 0 baseline (kept for index 1)

% Fine evaluation grid in log-space, but skip the baseline column for prevalence
fineGrid_log10C_pop = linspace(log10(pop_C(2)), log10(pop_C(end)), 400);
fineGrid_C_pop      = 10.^fineGrid_log10C_pop;

[nPop, ~, ~] = size(pop.apd90M);

for ct = 1:nCT
    fprintf('--- %s ---\n', ctNames{ct});

    % Per-cell baseline (C = 0) and dAPD90 across the rest of the grid
    base_M = pop.apd90M(:, 1, ct);  % nPop x 1
    base_F = pop.apd90F(:, 1, ct);
    dM_pop = squeeze(pop.apd90M(:, 2:end, ct)) - base_M;  % nPop x (nConc-1)
    dF_pop = squeeze(pop.apd90F(:, 2:end, ct)) - base_F;

    % Interpolate each cell's dAPD90 onto the fine concentration grid
    dM_fine_pop = zeros(nPop, length(fineGrid_C_pop));
    dF_fine_pop = zeros(nPop, length(fineGrid_C_pop));
    for j = 1:nPop
        dM_fine_pop(j, :) = pchip(log10(pop_C(2:end)), dM_pop(j, :), fineGrid_log10C_pop);
        dF_fine_pop(j, :) = pchip(log10(pop_C(2:end)), dF_pop(j, :), fineGrid_log10C_pop);
    end

    % Aggregate prevalence at each T, evaluated on the fine grid
    aggData = [dM_fine_pop; dF_fine_pop];  % (2*nPop) x nFine
    prevAgg = nan(length(T_sens), length(fineGrid_C_pop));
    for s = 1:length(T_sens)
        prevAgg(s, :) = mean(aggData > T_sens(s), 1) * 100;
    end
    prevalenceCurves{ct} = struct('C', fineGrid_C_pop, 'prevAgg', prevAgg, ...
        'T', T_sens);

    % For each (T, Mumford prevalence target) find closest matching C
    for s = 1:length(T_sens)
        for g = 1:length(mumfordPrev)
            [diffMin, idxBest] = min(abs(prevAgg(s, :) - mumfordPrev(g)));
            if isfinite(diffMin)
                mumfordXref(ct, g, s) = fineGrid_C_pop(idxBest);
            end
        end
    end

    % Print results
    for s = 1:length(T_sens)
        fprintf('  T = %d ms:\n', T_sens(s));
        for g = 1:length(mumfordPrev)
            ceilingHere = pkCeilingPerWater * mumfordWater(g);
            xx = mumfordXref(ct, g, s);
            if isnan(xx)
                fprintf('    Group %d (water=%d ug/L, target prev=%.1f%%): no match\n', ...
                    g, mumfordWater(g), mumfordPrev(g));
            else
                fprintf('    Group %d (water=%d ug/L, target prev=%.1f%%): cell=%.4f umol/L  (ceiling %.4f, %s)\n', ...
                    g, mumfordWater(g), mumfordPrev(g), xx, ceilingHere, ...
                    ternary(xx <= ceilingHere, 'PASSES', 'EXCEEDS'));
            end
        end
    end
end

%% ========================================================================
%% PART C -- Joint analysis: linear bridge fit per cell type
%% ========================================================================
fprintf('\n=== PART C: Linear-bridge fit per cell type (T = %d central) ===\n', T_central);

s_central = find(T_sens == T_central);

for ct = 1:nCT
    % Assemble 4 (water, cell) points: 1 Chen + 3 Mumford at central T
    waters = [chenWaterMean, mumfordWater];  % 1x4
    cells  = [xrefChen(ct), reshape(mumfordXref(ct, :, s_central), 1, [])];  % 1x4

    valid = ~isnan(cells);
    if sum(valid) < 2
        fprintf('--- %s --- insufficient valid points for linear fit (%d valid).\n', ...
            ctNames{ct}, sum(valid));
        continue
    end

    % Linear fit through origin (least squares):  cells = k * waters
    w = waters(valid)';
    c = cells(valid)';
    k = (w' * c) / (w' * w);
    kFit(ct) = k;

    fprintf('--- %s ---\n', ctNames{ct});
    fprintf('  k_%s = %.6f umol/L per ug/L water\n', lower(ctNames{ct}), k);
    fprintf('  Equivalent mass-mass ratio (As atomic mass 75 g/mol): %.4f\n', ...
        k * 75);
    fprintf('  Points used (%d/4 valid):\n', sum(valid));
    for q = 1:4
        if valid(q)
            fprintf('    water=%4d ug/L  ->  cell=%.4f umol/L  (predicted %.4f, residual %+.4f)\n', ...
                waters(q), cells(q), k * waters(q), cells(q) - k * waters(q));
        else
            fprintf('    water=%4d ug/L  ->  cell=NaN  (skipped)\n', waters(q));
        end
    end
end

%% ========================================================================
%% Figures
%% ========================================================================

%% Figure 1: R(C) per cell type with Chen crossings + envelope
fig1 = figure('Position', [100 100 1200 350]); clf
for ct = 1:nCT
    subplot(1, nCT, ct);
    rc = ratioCurves{ct};
    semilogx(rc.C, rc.R, 'k-', 'LineWidth', 1.5); hold on
    % Mask shading
    invalid = ~rc.mask;
    yl = [-2 15];
    fill([rc.C, fliplr(rc.C)], [yl(1)*ones(size(rc.C)) .* invalid - 1e3*(~invalid), ...
                                fliplr(yl(2)*ones(size(rc.C)) .* invalid - 1e3*(~invalid))], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    yline(chenRatio, 'r-', 'LineWidth', 1.5, 'Label', 'Chen 4.3');
    yline(chenRatioCIlow,  'r--');
    yline(chenRatioCIhigh, 'r--');
    xline(IC50_Kr, 'b:', 'Label', 'IC_{50}^{Kr}');
    if ~isnan(xrefChen(ct))
        plot(xrefChen(ct), chenRatio, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
    end
    set(gca, 'XScale', 'log', 'FontSize', 12);
    xlabel('[ATO] (\mumol/L)');
    ylabel('R(C) = (dAPD_F/dC) / (dAPD_M/dC)');
    title(ctNames{ct});
    ylim(yl);
end
sgtitle('Local slope ratio R(C); star = Chen X_{ref}^{Chen}');

%% Figure 2: Aggregate prevalence per cell type with Mumford targets marked
fig2 = figure('Position', [100 100 1200 350]); clf
mumfordColors = lines(length(mumfordPrev));
TStyles = {'-', '--'};
for ct = 1:nCT
    subplot(1, nCT, ct);
    pc = prevalenceCurves{ct};
    for s = 1:length(T_sens)
        semilogx(pc.C, pc.prevAgg(s, :), TStyles{s}, 'Color', 'k', ...
            'LineWidth', 1.5); hold on
    end
    for g = 1:length(mumfordPrev)
        yline(mumfordPrev(g), '-', 'Color', mumfordColors(g, :), 'LineWidth', 1.0);
        if ~isnan(mumfordXref(ct, g, s_central))
            plot(mumfordXref(ct, g, s_central), mumfordPrev(g), 'o', ...
                'MarkerSize', 9, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', mumfordColors(g, :));
        end
    end
    set(gca, 'XScale', 'log', 'FontSize', 12);
    xlabel('[ATO] (\mumol/L)');
    ylabel('Aggregate prevalence (%)');
    title(ctNames{ct});
    if ct == 1
        legend(arrayfun(@(t) sprintf('T=%d ms', t), T_sens, 'Uniform', false), ...
            'Location', 'northwest', 'FontSize', 9);
    end
end
sgtitle('Aggregate prevalence of \DeltaAPD_{90} > T; markers = Mumford 3.9/11.1/20.6%');

%% Figure 3: Bridge plot (water, cell) per cell type with linear-fit line
fig3 = figure('Position', [100 100 1200 350]); clf
for ct = 1:nCT
    subplot(1, nCT, ct);
    waters = [chenWaterMean, mumfordWater];
    cells  = [xrefChen(ct), reshape(mumfordXref(ct, :, s_central), 1, [])];

    % Plot points
    valid = ~isnan(cells);
    loglog(waters(1), cells(1), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r'); hold on
    loglog(waters(2:4), cells(2:4), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    % Linear fit
    if ~isnan(kFit(ct))
        wRange = logspace(log10(min(waters)*0.3), log10(max(waters)*3), 200);
        loglog(wRange, kFit(ct) * wRange, 'k-', 'LineWidth', 1.2);
    end
    set(gca, 'FontSize', 12);
    xlabel('Water arsenic (\mug/L)');
    ylabel('Cell free arsenite (\mumol/L)');
    title(sprintf('%s  (k = %.2e)', ctNames{ct}, kFit(ct)));
    if ct == 1
        legend({'Chen point', 'Mumford points', 'Linear fit'}, ...
            'Location', 'southeast', 'FontSize', 9);
    end
end
sgtitle('Water -> cardiomyocyte mapping; linear-through-origin fit per cell type');

%% Helpers
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
