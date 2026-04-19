% Single-cell concentration-response sweep: does the model predict greater
% dAPD90 in female than male models across the full arsenic concentration
% range?
%
% IKr and IKs are blocked via the Hill equation using IC50 values from
% Drolet et al. (2004). 20 log-spaced concentrations are constructed by
% merging two panels (0.01-10x each IC50).
%
% Parallel version: the (concentration, cell-type) grid is flattened and
% distributed across workers via parfor.
%% Setting parameters
clear

beats = 100;
showLastN = 2; % number of final beats to keep
ignoreFirst = beats - showLastN;
bcl = 1000;
options = [];

% IC50 values (u mol/L) — Drolet et al. 2004
IC50_Kr = 0.14;
IC50_Ks = 1.13;
h = 1; % Hill coefficient

% Concentration range: 10 log-spaced per channel, merged, plus 0 (baseline)
nConc = 20;
concKr = logspace(log10(0.01 * IC50_Kr), log10(10 * IC50_Kr), nConc/2);
concKs = logspace(log10(0.01 * IC50_Ks), log10(10 * IC50_Ks), nConc/2);
concentrations = [0, unique([concKr, concKs])];
nConc = length(concentrations);

cellTypes = [0, 1, 2]; % endo, epi, mid
cellTypeNames = {'Endocardial', 'Epicardial', 'Midmyocardial'};
nCT = length(cellTypes);

%% Start parallel pool (4 workers)
if isempty(gcp('nocreate'))
    parpool('local', 4);
end

%% Sweep concentrations (flattened grid, parallelized)
tic
nTotal = nConc * nCT;

% Precompute index mapping so each parfor iter knows its (i, ct)
[iIdx, ctIdx] = ind2sub([nConc, nCT], 1:nTotal);

apd90M_flat = zeros(nTotal, 1);
apd90F_flat = zeros(nTotal, 1);
traceM_flat = cell(nTotal, 1);
traceF_flat = cell(nTotal, 1);

parfor idx = 1:nTotal
    i = iIdx(idx);
    ct = ctIdx(idx);
    X = concentrations(i);

    paramM = struct();
    paramM.bcl = bcl;
    paramM.model = @model_ToRORd_Land_Male;
    paramM.cellType = cellTypes(ct);
    paramM.IKr_Multiplier = 1 / (1 + (X / IC50_Kr)^h);
    paramM.IKs_Multiplier = 1 / (1 + (X / IC50_Ks)^h);

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

    fprintf('[ATO] = %.4f uM | %s | IKr block = %.1f%% | IKs block = %.1f%% | Male: %.1f ms | Female: %.1f ms\n', ...
        X, cellTypeNames{ct}, (1 - paramM.IKr_Multiplier)*100, ...
        (1 - paramM.IKs_Multiplier)*100, apd90M_flat(idx), apd90F_flat(idx));
end

% Reshape flat results back to (nConc, nCT)
apd90M = reshape(apd90M_flat, [nConc, nCT]);
apd90F = reshape(apd90F_flat, [nConc, nCT]);
traceM = reshape(traceM_flat, [nConc, nCT]);
traceF = reshape(traceF_flat, [nConc, nCT]);

fprintf('Simulation time: %.1f seconds\n', toc);

%% Compute dAPD90 relative to baseline (first row = concentration 0)
dAPD90M = apd90M - apd90M(1, :);
dAPD90F = apd90F - apd90F(1, :);

%% Persist sweep outputs for downstream analysis (scriptConcentrationMapping.m)
save('singlecell_sweep.mat', 'apd90M', 'apd90F', 'concentrations', ...
    'cellTypeNames', 'traceM', 'traceF');

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
