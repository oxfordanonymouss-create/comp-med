% Population validation: 25 male + 25 female virtual cells per cell type
% with variability in 11 ionic parameters, simulated across the arsenic
% concentration range. Compares population-level metrics against clinical
% data from Mumford et al. (2007) and Chen et al. (2013).
%
% Parallel version: the (pop, concentration, cell-type) grid is flattened
% and distributed across workers via parfor.
%% Setting parameters
clear

beats = 100;
showLastN = 2;
ignoreFirst = beats - showLastN;
bcl = 1000;
options = [];
nPop = 25; % virtual cells per sex per cell type

% IC50 values (umol/L) — Drolet et al. 2004
IC50_Kr = 0.14;
IC50_Ks = 1.13;
h = 1; % Hill coefficient

% Concentration range: same as single-cell sweep
nConcPerChannel = 10;
concKr = logspace(log10(0.01 * IC50_Kr), log10(10 * IC50_Kr), nConcPerChannel);
concKs = logspace(log10(0.01 * IC50_Ks), log10(10 * IC50_Ks), nConcPerChannel);
concentrations = [0, unique([concKr, concKs])];
nConc = length(concentrations);

cellTypes = [0, 1, 2]; % endo, epi, mid
cellTypeNames = {'Endo', 'Epi', 'Mid'};
nCT = length(cellTypes);

%% Generate population via Latin hypercube sampling
paramNames = {'IKr_Multiplier', 'IKs_Multiplier', 'IK1_Multiplier', ...
    'INa_Multiplier', 'INaL_Multiplier', 'ICaL_Multiplier', ...
    'Ito_Multiplier', 'INaCa_Multiplier', 'INaK_Multiplier', ...
    'Jrel_Multiplier', 'Jup_Multiplier'};
nParams = length(paramNames);

% LHS in [0,1], scaled to [0.5, 1.5] (±50% around baseline)
lhsSamples = lhsdesign(nPop, nParams);
multipliers = 0.5 + lhsSamples; % nPop x nParams

%% Start parallel pool (4 workers)
if isempty(gcp('nocreate'))
    parpool('local', 8);
end

%% Simulation loop (flattened grid, parallelized)
tic
nTotal = nPop * nConc * nCT;

% Precompute index mapping so each parfor iter knows its (j, i, ct)
[jIdx, iIdx, ctIdx] = ind2sub([nPop, nConc, nCT], 1:nTotal);

apd90M_flat = zeros(nTotal, 1);
apd90F_flat = zeros(nTotal, 1);
tracesM_flat = cell(nTotal, 1);
tracesF_flat = cell(nTotal, 1);

parfor idx = 1:nTotal
    j = jIdx(idx);
    i = iIdx(idx);
    ct = ctIdx(idx);
    X = concentrations(i);

    paramM = struct();
    paramM.bcl = bcl;
    paramM.model = @model_ToRORd_Land_Male;
    paramM.cellType = cellTypes(ct);

    paramF = struct();
    paramF.bcl = bcl;
    paramF.model = @model_ToRORd_Land_Female;
    paramF.cellType = cellTypes(ct);

    % Set the 11 ionic multipliers from LHS
    for p = 1:nParams
        paramM.(paramNames{p}) = multipliers(j, p);
        paramF.(paramNames{p}) = multipliers(j, p);
    end

    % Drug block on top of cell-specific variability
    fracKr = 1 / (1 + (X / IC50_Kr)^h);
    fracKs = 1 / (1 + (X / IC50_Ks)^h);
    paramM.IKr_Multiplier = multipliers(j, 1) * fracKr;
    paramM.IKs_Multiplier = multipliers(j, 2) * fracKs;
    paramF.IKr_Multiplier = paramM.IKr_Multiplier;
    paramF.IKs_Multiplier = paramM.IKs_Multiplier;

    X0 = getStartingState('m_endo');

    [tM, XM] = modelRunner(X0, options, paramM, beats, ignoreFirst);
    cM = getCurrentsStructure(tM, XM, beats, paramM, 0);
    apd90M_flat(idx) = getAPD(cM.time, cM.V, 0.9, bcl);
    tracesM_flat{idx} = struct('time', cM.time, 'V', cM.V);

    [tF, XF] = modelRunner(X0, options, paramF, beats, ignoreFirst);
    cF = getCurrentsStructure(tF, XF, beats, paramF, 0);
    apd90F_flat(idx) = getAPD(cF.time, cF.V, 0.9, bcl);
    tracesF_flat{idx} = struct('time', cF.time, 'V', cF.V);

    fprintf('%s — Cell %d/%d, [ATO] = %.4f uM done\n', ...
        cellTypeNames{ct}, j, nPop, X);
end

% Reshape flat results back to (nPop, nConc, nCT)
apd90M = reshape(apd90M_flat, [nPop, nConc, nCT]);
apd90F = reshape(apd90F_flat, [nPop, nConc, nCT]);
tracesM = reshape(tracesM_flat, [nPop, nConc, nCT]);
tracesF = reshape(tracesF_flat, [nPop, nConc, nCT]);

fprintf('Simulation time: %.1f seconds\n', toc);

%% Compute dAPD90 relative to each cell's own baseline
dAPD90M = apd90M - apd90M(:, 1, :);
dAPD90F = apd90F - apd90F(:, 1, :);

%% Persist sweep outputs for downstream analysis
save('population_sweep.mat', 'apd90M', 'apd90F', 'concentrations', ...
    'cellTypeNames', 'tracesM', 'tracesF');

%% Helper functions
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

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
