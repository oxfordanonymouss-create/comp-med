% Verification script for the sex-specific ToR-ORd-Land model.
% Runs baseline male and female simulations for endo, epi, and mid cell
% types at BCL = 1000 ms. Plots APs side by side and compares APD90 values
% against Holmes et al. (2025).
%% Setting parameters
clear

cellTypes = [0, 1, 2]; % 0 = endo, 1 = epi, 2 = mid
cellTypeNames = {'Endo', 'Epi', 'Mid'};

% Male parameters
paramM.bcl = 1000;
paramM.model = @model_ToRORd_Land_Male;
paramM.IKr_Multiplier = 1;
paramM.cellType = 0;

paramsM(1:length(cellTypes)) = paramM;
for i = 1:length(cellTypes)
    paramsM(i).cellType = cellTypes(i);
end

% Female parameters
paramF.bcl = 1000;
paramF.model = @model_ToRORd_Land_Female;
paramF.IKr_Multiplier = 1;
paramF.cellType = 0;


paramsF(1:length(cellTypes)) = paramF;
for i = 1:length(cellTypes)
    paramsF(i).cellType = cellTypes(i);
end

options = [];
beats = 200; % the more beats, the more precise and closer to the verification results
ignoreFirst = beats - 1;

%% Simulation
n = length(cellTypes);
parfor (i = 1:n, 0)
    X0 = getStartingState('m_endo');
    [timeM{i}, XM{i}] = modelRunner(X0, options, paramsM(i), beats, ignoreFirst);
    currentsM{i} = getCurrentsStructure(timeM{i}, XM{i}, beats, paramsM(i), 0);
end
parfor (i = 1:n, 0)
    X0 = getStartingState('m_endo');
    [timeF{i}, XF{i}] = modelRunner(X0, options, paramsF(i), beats, ignoreFirst);
    currentsF{i} = getCurrentsStructure(timeF{i}, XF{i}, beats, paramsF(i), 0);
end

%% Compute APD90
apd90M = zeros(1, n);
apd90F = zeros(1, n);
for i = 1:n
    apd90M(i) = getAPD(currentsM{i}.time, currentsM{i}.V, 0.9);
    apd90F(i) = getAPD(currentsF{i}.time, currentsF{i}.V, 0.9);
end

%% Plot APs — one figure per cell type, male vs female
for i = 1:n
    figure(i); clf
    plot(currentsM{i}.time, currentsM{i}.V, 'b', 'LineWidth', 2); hold on
    plot(currentsF{i}.time, currentsF{i}.V, 'r', 'LineWidth', 2);
    title(sprintf('%s — Action Potential', cellTypeNames{i}));
    legend(sprintf('Male (APD_{90} = %.1f ms)', apd90M(i)), ...
           sprintf('Female (APD_{90} = %.1f ms)', apd90F(i)));
    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');
    xlim([-25 600]);
    ylim([-100 50]);
end

%% Plot APD90 bar chart
figure(n+1); clf
barData = [apd90M; apd90F]';
b = bar(barData);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
set(gca, 'XTickLabel', cellTypeNames);
legend('Male', 'Female');
xlabel('Cell type');
ylabel('APD_{90} (ms)');
title('Baseline APD_{90} — Male vs Female');

%% Print APD90
fprintf('%-6s  Male APD90   Female APD90   Diff\n', 'Cell');
for i = 1:n
    pctDiff = (apd90F(i) - apd90M(i)) / apd90M(i) * 100;
    fprintf('%-6s  %.1f ms      %.1f ms        %+.1f%%\n', cellTypeNames{i}, apd90M(i), apd90F(i), pctDiff);
end

%% Helper function
function apd = getAPD(time, V, level)
    baseline  = V(end);
    threshold = baseline + (1 - level) * (max(V) - baseline);
    above = V > threshold;
    startIdx = find(above, 1, 'first');
    endIdx = find(above & (1:length(above))' >= startIdx, 1, 'last');
    apd = time(endIdx) - time(startIdx);
end
