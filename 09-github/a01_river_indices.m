%% River indices, January 1993-December 2024
% Q contains individual monthly gauges, already screened at >=80% daily coverage.
% MRD contains river indices (m3/s); LRD is [lon lat], NRD names, tRD decimal years.
root = fileparts(mfilename('fullpath'));
folder = fullfile(root, 'outputs');
if ~isfolder(folder)
    mkdir(folder);
end
R = load(fullfile(root, 'data', 'river_monthly.mat'));
Q = R.Q;
months = R.months;
members = R.members;
NRD = string(R.NRD(:));
LRD = R.LRD;
tRD = R.tRD;
MRD = nan(numel(months), numel(NRD));
contributors = MRD;
% Average log anomalies, then restore the downstream gauge's mean log flow.
% Use each gauge's available 1993-2024 months; a single-gauge index is unchanged.
for j = 1:numel(NRD)
    group = find(string(members.river) == strtrim(NRD(j)));
    reference = group(ismember(lower(string(members.reference(group))), ["true", "1"]));
    z = Q(:, group); % Gauge means use available screened observations in 1993-2024.
    if numel(group) == 1
        index = z;
    else
        z(z <= 0) = NaN;
        logs = log10(z);
        center = mean(logs, 1, 'omitnan');
        anomaly = logs - center;
        combined = mean(anomaly, 2, 'omitnan');
        index = 10 .^ (combined + center(group == reference));
    end
    MRD(:, j) = index;
    contributors(:, j) = sum(isfinite(z), 2);
end
R = struct('MRD', MRD, 'NRD', NRD, 'LRD', R.LRD, 'tRD', R.tRD, ...
    'reference_sites', R.reference_sites, 'units', R.units, 'months', months, ...
    'coverage', R.coverage, 'contributors', contributors);
save(fullfile(folder, 'results.mat'), 'R', '-v7.3');
T = array2table(MRD, 'VariableNames', matlab.lang.makeValidName(NRD));
T = addvars(T, months, 'Before', 1, 'NewVariableNames', 'month');
writetable(T, fullfile(folder, 'river_indices.csv'));
fprintf('River indices: %d gauges, %d rivers, m3/s.\n', height(members), numel(NRD));
