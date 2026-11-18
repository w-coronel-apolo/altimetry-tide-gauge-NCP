%% Wind and discharge correction
% Select predictors using detrended records, then retain their trends
% when calculating the near-coastal correction.
% A stores corrected series/trends; Models stores fits and local PCA bases.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'functions'));
S = load(fullfile(root, 'outputs', 'results.mat'));
O = S.O;
R = S.R;
W = load(fullfile(root, 'data', 'forcing.mat'));
M = O.M;
t = O.t;
N = O.N;
L = O.L;
ALT = O.ALT;
MRD = R.MRD;
LRD = R.LRD;
NRD = R.NRD;
tRD = R.tRD;
raw = O.ALT_TG;
n = numel(N);
nt = numel(t);
fc = W.FC_monthly(W.tFC >= 1993 & W.tFC < 2025);
mapping = O.river_mapping;
[~, row] = ismember(station_key(N), station_key(mapping.Station));
[usedQ, qi] = ismember(station_key(mapping.RiverName(row)), station_key(NRD));
q = double(MRD);
qdt = double(detrend(MRD, 'omitnan'));
box = LRD(:, 1) >= -80 & LRD(:, 1) <= -74 & LRD(:, 2) >= 36 & LRD(:, 2) <= 40;
qsum = sum(MRD(:, box), 2, 'omitnan');
qsum(all(~isfinite(MRD(:, box)), 2)) = NaN;
qsumdt = double(detrend(qsum, 'omitnan'));
qsum = double(qsum);
fcdt = double(detrend(fc, 'omitnan'));
fc = double(fc);
ydt = detrend(raw, 'omitnan');
% Monthly anomalies, with and without detrending, for Q and Florida Current.
for m = 1:12
    q(m:12:end, :) = q(m:12:end, :) - mean(q(m:12:end, :), 1, 'omitnan');
    qdt(m:12:end, :) = qdt(m:12:end, :) - mean(qdt(m:12:end, :), 1, 'omitnan');
    qsum(m:12:end) = qsum(m:12:end) - mean(qsum(m:12:end), 'omitnan');
    qsumdt(m:12:end) = qsumdt(m:12:end) - mean(qsumdt(m:12:end), 'omitnan');
    fc(m:12:end) = fc(m:12:end) - mean(fc(m:12:end), 'omitnan');
    fcdt(m:12:end) = fcdt(m:12:end) - mean(fcdt(m:12:end), 'omitnan');
end
[lon, lat] = ndgrid(W.lon, W.lat);
U = reshape(W.U10, [], nt)';
V = reshape(W.V10, [], nt)';
%% Estuary groups and Florida Current stations
N = string(N(:));
between = @(first, last)(1:numel(N))' >= find(N == first) & (1:numel(N))' <= find(N == last);
groups.stlaw = between('ALEXANDRIA BAY', 'RIVIERE-AU-RENARD');
groups.delaware = between('CAPE MAY', 'NEWBOLD');
groups.chesapeake = between('KIPTOPEKE BEACH', 'CHESAPEAKE CITY');
groups.estuary = groups.stlaw | groups.delaware | groups.chesapeake;
groups.outside = ~groups.estuary;
groups.fc = between('FERNANDINA BEACH', 'VACA KEY');
A = struct('t', t, 'N', N, 'L', L, 'ALT_TG', raw, 'prediction', nan(nt, n), ...
    'corrected', nan(nt, n));
fields = {'usedQ', 'usedQsum', 'usedFC', 'selectedFC', 'nPC_U', 'nPC_V', 'candidate_count', ...
    'selected_count', 'adjusted_r2', 'q_col', 'n_grid', 'corr_q', 'corr_q_sum', 'corr_fc'};
for f = fields
    A.(f{1}) = nan(n, 1);
end
%% One local model per TG
% Wind box: +/-4 degrees; retain PCs explaining >=95% variance.
% The fit requires >=180 complete months, with actual sample times retained.
models = cell(n, 1);
predictors = cell(n, 1);
for j = 1:n
    grid = find(abs(lon - L(j, 1)) <= 4 & abs(lat - L(j, 2)) <= 4);
    u = wind_pca(U(:, grid));
    v = wind_pca(V(:, grid));
    nu = size(u.scores, 2);
    nv = size(v.scores, 2);
    X = [u.scores, v.scores];
    labels = ["U" + (1:nu), "V" + (1:nv)];
    cq = NaN;
    cs = NaN;
    cf = NaN;
    sumChoice = false;
    if usedQ(j)
        cq = corr(ydt(:, j), qdt(:, qi(j)), 'Rows', 'pairwise');
    end
    if groups.chesapeake(j)
        cs = corr(ydt(:, j), qsumdt, 'Rows', 'pairwise');
    end
    if usedQ(j)
        sumChoice = isfinite(cs) && (~isfinite(cq) || abs(cs) > abs(cq));
        if sumChoice
            X = [X, qsum];
        else
            X = [X, q(:, qi(j))];
        end
        labels = [labels, "Q"];
    end
    if groups.fc(j)
        X = [X, fc];
        labels = [labels, "FC"];
        cf = corr(ydt(:, j), fcdt, 'Rows', 'pairwise');
    end
    F = correction_fit(raw(:, j), X, t);
    A.prediction(:, j) = F.prediction;
    A.corrected(:, j) = F.corrected;
    qindex = NaN;
    if usedQ(j)
        qindex = qi(j);
    end
    values = [usedQ(j), sumChoice, groups.fc(j), any(labels(F.selected) == "FC"), nu, nv, ...
        size(X, 2), numel(F.selected), F.adjr2, qindex, numel(grid), cq, cs, cf];
    for k = 1:numel(fields)
        A.(fields{k})(j) = values(k);
    end
    predictors{j} = struct('X', X, 'names', labels, 'u', u, 'v', v, 'grid', grid);
    models{j} = F;
    fprintf('%3d/131 %-24s %d predictors\n', j, N(j), numel(F.selected));
end

%% OLS trends and HAC errors (mm/yr)
% Full records require 60 months; paired analyses require 180.
% HAC12/24/36 use calendar-month pairs, not adjacent nonmissing observations.
common = isfinite(raw) & isfinite(A.corrected);
paired = raw;
paired(~common) = NaN;
Y = {raw, A.corrected, paired, A.corrected};
labels = {'raw', 'corrected', 'paired_raw', 'paired_corrected'};
for k = 1:4
    label = labels{k};
    lags = [12 24 36];
    h = 2;
    minN = 60;
    if k >= 3
        lags = 24;
        h = 1;
        minN = 180;
    end
    rate = nan(n, 1);
    se = rate;
    nobs = rate;
    errors = nan(n, numel(lags));
    for j = 1:n
        F = trend_hac(t, Y{k}(:, j), lags, minN);
        rate(j) = F.rate;
        errors(j, :) = F.se;
        se(j) = F.se(h);
        nobs(j) = F.n;
    end
    A.([label '_rate']) = rate;
    A.([label '_se']) = se;
    A.([label '_hac_all']) = errors;
    A.([label '_nobs']) = nobs;
end
A.common_count = sum(common, 1)';
A.first_year = nan(n, 1);
A.last_year = nan(n, 1);
for j = 1:n
    times = t(common(:, j));
    A.first_year(j) = times(1);
    A.last_year(j) = times(end);
end
A.final_mask = A.common_count >= 180 & isfinite(A.paired_corrected_rate) & isfinite(A.paired_corrected_se);
A.groups = groups;
Models = struct('models', {models}, 'predictors', {predictors});
%% Wind-only and discharge-only models
% Figures 3/5 use shared months BEFORE detrending; S4/S5 detrend first.
% Keep these distinct because selecting months changes the fitted trends.
t = A.t;
raw = A.ALT_TG;
[nt, n] = size(raw);
groups = A.groups;
total = qsum; % Reuse the same deseasonalized Chesapeake sum and q matrix.
representatives = ["TROIS-RIVIERES", "NEUVILLE", "QUEBEC", "NEW LONDON", "KINGS POINT", ...
    "REEDY POINT", "PHILADELPHIA", "NEWBOLD", "MONEY POINT", "WASHINGTON DC", ...
    "BALTIMORE", "WILMINGTON", "APALACHICOLA", "BAY WAVELAND YACHT CLUB", "LAKE CHARLES"];
C.common_adj_r2 = nan(n, 3);
C.common_ordinary_r2 = nan(n, 3);
C.common_n_selected = nan(n, 3);
C.common_n_obs = zeros(n, 1);
C.common_used_q = false(n, 1);
C.common_used_q_sum = false(n, 1);
for profile = ["sensitivity", "timeseries"]
    C.(profile + "_prediction_uv") = nan(nt, n);
    C.(profile + "_prediction_q") = nan(nt, n);
    C.(profile + "_adj_r2") = nan(n, 2);
    C.(profile + "_n_obs") = zeros(n, 2);
    C.(profile + "_n_selected") = nan(n, 2);
end
for j = 1:n
    p = predictors{j};
    y = raw(:, j);
    uv = p.X(:, startsWith(p.names, 'U') | startsWith(p.names, 'V'));
    savedQ = p.X(:, p.names == "Q");
    localQ = zeros(nt, 0);
    commonQ = localQ;
    if A.usedQ(j)
        localQ = q(:, A.q_col(j));
        commonQ = localQ;
    end
    if groups.chesapeake(j)
        yd = detrend(y, 'omitnan', 'SamplePoints', t);
        rlocal = NaN;
        if ~isempty(localQ)
            rlocal = corr(yd, detrend(localQ, 'omitnan', 'SamplePoints', t), 'Rows', 'pairwise');
        end
        rsum = corr(yd, detrend(total, 'omitnan', 'SamplePoints', t), 'Rows', 'pairwise');
        if A.usedQ(j) && abs(rsum) > abs(rlocal)
            commonQ = total;
            C.common_used_q_sum(j) = true;
        end
        if isfinite(rsum) && (~isfinite(rlocal) || abs(rsum) > abs(rlocal))
            localQ = total;
        end
    end

    % Main Figures 3/5: common months first, then actual-time detrending.
    allX = [uv, commonQ];
    ok = isfinite(y) & all(isfinite(allX), 2);
    C.common_n_obs(j) = nnz(ok);
    C.common_used_q(j) = ~isempty(commonQ);
    if nnz(ok) >= 180
        Y = detrend(y(ok), 'omitnan', 'SamplePoints', t(ok));
        matrices = {uv, commonQ, allX};
        for k = 1:3
            if isempty(matrices{k})
                continue;
            end
            F = stepwise_fit(detrend(matrices{k}(ok, :), 'omitnan', 'SamplePoints', t(ok)), Y);
            C.common_adj_r2(j, k) = F.adjr2;
            C.common_ordinary_r2(j, k) = F.r2;
            C.common_n_selected(j, k) = numel(F.selected);
        end
    end

    % S5: preserve monthly time gaps, as in the full correction model.
    matrices = {uv, savedQ};
    labels = ["uv", "q"];
    for k = 1:2
        if isempty(matrices{k})
            continue;
        end
        F = correction_fit(y, matrices{k}, t);
        C.("sensitivity_prediction_" + labels(k))(:, j) = F.prediction;
        C.sensitivity_adj_r2(j, k) = F.adjr2;
        C.sensitivity_n_obs(j, k) = nnz(F.fit_mask);
        C.sensitivity_n_selected(j, k) = numel(F.selected);
    end

    % S4: detrend each full series in actual time before finding complete rows.
    if ~ismember(A.N(j), representatives)
        continue;
    end
    matrices = {uv, localQ};
    for k = 1:2
        X = matrices{k};
        if isempty(X)
            continue;
        end
        F = correction_fit(y, X, t);
        % Center over all predictable rows before imposing the observed-y mask.
        xp = X - mean(X, 1, 'omitnan');
        ok = all(isfinite(xp), 2);
        pred = [ones(nnz(ok), 1), xp(ok, F.selected)] * F.b;
        F.prediction(:) = NaN;
        F.prediction(ok) = pred - mean(pred);
        F.prediction(~isfinite(y)) = NaN;
        C.("timeseries_prediction_" + labels(k))(:, j) = F.prediction;
        C.timeseries_adj_r2(j, k) = F.adjr2;
        C.timeseries_n_obs(j, k) = nnz(F.fit_mask);
        C.timeseries_n_selected(j, k) = numel(F.selected);
    end
end
save(fullfile(root, 'outputs', 'results.mat'), 'A', 'Models', 'C', '-append');

%% Local calculations
function P = wind_pca(X)
    % Center wind columns, fill gaps linearly, and retain >=95% variance.
    % MATLAB pca determines the loading signs; keep columns unscaled.
    % Remove calendar-month means before PCA; inputs begin in January.
    X = double(X);
    for m = 1:12
        X(m:12:end, :) = X(m:12:end, :) - mean(X(m:12:end, :), 1, 'omitnan');
    end
    P.columns = find(any(isfinite(X), 1));
    X = X(:, P.columns);
    X = fillmissing(X, 'linear', 'EndValues', 'nearest');
    P.mean = mean(X);
    [coeff, score, e] = pca(X);
    P.explained = 100 * e / sum(e);
    k = find(cumsum(P.explained) >= 95, 1);
    P.loadings = coeff(:, 1:k);
    P.scores = score(:, 1:k);
    P.eigenvalues = e;
end

function F = correction_fit(y, X, t)
    yf = detrend(y, 'omitnan', 'SamplePoints', t);
    xf = detrend(X, 'omitnan', 'SamplePoints', t);
    ok = isfinite(yf) & all(isfinite(xf), 2);
    assert(nnz(ok) >= 180, 'Fewer than 180 complete fitting months.');
    F = stepwise_fit(xf(ok, :), yf(ok));
    F.fit_mask = ok;
    A = [ones(numel(y), 1), X(:, F.selected) - mean(X(:, F.selected), 1, 'omitnan')];
    F.prediction = A * F.b;
    F.prediction = F.prediction - mean(F.prediction, 'omitnan');
    F.corrected = y - F.prediction;
end

function F = stepwise_fit(X, y)
    % MATLAB partial-F selection; keep coefficient order aligned with X.
    mdl = stepwiselm(X, y, 'constant', 'Upper', 'linear', 'Criterion', 'sse', ...
        'PEnter', .05, 'PRemove', .10, 'Verbose', 0);
    selected = find(any(mdl.Formula.Terms(:, 1:size(X, 2)) ~= 0, 1));
    [~, order] = ismember(["(Intercept)", "x" + selected], string(mdl.CoefficientNames));
    F = struct('selected', selected, 'b', mdl.Coefficients.Estimate(order), ...
        'covariance', mdl.CoefficientCovariance(order, order), 'df', mdl.DFE, ...
        'sse', mdl.SSE, 'r2', mdl.Rsquared.Ordinary, 'adjr2', mdl.Rsquared.Adjusted);
end

function key = station_key(name)
    key = upper(string(name));
    key = erase(key, ["'", "."]);
    key = regexprep(key, '\s*\([^)]*\)', '');
    key = regexprep(key, ',.*$', '');
    key = strtrim(regexprep(key, '\s+', ' '));
end
