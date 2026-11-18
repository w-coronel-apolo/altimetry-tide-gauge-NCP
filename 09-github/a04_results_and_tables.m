%% Results and supplementary tables
% Compare with GNSS/GIA; calculate paired variability and trend sensitivities.
% Trends/errors: mm/yr, one HAC24 SE. Coefficient errors: classical OLS SE.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'functions'));
S = load(fullfile(root, 'outputs', 'results.mat'));
A = S.A;
O = S.O;
Models = S.Models;
R = S.R;
C = S.C;
K = load(fullfile(root, 'data', 'ancillary.mat'));
M = O.M;
L = O.L;
N = O.N;
t = O.t;
ALT = O.ALT;
n = numel(N);
MRD = R.MRD;
LRD = R.LRD;
NRD = R.NRD;
tRD = R.tRD;

%% GNSS quality and GIA at the TG coordinates
% Quality thresholds come from the global NGL table, not this regional subset.
D = struct;
T = K.gnss;
[~, j] = ismember(station_key(A.N), station_key(T.Station));
T = T(j, :);
globalMetrics = K.gnss_global_metrics;
metrics = {'su', 'eta', 'zeta', 'ups'};
values = T{:, metrics};
good = zeros(1, 4);
poor = good;
qgood = [66.6667 33.3333 66.6667 66.6667];
qpoor = [95 5 95 95];
for k = 1:4
    good(k) = prctile(globalMetrics.(metrics{k}), qgood(k));
    poor(k) = prctile(globalMetrics.(metrics{k}), qpoor(k));
end
available = isfinite(T.rate);
complete = available & all(isfinite(values), 2);
quality = repmat("Missing", height(T), 1);
quality(available) = "Unclassified";
quality(complete) = "Medium";
quality(complete & values(:, 1) < good(1) & values(:, 2) > good(2) & values(:, 3) < good(3) & values(:, 4) < good(4)) = "Good";
quality(complete & (values(:, 1) >= poor(1) | values(:, 2) <= poor(2) | values(:, 3) >= poor(3) | values(:, 4) >= poor(4))) = "Poor";
separation = nan(height(T), 1);
for j = 1:height(T)
    separation(j) = distance_km(A.L(j, :), [T.lon(j), T.lat(j)]);
end
D.gnss = struct('rate', T.rate, 'se', T.su, 'metrics', values, 'quality', quality, ...
    'good', quality == "Good", 'matched', available, 'coordinates', [T.lon T.lat], ...
    'distance_km', separation, 'threshold_good', good, 'threshold_poor', poor);
assert(~available(A.N == "NEW HAVEN"), 'New Haven must remain unmatched.');
portland = A.N == "PORTLAND";
assert(T.lon(portland) > -71 && T.lon(portland) < -69 && T.lat(portland) > 43 && T.lat(portland) < 45);
D.groups = A.groups;
fields = {'mean', 'std'};
for k = 1:2
    grid = K.gia.(fields{k});
    x = double(grid.x);
    y = double(grid.y);
    z = double(grid.z);
    x(x > 180) = x(x > 180) - 360;
    [x, order] = sort(x);
    z = z(order, :);
    [x, uniqueRows] = unique(x, 'stable');
    z = z(uniqueRows, :);
    F = griddedInterpolant({x, y}, z, 'linear', 'none');
    D.gia.(fields{k}) = F(A.L(:, 1), A.L(:, 2));
end
D.gia.gulf = A.L(:, 1) <= -81.5 & A.L(:, 2) <= 30.75;
D.gia.atlantic = ~D.gia.gulf;
D.gia.atlantic_profile = gia_profile(A.L(D.gia.atlantic, 2), D.gia.mean(D.gia.atlantic), D.gia.std(D.gia.atlantic));
D.gia.gulf_profile = gia_profile(-A.L(D.gia.gulf, 1), D.gia.mean(D.gia.gulf), D.gia.std(D.gia.gulf));
D.cape_may = cape_may_data(A, O, K, root);

%% Paired-month variability and TG/ALT coherence
r = nan(n, 2);
sd = r;
nc = zeros(n, 1);
for j = 1:n
    ok = isfinite(A.ALT_TG(:, j)) & isfinite(A.corrected(:, j));
    nc(j) = nnz(ok);
    sd(j, :) = std(detrend([A.ALT_TG(ok, j), A.corrected(ok, j)], 'omitnan', 'SamplePoints', t(ok)), 0, 1);
    adjusted = O.ALT(:, j) - A.prediction(:, j);
    y = detrend([O.Mds(ok, j), O.ALT(ok, j), adjusted(ok)], 'omitnan', 'SamplePoints', t(ok));
    r(j, :) = corr(y(:, 1), y(:, 2:3), 'Rows', 'pairwise');
end
delta = r(:, 2) - r(:, 1);
status = repmat("unchanged", n, 1);
status(delta > 1e-10) = "increased";
status(delta < -1e-10) = "decreased";
coherence = table((0:n - 1)', A.N, nc, r(:, 1), r(:, 2), delta, status, A.selected_count, ...
    sd(:, 1), sd(:, 2), 100 * (1 - sd(:, 2) ./ sd(:, 1)), ...
    'VariableNames', {'station_index', 'station', 'n_common', 'r_raw', 'r_corrected', 'delta_r', 'status', ...
    'n_selected', 'sd_raw', 'sd_corrected', 'sd_reduction_percent'});
D.coherence = coherence;

%% Expanding 10-28 year windows, anchored in January 1993
rows = cell(19 * n, 10);
row = 0;
for years = 10:28
    count = 12 * years;
    for j = 1:n
        row = row + 1;
        raw = A.ALT_TG(1:count, j);
        corrected = A.corrected(1:count, j);
        ok = isfinite(raw) & isfinite(corrected);
        errors = [NaN NaN];
        rates = errors;
        ratio = NaN;
        if nnz(ok) >= ceil(.8 * count)
            times = t(1:count);
            a = trend_hac(times(ok), raw(ok), 24, 3);
            b = trend_hac(times(ok), corrected(ok), 24, 3);
            errors = [a.se b.se];
            rates = [a.rate b.rate];
            ratio = b.se / a.se;
        end
        rows(row, :) = {A.N(j), j - 1, years, nnz(ok), logical(A.usedQ(j)), errors(1), errors(2), rates(1), rates(2), ratio};
    end
end
windows = cell2table(rows, 'VariableNames', {'station', 'station_index', 'years', 'n_common', 'has_q', ...
    'raw_se', 'corrected_se', 'raw_rate', 'corrected_rate', 'ratio'});
D.windows = windows;

%% Full-model and individual-forcing trend adjustments on identical months
D.adjustment = [A.paired_raw_rate - A.paired_corrected_rate, nan(n, 2)];
predictions = {C.sensitivity_prediction_uv, C.sensitivity_prediction_q};
for k = 1:2
    for j = 1:n
        corrected = A.ALT_TG(:, j) - predictions{k}(:, j);
        ok = isfinite(A.ALT_TG(:, j)) & isfinite(corrected);
        if nnz(ok) < 180
            continue;
        end
        raw = trend_hac(t(ok), A.ALT_TG(ok, j), 0, 180);
        fit = trend_hac(t(ok), corrected(ok), 0, 180);
        D.adjustment(j, k + 1) = raw.rate - fit.rate;
    end
end


save(fullfile(root, 'outputs', 'results.mat'), 'D', '-append');
%% Tables S1 and S2: station rates, model summaries and coefficients
n = numel(A.N);
dates = string(datetime(1993, 1, 1) + calmonths(round((A.t - 1993) * 12)), 'yyyy-MM');
stationRows = cell(n, 1);
modelRows = cell(n, 1);
coefficientRows = {};
ches = R.LRD(:, 1) >= -80 & R.LRD(:, 1) <= -74 & R.LRD(:, 2) >= 36 & R.LRD(:, 2) <= 40;
for j = 1:n
    row = struct('station_name', A.N(j), 'coastal_index', j, ...
        'longitude_deg_east', A.L(j, 1), 'latitude_deg_north', A.L(j, 2));
    paired = isfinite(A.ALT_TG(:, j)) & isfinite(A.corrected(:, j));
    y = A.ALT_TG(:, j);
    y(~paired) = NaN;
    series = {O.Mds(:, j), A.ALT_TG(:, j), A.corrected(:, j), y};
    labels = ["tg", "alt_tg", "alt_tgncp", "paired_alt_tg"];
    F = trend_hac(A.t, O.Mds(:, j), 24, 180); % ALT-TG trends are already saved by a03.
    rates = [F.rate, A.raw_rate(j), A.paired_corrected_rate(j), A.paired_raw_rate(j)];
    errors = [F.se, A.raw_se(j), A.paired_corrected_se(j), A.paired_raw_se(j)];
    for k = 1:4
        ok = isfinite(series{k});
        label = labels(k);
        if nnz(ok) < 180
            rates(k) = NaN;
            errors(k) = NaN;
        end
        row.(label + "_rate_mm_yr") = rates(k);
        row.(label + "_hac24_se_mm_yr") = errors(k);
        row.(label + "_start_month") = dates(find(ok, 1));
        row.(label + "_end_month") = dates(find(ok, 1, 'last'));
        row.(label + "_n_months") = nnz(ok);
    end
    row.paired_raw_minus_corrected_mm_yr = A.paired_raw_rate(j) - A.corrected_rate(j);
    row.gnss_rate_mm_yr = D.gnss.rate(j);
    row.gnss_su_mm_yr = D.gnss.se(j);
    row.gnss_eta = D.gnss.metrics(j, 2);
    row.gnss_zeta_mm_yr = D.gnss.metrics(j, 3);
    row.gnss_ups_mm_yr = D.gnss.metrics(j, 4);
    row.gnss_quality = D.gnss.quality(j);
    stationRows{j} = row;

    P = Models.predictors{j};
    F = Models.models{j};
    terms = P.names(F.selected);
    ok = F.fit_mask;
    source = "";
    if A.usedQ(j)
        source = R.NRD(A.q_col(j));
    end
    if A.usedQsum(j)
        source = strjoin(R.NRD(ches), ';');
    end
    summary = struct('station_name', A.N(j), 'coastal_index', j, 'n_fit_months', nnz(ok), ...
        'fit_start_month', dates(find(ok, 1)), 'fit_end_month', dates(find(ok, 1, 'last')), ...
        'adjusted_r_squared', F.adjr2, 'n_candidate_predictors', numel(P.names), ...
        'n_selected_predictors', numel(terms), 'selected_predictors', strjoin(terms, ';'), ...
        'n_selected_u_pcs', nnz(startsWith(terms, 'U')), 'n_selected_v_pcs', nnz(startsWith(terms, 'V')), ...
        'q_offered', logical(A.usedQ(j)), 'q_selected', any(terms == "Q"), ...
        'fc_offered', logical(A.usedFC(j)), 'fc_selected', any(terms == "FC"), ...
        'discharge_sources', source, 'discharge_sum_used', logical(A.usedQsum(j)), ...
        'n_paired_months', nnz(paired), 'raw_detrended_sd_mm', D.coherence.sd_raw(j), ...
        'corrected_detrended_sd_mm', D.coherence.sd_corrected(j), 'sd_reduction_percent', D.coherence.sd_reduction_percent(j));
    modelRows{j} = summary;
    names = ["Intercept", terms];
    matlabNames = ["(Intercept)", "x" + F.selected];
    for k = 1:numel(names)
        term = names(k);
        family = term;
        unit = "Sv";
        if term == "Intercept"
            family = "intercept";
            unit = "1";
        elseif startsWith(term, 'U')
            family = "wind_U";
            unit = "m/s";
        elseif startsWith(term, 'V')
            family = "wind_V";
            unit = "m/s";
        elseif term == "Q"
            unit = "m^3/s";
        end
        coefficientUnit = "mm/(" + unit + ")";
        if unit == "1"
            coefficientUnit = "mm";
        end
        coefficientRows{end + 1} = struct('station_name', A.N(j), 'coastal_index', j, ...
            'term', term, 'matlab_term', matlabNames(k), 'family', family, 'predictor_units', unit, ...
            'coefficient', F.b(k), 'coefficient_units', coefficientUnit, ...
            'conditional_ols_se', sqrt(max(F.covariance(k, k), 0)), 'residual_df', F.df);
    end
end
S1 = struct2table(vertcat(stationRows{:}));
S2 = struct2table(vertcat(modelRows{:}));
coefficients = struct2table(vertcat(coefficientRows{:}));
folder = fullfile(root, 'outputs');
writetable(S1, fullfile(folder, 'Table_S1_station_rates.csv'));
writetable(S2, fullfile(folder, 'Table_S2_model_summary.csv'));
writetable(coefficients, fullfile(folder, 'Table_S3_regression_coefficients.csv'));
fprintf('Supplementary tables: %d stations, %d coefficient rows.\n', height(S1), height(coefficients));
fprintf('Scientific comparisons and sensitivity calculations completed.\n');

%% Local calculations
function P = gia_profile(x, m, s)
    [x, ~, group] = unique(x);
    m = accumarray(group, m, [], @mean);
    s = accumarray(group, s, [], @mean);
    grid = linspace(min(x), max(x), 500)';
    P = [grid, interp1(x, m, grid, 'pchip'), max(interp1(x, s, grid, 'linear'), 0)];
end

function S = cape_may_data(A, O, K, root)
    j = find(A.N == "CAPE MAY");
    T = K.cape_may_gnss;
    mjd = T.("__MJD");
    height = T.("__height(m)");
    site = string(T.site);
    ok = isfinite(mjd) & isfinite(height);
    assert(all(strcmpi(strtrim(site(ok)), 'CMAP')));
    mjd = mjd(ok);
    height = height(ok) * 1000;
    [mjd, order] = sort(mjd);
    height = height(order);
    dates = datetime(1858, 11, 17) + days(mjd);
    key = (year(dates) - 1993) * 12 + month(dates) - 1;
    [monthlyKey, ~, group] = unique(key);
    monthly = accumarray(group, height, [], @mean);
    for m = 0:11
        rows = mod(monthlyKey, 12) == m;
        monthly(rows) = monthly(rows) - mean(monthly(rows));
    end
    [shared, ia, ig] = intersect(round((A.t - 1993) * 12), monthlyKey);
    ok = isfinite(A.ALT_TG(ia, j)) & isfinite(monthly(ig));
    shared = shared(ok);
    ia = ia(ok);
    ig = ig(ok);
    S.t = 1993 + shared / 12;
    y = [monthly(ig), A.ALT_TG(ia, j)];
    y = y - mean(y);
    S.series = y;
    S.rates = zeros(2, 1);
    S.se = S.rates;
    S.n = numel(shared);
    for k = 1:2
        F = trend_hac(S.t, y(:, k), 24, 30);
        S.rates(k) = F.rate;
        S.se(k) = F.se;
    end
    S.plot_t = (1993 + (min(shared):max(shared)) / 12)';
    S.plot_y = nan(numel(S.plot_t), 2);
    [~, row] = ismember(shared, (min(shared):max(shared))');
    S.plot_y(row, :) = y;
    ok = key >= min(shared) & key <= max(shared);
    d = dates(ok);
    S.daily_t = year(d) + days(d - datetime(year(d), 1, 1)) ./ days(datetime(year(d) + 1, 1, 1) - datetime(year(d), 1, 1));
    S.daily_height = height(ok) - mean(height(ok));
    cells = O.selected_indices_matrix(:, j);
    cells = cells(isfinite(cells));
    [ix, iy] = ind2sub([numel(O.lon), numel(O.lat)], cells);
    S.alt_cells = nan(numel(A.t), numel(cells));
    G = matfile(fullfile(root, 'outputs', 'grids.mat'));
    for k = 1:numel(cells)
        S.alt_cells(:, k) = double(reshape(G.ALTds(ix(k), iy(k), :), [], 1));
    end
    S.station = j;
    S.longitude = O.lon(ix);
    S.latitude = O.lat(iy);
    S.r = O.top_correlations_matrix(1:numel(cells), j);
end

function key = station_key(name)
    key = upper(string(name));
    key = erase(key, ["'", "."]);
    key = regexprep(key, '\s*\([^)]*\)', '');
    key = regexprep(key, ',.*$', '');
    key = strtrim(regexprep(key, '\s+', ' '));
end

function d = distance_km(a, b)
    % Great-circle distance; a is one lon/lat point, b contains rows of points.
    a = deg2rad(a);
    b = deg2rad(b);
    h = sin((b(:, 2) - a(2)) / 2) .^ 2 + cos(a(2)) * cos(b(:, 2)) .* sin((b(:, 1) - a(1)) / 2) .^ 2;
    d = 2 * 6371.0088 * asin(sqrt(max(0, min(1, h))));
    d(~isfinite(h)) = NaN;
end
