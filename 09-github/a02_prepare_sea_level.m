%% Sea level and offshore composites
% M: monthly TGs (mm), columns follow N and L=[lon lat]; t: decimal years.
% ALT: SLA+DAC composite for each TG, in the same order and units as M.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'functions'));
D = load(fullfile(root, 'data', 'sea_level.mat'));
M = D.M;
L = D.L;
N = string(D.N(:));
t = D.t(:);
% Remove calendar-month means, keeping the original NaN rows.
Mds = double(M);
Mdtds = double(detrend(M, 'omitnan', 'SamplePoints', t));
for m = 1:12
    Mds(m:12:end, :) = Mds(m:12:end, :) - mean(Mds(m:12:end, :), 1, 'omitnan');
    Mdtds(m:12:end, :) = Mdtds(m:12:end, :) - mean(Mdtds(m:12:end, :), 1, 'omitnan');
end

%% Gridded SLA+DAC: remove seasonal cycles and correlate with the TGs
input = fullfile(root, 'data', 'altimetry_monthly.h5');
D.lon = h5read(input, '/lon');
D.lat = h5read(input, '/lat');
nx = numel(D.lon);
ny = numel(D.lat);
nt = numel(t);
G = matfile(fullfile(root, 'outputs', 'grids.mat'), 'Writable', true);
G.ALTds(nx, ny, nt) = single(0);
G.R_grid(nx, ny, 131) = single(0);
for first = 1:24:ny
    count = min(24, ny - first + 1);
    A = double(h5read(input, '/ALT', [1 first 1], [nx count nt]));
    A = reshape(A, [], nt)';
    ds = A;
    usable = sum(isfinite(A), 1) >= 3; % Exclude empty or insufficient grid records.
    detrended = detrend(A(:, usable), 'omitnan', 'SamplePoints', t);
    for m = 1:12
        ds(m:12:end, :) = ds(m:12:end, :) - mean(ds(m:12:end, :), 1, 'omitnan');
        detrended(m:12:end, :) = detrended(m:12:end, :) - mean(detrended(m:12:end, :), 1, 'omitnan');
    end
    ds = single(ds);
    dt = nan(size(A), 'single');
    dt(:, usable) = single(detrended);
    r = single(corr(double(dt), Mdtds, 'Rows', 'pairwise'));
    G.ALTds(:, first:first + count - 1, :) = reshape(ds', nx, count, nt);
    G.R_grid(:, first:first + count - 1, :) = reshape(r, nx, count, 131);
end

%% Select up to 12 offshore cells and form correlation-weighted ALT
% Delaware TGs share Cape May cells/weights. St. Lawrence and Berwick use
% their original geographic windows and correlation thresholds.
[lon, lat] = ndgrid(D.lon, D.lat);
lon = lon(:);
lat = lat(:);
names = N;
index = @(name)find(names == name);
rim = index('RIMOUSKI');
stf = index('ST-FRANCOIS');
cape = index('CAPE MAY');
newbold = index('NEWBOLD');
berwick = index('BERWICK');
selected = nan(12, 131);
weights = selected;
actual = selected;
distance = selected;
ALT = nan(nt, 131);
for j = 1:131
    r = double(reshape(G.R_grid(:, :, j), [], 1));
    if j > cape && j <= newbold
        cells = sharedCells;
        w = sharedWeights;
    else
        box = true(size(r));
        threshold = .5;
        if j == berwick
            box = abs(lon - L(j, 1)) <= 2 & abs(lat - L(j, 2)) <= 2;
            threshold = 0;
        elseif j == cape
            box = lon >= L(j, 1) - 1 & lon <= L(j, 1) + 2 & abs(lat - L(j, 2)) <= 2;
            threshold = 0;
        elseif j <= rim
            box = lon >= L(rim, 1) & lon <= L(rim, 1) + 5 & abs(lat - L(rim, 2)) <= 2.5;
            threshold = .4;
            if j <= stf
                threshold = 0;
            end
        end
        candidates = find(box & isfinite(r) & r >= threshold);
        [~, order] = sort(r(candidates), 'descend');
        cells = candidates(order(1:min(12, numel(order))));
        w = r(cells);
        if j == cape
            sharedCells = cells;
            sharedWeights = w;
        end
    end
    assert(~isempty(cells) && abs(sum(w)) > 1e-8, 'No valid offshore ALT.');
    n = numel(cells);
    selected(1:n, j) = cells;
    weights(1:n, j) = w;
    actual(1:n, j) = r(cells);
    distance(1:n, j) = distance_km(L(j, :), [lon(cells), lat(cells)]);
    [ix, iy] = ind2sub([nx ny], cells);
    series = nan(nt, n);
    for k = 1:n
        series(:, k) = double(reshape(G.ALTds(ix(k), iy(k), :), nt, 1));
    end
    good = isfinite(series);
    series(~good) = 0;
    ALT(:, j) = (series * w) ./ (double(good) * w);
end
D.M = M;
D.L = L;
D.N = N;
D.t = t;
D.Mds = Mds;
D.Mdtds = Mdtds;
D.ALT = ALT;
D.ALT_TG = ALT - Mds;
D.ALT_TGdt = detrend(D.ALT_TG, 'omitnan', 'SamplePoints', t);
D.selected_indices_matrix = selected;
D.top_correlations_matrix = weights;
D.actual_correlations_matrix = actual;
D.selected_distance_km_matrix = distance;
O = D;
save(fullfile(root, 'outputs', 'results.mat'), 'O', '-append');
fprintf('Sea level: 131 TGs, 1993-2024; correlations and composites rebuilt.\n');

%% Distance calculation
function d = distance_km(a, b)
    % Great-circle distance; a is one lon/lat point, b contains rows of points.
    a = deg2rad(a);
    b = deg2rad(b);
    h = sin((b(:, 2) - a(2)) / 2) .^ 2 + cos(a(2)) * cos(b(:, 2)) .* sin((b(:, 1) - a(1)) / 2) .^ 2;
    d = 2 * 6371.0088 * asin(sqrt(max(0, min(1, h))));
    d(~isfinite(h)) = NaN;
end
