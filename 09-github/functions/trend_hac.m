function F = trend_hac(t, y, lags, minN)
    % OLS trend; calendar-pair Newey-West errors with Bartlett weights.
    % t is decimal years; y is mm, so the slope and its standard error are mm/yr.
    % Match residuals by calendar lag and apply the n/(n-2) finite-sample adjustment.
    if nargin < 3
        lags = [12 24 36];
    end
    if nargin < 4
        minN = 60;
    end
    ok = isfinite(t) & isfinite(y);
    t = t(ok);
    y = y(ok);
    n = numel(y);
    F = struct('rate', NaN, 'se', nan(1, numel(lags)), 'n', n, 'ordinary_se', NaN);
    if n < max(3, minN)
        return;
    end
    months = round(12 * (t - t(1)));
    assert(max(abs(12 * (t - t(1)) - months)) < .05 && numel(unique(months)) == n);
    X = [ones(n, 1), t - mean(t)];
    b = X \ y;
    e = y - X * b;
    B = inv(X' * X);
    Z = X .* e;
    G = zeros(2, 2, max(lags) + 1);
    G(:, :, 1) = Z' * Z;
    for k = 1:max(lags)
        [found, j] = ismember(months - k, months);
        G(:, :, k + 1) = Z(found, :)' * Z(j(found), :);
    end
    for j = 1:numel(lags)
        S = G(:, :, 1);
        for k = 1:lags(j)
            S = S + (1 - k / (lags(j) + 1)) * (G(:, :, k + 1) + G(:, :, k + 1)');
        end
        V = n / (n - 2) * B * S * B;
        F.se(j) = sqrt(max(V(2, 2), 0));
    end
    F.rate = b(2);
    F.ordinary_se = sqrt((e' * e) / (n - 2) * B(2, 2));
end
