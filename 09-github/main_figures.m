function figures = main_figures(numbers)
    % Draw all 9 figures, or select one: main_figures(5).
    % Each case below is one figure; plotting helpers follow at the bottom.
    % Open the MATLAB figures for review and editing.
    if nargin == 0
        numbers = 1:9;
    end
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'functions'), fullfile(root, 'external', 'm_map'));
    S = load(fullfile(root, 'outputs', 'results.mat'));
    K = load(fullfile(root, 'data', 'ancillary.mat'));
    cmap = interp1([0 .5 1], [.13 .28 .54; 1 1 1; .68 .18 .10], linspace(0, 1, 256)); % Blue-white-red rate colors.
    figures = gobjects(numel(numbers), 1);
    for figIndex = 1:numel(numbers)
        switch numbers(figIndex)
        case 1 % Coastal signal
            A = S.A;
            O = S.O;
            f = paper_figure('Figure 1', 18, 14);
            ax = axes(f, 'Position', [.06 .60 .90 .34]);
            hold on;
            axis([0 10 0 3]);
            axis off
            patch([0 0 2.5 4 10], [.2 1.5 1.5 .7 .2], [.77 .71 .61], 'EdgeColor', 'none');
            patch([2.5 4 10 10 4 2.5], [1.8 1.45 1.45 .2 .5 .9], [.66 .84 .90], 'EdgeColor', 'none');
            plot([2.5 4 10], [1.8 1.45 1.45], 'Color', [0 .45 .74], 'LineWidth', 1.2);
            plot([2.5 10], [1.45 1.45], 'k--');
            plot([2.6 2.6], [1.1 2], 'k', 'LineWidth', 3);
            plot([7.6 8.1 8.6], [2.3 2.3 2.3], '-s', 'Color', [.3 .3 .3], ...
                'MarkerSize', 6, 'MarkerFaceColor', [.8 .8 .8]);
            plot([8.1 8.1], [2.15 1.5], ':', 'Color', [.4 .4 .4]);
            text(.2, 2.5, 'River discharge');
            quiver(1.7, 2.4, .8, -.5, 0, 'Color', [.85 .33 .1]);
            text(4.6, 2.7, 'Wind-driven coastal setup');
            quiver(5, 2.5, -1.5, -.65, 0, 'Color', [0 .45 .74]);
            text(1, .55, 'Tide gauge');
            text(7, .7, 'Offshore altimetry');
            text(4, .95, 'Coastal signal remaining in ALT-TG', 'FontSize', 7);
            text(ax, 0, 1, 'a', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12);
            G = matfile(fullfile(root, 'outputs', 'grids.mat'));
            bathy = K.bathy;
            names = ["WASHINGTON DC", "NEWBOLD"];
            bounds = {[-79 -73], [36 40]; [-77 -71], [37 41]};
            for k = 1:2
                j = find(A.N == names(k));
                ax = axes(f, 'Position', [.06 + (k - 1) * .49 .12 .40 .39]);
                map_axes(ax, bounds{k, 1}, bounds{k, 2});
                m_pcolor(O.lon, O.lat, double(G.R_grid(:, :, j))');
                shading flat;
                clim([-.3 1]);
                colormap(ax, parula);
                m_contour(bathy.lon_sub, bathy.lat_sub, bathy.Z_sub, [-400 -400], 'Color', [.4 .4 .4], ...
                    'LineWidth', .4);
                m_gshhs('i', 'patch', [.16 .16 .16], 'edgecolor', 'none');
                m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
                delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
                m_scatter(A.L(j, 1), A.L(j, 2), 24, 'w', 'filled', 'MarkerEdgeColor', 'k');
                title(names(k), 'FontSize', 8);
                panel_label(ax, char('a' + k));
            end
            cb = colorbar(ax, 'southoutside');
            cb.Position = [.20 .055 .60 .018];
            cb.Label.String = 'Correlation, r';

        case 2 % Sea-level trends, selected cells, and correlations
            A = S.A;
            O = S.O;
            D = S.D;
            R = S.R;
            G = matfile(fullfile(root, 'outputs', 'grids.mat'));
            nx = numel(O.lon);
            ny = numel(O.lat);
            rate = nan(nx, ny);
            for first = 1:24:ny
                rows = first:min(first + 23, ny);
                X = double(G.ALTds(:, rows, :));
                rate(:, rows) = reshape(grid_slopes(A.t, reshape(X, [], numel(A.t))'), nx, numel(rows));
            end
            tgRate = grid_slopes(A.t, O.Mds);
            cells = unique(O.selected_indices_matrix(:));
            cells = cells(isfinite(cells));
            [ix, iy] = ind2sub([nx ny], cells);
            riverIds = A.q_col(logical(A.usedQ));
            if any(A.usedQsum)
                riverIds = [riverIds; find(R.LRD(:, 1) >= -80 & R.LRD(:, 1) <= -74 & ...
                    R.LRD(:, 2) >= 36 & R.LRD(:, 2) <= 40)];
            end
            riverIds = unique(riverIds);
            f = paper_figure('Figure 2', 18, 15);
            ax = axes(f, 'Position', [.06 .58 .39 .35]);
            map_axes(ax, [-99 -50], [23 53]);
            m_pcolor(O.lon, O.lat, rate');
            shading flat;
            clim([-3 10]);
            colormap(ax, cmap);
            m_gshhs('l', 'patch', [.16 .16 .16], 'edgecolor', 'none');
            m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
            delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
            hq = m_scatter(R.LRD(riverIds, 1), R.LRD(riverIds, 2), 16, 'd', 'filled', ...
                'MarkerFaceColor', [.93 .69 .13]);
            ht = m_scatter(A.L(:, 1), A.L(:, 2), 14, tgRate, 'filled', 'MarkerEdgeColor', 'k');
            legend([hq ht], {'River-index reference gauges', 'Tide gauges'}, 'Location', 'southeast', ...
                'FontSize', 6, 'Box', 'off');
            cb = colorbar(ax, 'Position', [.47 .60 .015 .29]);
            cb.Label.String = 'Trend (mm yr^{-1})';
            panel_label(ax, 'a');
            ax = axes(f, 'Position', [.06 .10 .39 .35]);
            map_axes(ax, [-99 -50], [23 53]);
            m_gshhs('l', 'patch', [.16 .16 .16], 'edgecolor', 'none');
            m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
            delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
            m_scatter(O.lon(ix), O.lat(iy), 7, [.85 .33 .1], 'filled');
            title('Selected ALT cells', 'FontSize', 8);
            panel_label(ax, 'b');
            ax = axes(f, 'Position', [.59 .09 .18 .84]);
            hold on;
            shade_stations(D.groups.estuary, [0 1]);
            r = mean(O.actual_correlations_matrix, 1, 'omitnan')';
            y = (1:131)';
            q = logical(A.usedQ);
            plot(r, y, 'k-', 'LineWidth', .6);
            h0 = plot(r(~q), y(~q), 'ko', 'MarkerSize', 2.5, 'MarkerFaceColor', 'k');
            h1 = plot(r(q), y(q), 'kd', 'MarkerSize', 3.5, 'MarkerFaceColor', [.93 .69 .13]);
            set(ax, 'YDir', 'reverse', 'YAxisLocation', 'right', 'YTick', 1:3:131, 'YTickLabel', A.N(1:3:131), ...
                'FontSize', 6);
            xlim([0 1]);
            ylim([.5 131.5]);
            xlabel('Mean selected-cell r');
            panel_label(ax, 'c');
            legend([h1 h0], {'Assigned discharge', 'No assigned discharge'}, 'Location', 'southwest', ...
                'FontSize', 6, 'Box', 'off');

        case 3 % Network model performance
            A = S.A;
            C = S.C;
            f = paper_figure('Figure 3', 13, 22);
            ax = axes(f, 'Position', [.11 .07 .58 .90]);
            hold on
            shade_stations(logical(A.usedQ), [-.015 1]);
            y = (1:131)';
            hc = barh(y, C.common_adj_r2(:, 3), .88, 'FaceColor', [0 .447 .698], 'EdgeColor', 'none');
            hu = plot(C.common_adj_r2(:, 1), y, 's', 'Color', [.34 .71 .91], 'MarkerFaceColor', [.34 .71 .91], ...
                'MarkerSize', 4);
            hq = plot(C.common_adj_r2(:, 2), y, 'd', 'Color', [.93 .69 .13], 'MarkerFaceColor', [.93 .69 .13], ...
                'MarkerSize', 4);
            set(ax, 'YDir', 'reverse', 'YAxisLocation', 'right', 'YTick', 1:2:131, 'YTickLabel', A.N(1:2:131), ...
                'FontSize', 6);
            xlim([-.015 1]);
            ylim([.25 131.75]);
            xlabel('Adjusted R^2', 'FontSize', 8);
            legend([hc hu hq], {'UV+Q', 'UV', 'Q'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);

        case 4 % Reedy Point composites
            A = S.A;
            W = load(fullfile(root, 'data', 'forcing.mat'));
            j = find(A.N == "REEDY POINT");
            event = -detrend(A.ALT_TG(:, j), 'omitnan', 'SamplePoints', A.t);
            mu = mean(event, 'omitnan');
            sigma = std(event, 'omitnan');
            limits = [mu + sigma, mu - sigma];
            events = {event >= limits(1), event <= limits(2)};
            colors = [.85 .33 .10; 0 .45 .74];
            f = paper_figure('Figure 4', 18, 16);
            ax = axes(f, 'Position', [.08 .71 .88 .25]);
            hold on
            for k = 1:2
                base = limits(k) * ones(size(event));
                y = base;
                y(events{k}) = event(events{k});
                y(~isfinite(event)) = NaN;
                starts = find(isfinite(y) & [true; ~isfinite(y(1:end - 1))]);
                stops = find(isfinite(y) & [~isfinite(y(2:end)); true]);
                for r = 1:numel(starts)
                    ii = starts(r):stops(r);
                    fill([A.t(ii); flipud(A.t(ii))], [y(ii); flipud(base(ii))], colors(k, :), 'EdgeColor', 'none');
                end
            end
            plot(A.t, event, 'k', 'LineWidth', .35);
            hp = yline(limits(1), ':', 'Color', colors(1, :), 'LineWidth', 1);
            hn = yline(limits(2), ':', 'Color', colors(2, :), 'LineWidth', 1);
            legend([hp hn], {'> +1\sigma', '< -1\sigma'}, 'Orientation', 'horizontal', ...
                'Location', 'northeast', 'Box', 'off');
            xlim([1993 2025]);
            xlabel('Year');
            ylabel('TG - ALT anomaly (mm)');
            panel_label(ax, 'a');
            nx = numel(W.lon);
            ny = numel(W.lat);
            fields = {'SLP', 'U10', 'V10'};
            composite = cell(2, 3);
            for c = 1:3
                x = reshape(double(W.(fields{c})), [], numel(A.t))';
                for m = 1:12
                    x(m:12:end, :) = x(m:12:end, :) - mean(x(m:12:end, :), 1, 'omitnan');
                end
                for k = 1:2
                    composite{k, c} = reshape(mean(x(events{k}, :), 1, 'omitnan'), nx, ny);
                end
            end
            [lon, lat] = ndgrid(W.lon, W.lat);
            ix = 1:3:nx;
            iy = 1:3:ny;
            eventNames = ["Positive", "Negative"];
            for k = 1:2
                ax = axes(f, 'Position', [.08 + (k - 1) * .49 .17 .39 .44]);
                map_axes(ax, [-80 -65], [30 45]);
                m_contourf(W.lon, W.lat, composite{k, 1}', -250:25:250, 'LineStyle', 'none');
                clim([-250 250]);
                colormap(ax, cmap);
                m_gshhs('i', 'patch', [.16 .16 .16], 'edgecolor', 'none');
                m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
                delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
                wind_arrows(lon(ix, iy), lat(ix, iy), composite{k, 2}(ix, iy), composite{k, 3}(ix, iy), ...
                    [.92 .92 .92]);
                m_scatter(A.L(j, 1), A.L(j, 2), 24, 'w', 'filled', 'MarkerEdgeColor', 'k');
                m_patch([-71 -67.5 -67.5 -71], [42.8 42.8 44.2 44.2], 'w', 'EdgeColor', 'none');
                wind_arrows(-70.4, 43.2, 2, 0, [0 0 0]);
                m_text(-70.5, 43.7, '2 m s^{-1}', 'FontSize', 7, 'Color', 'k');
                title(sprintf('%s events (n = %d)', eventNames(k), nnz(events{k})), 'FontSize', 8);
                panel_label(ax, char('a' + k));
            end
            cb = colorbar(ax, 'southoutside');
            cb.Position = [.16 .07 .70 .022];
            cb.Label.String = 'Sea-level-pressure anomaly (Pa)';

        case 5 % Estuary river miles
            A = S.A;
            D = S.D;
            C = S.C;
            R = S.R;
            f = paper_figure('Figure 5', 18.5, 16.5);
            groups = {D.groups.delaware, D.groups.stlaw};
            labels = {'Delaware Estuary', 'St. Lawrence Estuary'};
            lonLimits = {[-75.75 -74.5], [-77 -64]};
            latLimits = {[38.7 40.2], [43 56]};
            offsets = [.05 0; 0 .4];
            for k = 1:2
                ids = find(groups{k});
                T = K.river_miles{k};
                if k == 1
                    coordinates = [T.Lon T.Lat];
                    miles = T.RiverMile;
                else
                    coordinates = [T.lon T.lat];
                    miles = T.mimarker;
                end
                valid = all(isfinite(coordinates), 2) & isfinite(miles);
                coordinates = coordinates(valid, :);
                miles = miles(valid);
                mile = zeros(numel(ids), 1);
                for j = 1:numel(ids)
                    [~, nearest] = min(distance_km(A.L(ids(j), :), coordinates));
                    mile(j) = miles(nearest);
                end
                [mile, order] = sort(mile);
                ids = ids(order);
                r = C.common_adj_r2(ids, 1:2);
                r(abs(r) < 1e-8) = 0;
                bottom = .585 - (k - 1) * .48;
                ax = axes(f, 'Position', [.075 bottom .355 .335]);
                hold on
                plot(r(:, 1), mile, '-s', 'Color', [.38 .71 .91], 'MarkerFaceColor', [.38 .71 .91], ...
                    'MarkerSize', 5);
                plot(r(:, 2), mile, '-d', 'Color', [.93 .69 .13], 'MarkerFaceColor', [.93 .69 .13], ...
                    'MarkerSize', 5);
                xlim([-.02 1]);
                ylim([min(0, min(mile) - 5), max(mile) + max(5, .03 * range(mile))]);
                xlabel('Adjusted R^2');
                ylabel('River mile');
                title(labels{k}, 'FontSize', 8);
                legend({'UV', 'Q'}, 'Location', 'best', 'Box', 'off');
                panel_label(ax, char('a' + 2 * (k - 1)));
                ax = axes(f, 'Position', [.525 bottom .340 .335]);
                map_axes(ax, lonLimits{k}, latLimits{k});
                m_gshhs('h', 'patch', [.16 .16 .16], 'edgecolor', 'none');
                m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
                delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
                hu = m_scatter(A.L(ids, 1), A.L(ids, 2), 42, r(:, 1), 's', 'filled', 'MarkerEdgeColor', 'k');
                hq = m_scatter(A.L(ids, 1) + offsets(k, 1), A.L(ids, 2) + offsets(k, 2), 42, r(:, 2), 'd', ...
                    'filled', 'MarkerEdgeColor', 'k');
                handles = [hu hq];
                names = {'UV', 'Q'};
                if k == 2
                    gauges = unique(A.q_col(ids));
                    hg = m_plot(R.LRD(gauges, 1), R.LRD(gauges, 2), 'p', 'Color', [.18 .75 .94], ...
                        'MarkerFaceColor', [.18 .75 .94], 'LineStyle', 'none');
                    handles = [handles hg(1)];
                    names = [names {'River gauge'}];
                end
                colormap(ax, parula);
                clim([0 1]);
                legend(handles, names, 'Location', 'southeast', 'FontSize', 7, 'Box', 'on', 'Color', 'w');
                panel_label(ax, char('b' + 2 * (k - 1)));
            end
            cb = colorbar(ax, 'eastoutside');
            ax.Position = [.525 .105 .340 .335];
            cb.Position = [.90 .175 .018 .675];
            cb.Label.String = 'Adjusted R^2';

        case 6 % Variability
            A = S.A;
            D = S.D;
            f = paper_figure('Figure 6', 18, 13);
            layout = tiledlayout(f, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
            colors = [.65 .34 .16; 0 .51 .70];
            examples = ["PHILADELPHIA", "ATLANTIC CITY"];
            for k = 1:2
                j = find(A.N == examples(k));
                ax = nexttile(layout, 1 + (k - 1) * 3, [1 3]);
                hold on
                plot(A.t, A.ALT_TG(:, j), 'Color', colors(1, :), 'LineWidth', .65);
                plot(A.t, A.corrected(:, j), 'Color', colors(2, :), 'LineWidth', .65);
                text(.97, .94, ...
                    sprintf('\\sigma: %.0f \\rightarrow %.0f mm', D.coherence.sd_raw(j), ...
                    D.coherence.sd_corrected(j)), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 7);
                xlim([1993 2025]);
                xlabel('Year');
                ylabel('ALT-TG (mm)');
                title(examples(k), 'FontSize', 8);
                if k == 1
                    legend({'ALT-TG', 'ALT-TG corrected'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', 7);
                end
                panel_label(ax, char('a' + k - 1));
            end
            groups = {1:15, 16:115, 116:131};
            tiles = [7 8 12];
            spans = [1 4 1];
            limits = [400 100 300];
            for k = 1:3
                ids = groups{k};
                ax = nexttile(layout, tiles(k), [1 spans(k)]);
                hold on
                for j = ids(D.groups.estuary(ids))
                    patch([j - .5 j + .5 j + .5 j - .5], [0 0 limits(k) limits(k)], [.88 .88 .88], ...
                        'EdgeColor', 'none');
                end
                h = bar(ids, [D.coherence.sd_raw(ids), D.coherence.sd_corrected(ids)], 1);
                for j = 1:2
                    h(j).FaceColor = colors(j, :);
                    h(j).EdgeColor = 'none';
                end
                xlim([ids(1) - .75 ids(end) + .75]);
                ylim([0 limits(k)]);
                xticks(round(linspace(ids(1), ids(end), 3 + (k == 2) * 2)));
                if k == 1
                    ylabel('Standard deviation (mm)');
                    panel_label(ax, 'c');
                end
                if k == 2
                    xlabel('Tide-gauge index (coastal order)');
                end
            end

        case 7 % Rates and GNSS
            A = S.A;
            D = S.D;
            f = paper_figure('Figure 7', 18.5, 22);
            g = D.gnss;
            effect = A.corrected_rate - A.raw_rate;
            mapData = [g.rate - A.raw_rate, g.rate - A.corrected_rate, effect];
            mapNames = {'NGL - raw ALT-TG', 'NGL - ALT-TG_{NCP}', 'ALT-TG_{NCP} - raw ALT-TG'};
            for k = 1:3
                bottom = .69 - (k - 1) * .33;
                position = [.06 bottom .40 .28];
                ax = axes(f, 'Position', position);
                map_axes(ax, [-99 -50], [23 53]);
                m_gshhs('l', 'patch', [.16 .16 .16], 'edgecolor', 'none');
                m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
                delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
                mask = true(131, 1);
                limits = [-1 1];
                if k < 3
                    other = g.matched & ~g.good;
                    m_scatter(A.L(other, 1), A.L(other, 2), 15, 'o', 'MarkerEdgeColor', [.55 .55 .55]);
                    mask = g.good;
                    limits = [-5 5];
                end
                m_scatter(A.L(mask, 1), A.L(mask, 2), 22, mapData(mask, k), 'filled', 'MarkerEdgeColor', 'k', ...
                    'LineWidth', .25);
                clim(limits);
                colormap(ax, cmap);
                title(mapNames{k}, 'FontSize', 8);
                panel_label(ax, char('a' + k - 1));
                cb = colorbar(ax, 'eastoutside');
                ax.Position = position;
                cb.Position = [.495 bottom + .025 .014 .23];
                cb.Label.String = 'mm yr^{-1}';
                cb.FontSize = 7;
                labels = string(cb.TickLabels);
                labels(1) = "\leq " + labels(1);
                labels(end) = "\geq " + labels(end);
                cb.TickLabels = cellstr(labels);
            end
            edges = round(linspace(0, 131, 5));
            limits = [-15 15; -5 1; -4 1; -28 8];
            colors = [.85 .20 .16; .55 .55 .55; .93 .62 .10; .05 .48 .72];
            height = (.90 - 3 * .025) / 4;
            for k = 1:4
                ids = (edges(k) + 1:edges(k + 1))';
                y = (1:numel(ids))';
                ax = axes(f, 'Position', [.575 .925 - k * height - (k - 1) * .025 .22 height]);
                hold on
                hb = barh(y, effect(ids), .8, 'FaceColor', [.58 .58 .58], 'EdgeColor', 'none');
                good = g.good(ids);
                other = g.matched(ids) & ~good;
                hg = errorbar(g.rate(ids(good)), y(good) + .28, g.se(ids(good)), 'horizontal', '.', ...
                    'Color', colors(1, :), 'CapSize', 3);
                ho = errorbar(g.rate(ids(other)), y(other) + .28, g.se(ids(other)), 'horizontal', 'o', ...
                    'Color', colors(2, :), 'CapSize', 3, 'MarkerSize', 3);
                hr = errorbar(A.raw_rate(ids), y, A.raw_se(ids), 'horizontal', '.', 'Color', colors(3, :), ...
                    'CapSize', 3);
                hc = errorbar(A.corrected_rate(ids), y - .28, A.corrected_se(ids), 'horizontal', '.', ...
                    'Color', colors(4, :), 'CapSize', 3);
                xline(0, 'k:');
                xlim(limits(k, :));
                ylim([.25 numel(ids) + .75]);
                ticks = 1:3:numel(ids);
                if numel(ids) - ticks(end) >= 2
                    ticks(end + 1) = numel(ids);
                end
                set(ax, 'YDir', 'reverse', 'YAxisLocation', 'right', 'YTick', ticks, ...
                    'YTickLabel', A.N(ids(ticks)), 'FontSize', 5.5, 'LineWidth', .5);
                % Explicit right ticks remain visible in Illustrator EPS imports.
                for ytick = ticks
                    plot(limits(k, 2) - [.018 0] * diff(limits(k, :)), [ytick ytick], 'k', 'LineWidth', .5, ...
                        'HandleVisibility', 'off');
                end
                if k == 1
                    panel_label(ax, 'd');
                    legend([hb hg ho hr hc], ...
                        {'NCP adjustment', 'NGL: good', 'NGL: other', 'Raw ALT-TG', 'Corrected ALT-TG'}, ...
                        'Position', [.575 .945 .37 .045], 'NumColumns', 3, 'FontSize', 5.5, 'Box', 'off');
                end
                if k == 4
                    xlabel('VLM rate (mm yr^{-1})', 'FontSize', 8);
                end
            end

        case 8 % GNSS comparison by discharge assignment
            A = S.A;
            D = S.D;
            f = paper_figure('Figure 8', 18, 9);
            limits = [-6 7];
            for k = 1:2
                mask = D.gnss.good & logical(A.usedQ) == (k == 2);
                x = D.gnss.rate(mask);
                raw = A.raw_rate(mask);
                corrected = A.corrected_rate(mask);
                ax = subplot(1, 2, k, 'Parent', f);
                hold(ax, 'on');
                plot(ax, limits, limits, 'k--', 'LineWidth', 1);
                hr = scatter(ax, x, raw, 30, [.93 .62 .1], 'o', 'filled', 'MarkerEdgeColor', 'k');
                hc = scatter(ax, x, corrected, 25, [.38 .62 .80], 's', 'filled', 'MarkerEdgeColor', 'k');
                lg = legend(ax, [hr hc], ...
                    {sprintf('ALT-TG: RMS = %.2f, r = %.2f', rms(x - raw), corr(x, raw, 'Rows', 'pairwise')), ...
                    sprintf('ALT-TG_{NCP}: RMS = %.2f, r = %.2f', rms(x - corrected), ...
                        corr(x, corrected, 'Rows', 'pairwise'))}, ...
                    'Location', 'southwest', 'Box', 'off', 'Color', 'none', 'FontSize', 6);
                xlim(ax, limits);
                ylim(ax, limits);
                axis(ax, 'square');
                xlabel(ax, 'NGL-derived VLM rate (mm yr^{-1})');
                if k == 1
                    ylabel(ax, 'ALT-TG-derived VLM rate (mm yr^{-1})');
                end
                titles = {'No assigned discharge', 'Assigned discharge'};
                title(ax, titles{k}, 'FontSize', 8);
                panel_label(ax, char('a' + k - 1));
                drawnow;
                lg.Position(2) = lg.Position(2) - .035;
            end

        case 9 % VLM and GIA
            A = S.A;
            D = S.D;
            f = paper_figure('Figure 9', 17, 16);
            rate = A.paired_corrected_rate;
            error = 1.95996398454005 * A.paired_corrected_se;
            colors = cmap;
            limits = [-10 10];
            layout = K.layout;
            positions = [.075 .430 .574 .490; ...
                layout.stLawrencePosition + [.035 0 0 0]; ...
                layout.midAtlanticPosition + [.035 0 0 0]];
            lonLimits = {[-100.5 -52], [-77 -57], [-78.5 -73]};
            latLimits = {[24 51.5], [43 50.8], [36 41]};
            for k = 1:3
                if k == 2
                    positions(k, 2) = mapBox(2) + mapBox(4) - positions(k, 4) - .012;
                end
                if k == 3
                    positions(k, 2) = mapBox(2) + .01;
                end
                ax = axes(f, 'Position', positions(k, :));
                map_axes(ax, lonLimits{k}, latLimits{k});
                m_gshhs('l', 'patch', [.16 .16 .16], 'edgecolor', 'none');
                if k == 1
                    m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
                else
                    m_grid('linestyle', 'none', 'tickdir', 'out', 'xticklabels', [], 'yticklabels', []);
                end
                delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
                m_scatter(A.L(:, 1), A.L(:, 2), 24, rate, 'filled', 'MarkerEdgeColor', [.3 .3 .3], ...
                    'LineWidth', .35);
                clim(limits);
                colormap(ax, colors);
                if k == 1
                    panel_label(ax, 'a');
                    drawnow;
                    mapBox = plot_box(ax);
                    cb = colorbar(ax, 'northoutside');
                    ax.Position = positions(k, :);
                    cb.Position = [mapBox(1) mapBox(2) + mapBox(4) + .028 mapBox(3) .016];
                    cb.Label.String = 'VLM rate (mm yr^{-1})';
                    labels = string(cb.TickLabels);
                    labels(1) = "\leq " + labels(1);
                    labels(end) = "\geq " + labels(end);
                    cb.TickLabels = cellstr(labels);
                end
            end
            for k = 1:2
                if k == 1
                    mask = D.gia.atlantic;
                    coordinate = A.L(:, 2);
                    profile = D.gia.atlantic_profile;
                    ax = axes(f, 'Position', [.755 mapBox(2) .205 mapBox(4)]);
                    hold on
                    fill([profile(:, 2) - profile(:, 3); flipud(profile(:, 2) + profile(:, 3))], ...
                        [profile(:, 1); flipud(profile(:, 1))], [.86 .86 .86], 'EdgeColor', 'none');
                    plot(profile(:, 2), profile(:, 1), 'k', 'LineWidth', 1);
                    errorbar(rate(mask), coordinate(mask), error(mask), 'horizontal', 'LineStyle', 'none', ...
                        'Color', [.6 .6 .6], 'CapSize', 3);
                    scatter(rate(mask), coordinate(mask), 22, rate(mask), 'filled', 'MarkerEdgeColor', [.3 .3 .3]);
                    ylim([24 51.5]);
                    xlabel('VLM rate (mm yr^{-1})');
                    ylabel('Latitude (\circN)');
                    panel_label(ax, 'b');
                else
                    mask = D.gia.gulf;
                    coordinate = -A.L(:, 1);
                    profile = D.gia.gulf_profile;
                    ax = axes(f, 'Position', [mapBox(1) mapBox(2) - .345 mapBox(3) .275]);
                    hold on
                    hg = fill([profile(:, 1); flipud(profile(:, 1))], ...
                        [profile(:, 2) - profile(:, 3); flipud(profile(:, 2) + profile(:, 3))], ...
                        [.86 .86 .86], 'EdgeColor', 'none');
                    hm = plot(profile(:, 1), profile(:, 2), 'k', 'LineWidth', 1);
                    errorbar(coordinate(mask), rate(mask), error(mask), 'LineStyle', 'none', ...
                        'Color', [.6 .6 .6], 'CapSize', 3);
                    hv = scatter(coordinate(mask), rate(mask), 24, rate(mask), 'filled', ...
                        'MarkerEdgeColor', [.3 .3 .3]);
                    xlim([52 100.5]);
                    set(ax, 'XDir', 'reverse');
                    xlabel('Longitude (\circW)');
                    ylabel('VLM rate (mm yr^{-1})');
                    legend([hv hm hg], ...
                        {'ALT-TG_{NCP} (95% CI)', 'ICE-7G ensemble mean', 'ICE-7G \pm1 ensemble SD'}, ...
                        'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
                    panel_label(ax, 'c');
                end
                colormap(ax, colors);
                clim(limits);
            end

        otherwise
            error('Unknown figure number.');
        end
        figures(figIndex) = f;
    end
end


%% Figure helpers
function wind_arrows(lon, lat, u, v, color)
    scale = .30; % degrees latitude per m/s, identical in both panels and the key
    [x, y] = m_ll2xy(lon, lat);
    [xe, yn] = m_ll2xy(lon + scale * u ./ cosd(lat), lat + scale * v);
    quiver(x, y, xe - x, yn - y, 0, 'Color', color, 'LineWidth', .55, 'MaxHeadSize', .8);
end

function position = plot_box(ax)
    % Physical plot box after M_Map's equal-aspect sizing.
    position = ax.Position;
    f = ancestor(ax, 'figure');
    size = f.Position(3:4);
    width = position(3) * size(1);
    height = position(4) * size(2);
    ratio = diff(ax.XLim) / diff(ax.YLim);
    if width / height > ratio
        newWidth = height * ratio / size(1);
        position(1) = position(1) + (position(3) - newWidth) / 2;
        position(3) = newWidth;
    else
        newHeight = width / ratio / size(2);
        position(2) = position(2) + (position(4) - newHeight) / 2;
        position(4) = newHeight;
    end
end

function f = paper_figure(name, width, height)
    f = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 width height], 'PaperUnits', 'centimeters', 'PaperSize', [width height], ...
        'PaperPosition', [0 0 width height], 'Renderer', 'painters', 'InvertHardcopy', 'off', ...
        'DefaultAxesFontSize', 8, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial', ...
        'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k', 'DefaultAxesColor', 'w', ...
        'DefaultTextColor', 'k', 'DefaultAxesBox', 'on', 'DefaultAxesLayer', 'top');
    set(f, 'DefaultLegendTextColor', 'k', 'DefaultLegendColor', 'w', 'DefaultColorbarColor', 'k');
end

function map_axes(ax, lon, lat)
    axes(ax);
    hold on
    m_proj('Miller', 'lon', lon, 'lat', lat);
    m_patch(lon([1 2 2 1]), lat([1 1 2 2]), [.88 .88 .88], 'EdgeColor', 'none');
end

function panel_label(ax, label)
    text(ax, -.05, 1.04, label, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12, 'Clipping', 'off');
end


function shade_stations(mask, limits)
    for j = find(mask(:))'
        patch(limits([1 2 2 1]), [j - .5 j - .5 j + .5 j + .5], [.88 .88 .88], ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function rate = grid_slopes(t, Y, minN)
    % Vectorized OLS slopes with separate finite-month masks for each grid cell.
    if nargin < 3
        minN = 180;
    end
    t = t - mean(t);
    valid = isfinite(Y);
    Y(~valid) = 0;
    tx = double(valid) .* t;
    n = sum(valid);
    denominator = n .* sum(tx .^ 2) - sum(tx) .^ 2;
    rate = (n .* sum(tx .* Y) - sum(tx) .* sum(Y)) ./ denominator;
    rate(n < minN | denominator <= 0) = NaN;
end

function d = distance_km(a, b)
    % Great-circle distance; a is one lon/lat point, b contains rows of points.
    a = deg2rad(a);
    b = deg2rad(b);
    h = sin((b(:, 2) - a(2)) / 2) .^ 2 + cos(a(2)) * cos(b(:, 2)) .* sin((b(:, 1) - a(1)) / 2) .^ 2;
    d = 2 * 6371.0088 * asin(sqrt(max(0, min(1, h))));
    d(~isfinite(h)) = NaN;
end
