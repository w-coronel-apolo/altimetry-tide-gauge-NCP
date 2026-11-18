function figures = supplementary_figures(numbers)
    % Draw all 6 figures, or select one: supplementary_figures(5).
    % Each case below is one figure; plotting helpers follow at the bottom.
    % Open the MATLAB figures for review and editing.
    if nargin == 0
        numbers = 1:6;
    end
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'functions'), fullfile(root, 'external', 'm_map'));
    S = load(fullfile(root, 'outputs', 'results.mat'));
    K = load(fullfile(root, 'data', 'ancillary.mat'));
    figures = gobjects(numel(numbers), 1);
    for figIndex = 1:numel(numbers)
        switch numbers(figIndex)
        case 1 % Coherence
            A = S.A;
            D = S.D;
            f = paper_figure('Figure S1', 12, 22);
            ax = axes(f, 'Position', [.10 .06 .57 .91]);
            hold on;
            y = (1:131)';
            plot(D.coherence.r_raw, y, '-o', 'Color', [.15 .15 .15], 'MarkerFaceColor', [.15 .15 .15], ...
                'MarkerSize', 2.5, 'LineWidth', .75);
            plot(D.coherence.r_corrected, y, '-o', 'Color', [.8 .15 .15], 'MarkerFaceColor', [.8 .15 .15], ...
                'MarkerSize', 2.5, 'LineWidth', .75);
            set(ax, 'YDir', 'reverse', 'YAxisLocation', 'right', 'YTick', 1:2:131, 'YTickLabel', A.N(1:2:131), ...
                'FontSize', 6);
            xlim([min(0, floor(10 * min([D.coherence.r_raw; D.coherence.r_corrected])) / 10) 1]);
            ylim([.25 131.75]);
            xlabel('Correlation', 'FontSize', 8);
            legend({'TG vs ALT', 'TG vs ALT_{NCP}'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', 7);

        case 2 % Predictor selection
            A = S.A;
            Models = S.Models;
            f = paper_figure('Figure S2', 18, 13);
            layout = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
            ax = nexttile(layout, [1 2]);
            hold on
            plot(1:131, A.candidate_count, 'Color', [.5 .5 .5]);
            plot(1:131, A.selected_count, '-o', 'Color', [0 .4 .7], 'MarkerFaceColor', [0 .4 .7], 'MarkerSize', 2);
            xlim([1 131]);
            ylabel('Number of predictors');
            xlabel('Tide-gauge index (coastal order)');
            legend({'Candidate predictors', 'Selected predictors'}, 'Location', 'northeast', 'Box', 'off');
            panel_label(ax, 'a');
            ax = nexttile(layout);
            hold on;
            edges = -.5:1:max(A.candidate_count) + .5;
            histogram(A.candidate_count, edges, 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none');
            histogram(A.selected_count, edges, 'FaceColor', [0 .45 .74], 'EdgeColor', 'none');
            xlabel('Number of predictors');
            ylabel('Number of tide gauges');
            legend({'Candidate', 'Selected'}, 'Location', 'northwest', 'Box', 'off');
            panel_label(ax, 'b');
            counts = zeros(4, 2);
            families = ["U", "V", "Q", "FC"];
            for j = 1:131
                names = Models.predictors{j}.names;
                selected = names(Models.models{j}.selected);
                for k = 1:4
                    if k <= 2
                        offered = startsWith(names, families(k));
                        chosen = startsWith(selected, families(k));
                    else
                        offered = names == families(k);
                        chosen = selected == families(k);
                    end
                    counts(k, :) = counts(k, :) + [any(offered), any(chosen)];
                end
            end
            ax = nexttile(layout);
            h = bar(counts);
            h(1).FaceColor = [.6 .78 .88];
            h(2).FaceColor = [.93 .39 .12];
            xticklabels({'Zonal wind PCs', 'Meridional wind PCs', 'River discharge', 'Florida Current'});
            xtickangle(25);
            ylabel('Number of tide gauges');
            legend({'Offered', 'Selected'}, 'Box', 'off');
            panel_label(ax, 'c');

        case 3 % Cape May GNSS
            A = S.A;
            O = S.O;
            D = S.D;
            CM = D.cape_may;
            j = CM.station;
            f = paper_figure('Figure S3', 18.5, 13.5);
            ax = axes(f, 'Position', [.08 .49 .76 .46]);
            map_axes(ax, [-78 -72], [37 41]);
            m_gshhs('i', 'patch', [.16 .16 .16], 'edgecolor', 'none');
            m_grid('linestyle', 'none', 'tickdir', 'out', 'fontsize', 7);
            delete(findall(gca, 'Type', 'patch', 'Tag', 'm_grid_color')); % Keep EPS map layers visible.
            hg = m_scatter(CM.longitude, CM.latitude, 30, CM.r, 'filled', 'MarkerEdgeColor', 'k');
            ht = m_scatter(A.L(j, 1), A.L(j, 2), 36, 'w', 'filled', 'MarkerEdgeColor', 'k');
            colormap(ax, parula);
            clim([.9 1]);
            cb = colorbar(ax, 'Position', [.87 .57 .025 .34]);
            cb.Label.String = 'Pearson r';
            legend([ht hg], {'Cape May TG', 'Selected ALT grid'}, 'Location', 'southwest', 'Box', 'off', ...
                'TextColor', 'w', 'FontSize', 8);
            panel_label(ax, 'a');
            ax = axes(f, 'Position', [.12 .70 .31 .19]);
            hold on
            hc = plot(A.t, CM.alt_cells, 'Color', [.85 .33 .1], 'LineWidth', .4);
            ha = plot(A.t, O.ALT(:, j) - mean(O.ALT(:, j), 'omitnan'), 'k', 'LineWidth', .8);
            ht = plot(A.t, O.Mds(:, j) - mean(O.Mds(:, j), 'omitnan'), 'Color', [.18 .75 .94], 'LineWidth', .8);
            xlim([1993 2025]);
            ylabel('Water level (mm)');
            set(ax, 'FontSize', 6);
            legend([hc(1) ha ht], {'ALT grid', 'Composite', 'TG'}, 'Location', 'southwest', 'Box', 'off', ...
                'FontSize', 5);
            ax = axes(f, 'Position', [.13 .10 .73 .27]);
            hold on;
            colors = [.49 .18 .56; .65 .34 .16];
            hr = scatter(CM.daily_t, CM.daily_height, 8, [.86 .79 .89], 'filled');
            h = plot(CM.plot_t, CM.plot_y, 'LineWidth', .8);
            for k = 1:2
                h(k).Color = colors(k, :);
                plot(CM.t, CM.rates(k) * (CM.t - mean(CM.t)), '--', 'Color', colors(k, :), 'LineWidth', 1);
            end
            xlim([min(CM.t) - .1 max(CM.daily_t) + .1]);
            xlabel('Year');
            ylabel('Height (mm)');
            panel_label(ax, 'b');
            values = [CM.daily_height; CM.series(:)];
            pad = max(10, .1 * range(values));
            ylim([min(values) - pad, max(values) + 3 * pad]);
            legend([hr h'], {'Raw GNSS', 'GNSS monthly', 'Raw ALT-TG'}, 'Location', 'northeast', 'Box', 'off', ...
                'FontSize', 7);
            text(.02, .97, ...
                sprintf('GNSS: %.2f +/- %.2f mm yr^{-1}\nALT-TG: %.2f +/- %.2f mm yr^{-1}\n+/- 1 HAC24 SE', ...
                    CM.rates(1), CM.se(1), CM.rates(2), CM.se(2)), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, 'BackgroundColor', 'w');

        case 4 % Component time series
            A = S.A;
            C = S.C;
            names = ["TROIS-RIVIERES", "NEUVILLE", "QUEBEC", "NEW LONDON", "KINGS POINT", "REEDY POINT", ...
                "PHILADELPHIA", "NEWBOLD", "MONEY POINT", "WASHINGTON DC", "BALTIMORE", "WILMINGTON", ...
                "APALACHICOLA", "BAY WAVELAND YACHT CLUB", "LAKE CHARLES"];
            [found, ids] = ismember(names, A.N);
            assert(all(found));
            rowLimits = [-500 500; -200 200; -500 400; -500 400; -250 250];
            shown = cell(15, 1);
            for k = 1:15
                j = ids(k);
                y = [A.ALT_TG(:, j), C.timeseries_prediction_uv(:, j), C.timeseries_prediction_q(:, j)];
                amplitude = max(abs(y), [], 'all', 'omitnan');
                shown{k} = y + [0 .45 * amplitude -.45 * amplitude];
            end
            % Expand whole rows together, including Newbold and Lake Charles.
            for row = [3 5]
                y = cat(2, shown{(row - 1) * 3 + (1:3)});
                pad = .025 * diff(rowLimits(row, :));
                rowLimits(row, 1) = min(rowLimits(row, 1), 50 * floor((min(y, [], 'all', 'omitnan') - pad) / 50));
                rowLimits(row, 2) = max(rowLimits(row, 2), 50 * ceil((max(y, [], 'all', 'omitnan') + pad) / 50));
            end
            f = paper_figure('Figure S4', 18.5, 20);
            layout = tiledlayout(f, 5, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            colors = [.15 .15 .15; .34 .71 .91; .93 .69 .13];
            for k = 1:15
                row = ceil(k / 3);
                column = mod(k - 1, 3) + 1;
                ax = nexttile(layout);
                hold on
                h = plot(A.t, shown{k}, 'LineWidth', .55);
                for c = 1:3
                    h(c).Color = colors(c, :);
                end
                xlim([1993 2025]);
                ylim(rowLimits(row, :));
                if k == 1
                    ylim([-2000 2000]);
                end
                assert(all(shown{k}(isfinite(shown{k})) >= ax.YLim(1) & ...
                    shown{k}(isfinite(shown{k})) <= ax.YLim(2)));
                title(names(k), 'FontWeight', 'normal', 'FontSize', 7, 'Interpreter', 'none');
                text(.02, .04, sprintf('R^2_{adj}: UV %.2f; Q %.2f', C.timeseries_adj_r2(ids(k), :)), ...
                    'Units', 'normalized', 'FontSize', 6);
                if column == 1
                    ylabel('ALT-TG (mm)');
                end
                if row == 5
                    xlabel('Year');
                else
                    xticklabels([]);
                end
                if k == 1
                    legend(h, {'Observed', 'UV (offset)', 'Q (offset)'}, 'FontSize', 5, 'Box', 'off', ...
                        'Location', 'northwest');
                end
            end

        case 5 % Trend adjustments
            A = S.A;
            D = S.D;
            f = paper_figure('Figure S5', 18.5, 22);
            split = find(A.N == "RIMOUSKI");
            groups = {1:split, split + 1:131};
            colors = [0 .45 .70; .34 .71 .91; .93 .69 .13];
            titles = {'Combined correction', 'Wind-only correction', 'Discharge-only correction'};
            for k = 1:3
                for row = 1:2
                    ids = groups{row};
                    limits = [-10 10];
                    position = [.25 + (k - 1) * .26 .81 .22 .14];
                    if row == 2
                        limits = [-5 5];
                        position = [.25 + (k - 1) * .26 .07 .22 .68];
                    end
                    if k == 2
                        limits = [-1 1];
                    end
                    ax = axes(f, 'Position', position);
                    hold on
                    barh(ids, D.adjustment(ids, k), .8, 'FaceColor', colors(k, :), 'EdgeColor', 'none');
                    xline(0, 'k:');
                    set(ax, 'YDir', 'reverse', 'FontSize', 6);
                    xlim(limits);
                    ylim([ids(1) - .6 ids(end) + .6]);
                    if row == 1
                        ticks = ids;
                        title(titles{k}, 'FontSize', 8);
                        panel_label(ax, char('a' + k - 1));
                    else
                        ticks = unique([ids(1):2:ids(end), ids(end)]);
                        xlabel('Raw - corrected (mm yr^{-1})', 'FontSize', 7);
                    end
                    yticks(ticks);
                    if k == 1
                        yticklabels(A.N(ticks));
                    else
                        yticklabels([]);
                    end
                end
            end

        case 6 % Window uncertainty
            D = S.D;
            W = D.windows;
            f = paper_figure('Figure S6', 18.5, 15);
            colors = [.65 .34 .16; 0 .51 .70];
            names = {'No assigned discharge', 'Assigned discharge'};
            layout = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            for column = 1:2
                group = W.has_q == (column == 2);
                for row = 1:2
                    ax = nexttile(layout, (row - 1) * 2 + column);
                    hold on
                    fields = ["raw_se", "corrected_se"];
                    if row == 2
                        fields = "ratio";
                    end
                    handles = gobjects(numel(fields), 1);
                    for k = 1:numel(fields)
                        values = W.(fields(k));
                        color = colors(k, :);
                        if row == 2
                            color = colors(2, :);
                        end
                        scatter(W.years(group), values(group), 2, .75 + .25 * color, 'filled', ...
                            'MarkerEdgeColor', 'none');
                        summary = nan(19, 3);
                        for years = 10:28
                            summary(years - 9, :) = prctile(values(group & W.years == years), [2.5 50 97.5]);
                        end
                        x = (10:28)';
                        fill([x; flipud(x)], [summary(:, 1); flipud(summary(:, 3))], .85 + .15 * color, ...
                            'EdgeColor', 'none');
                        handles(k) = plot(x, summary(:, 2), 'Color', color, 'LineWidth', 1.2);
                    end
                    uistack(handles, 'top');
                    xlim([9.5 28.5]);
                    panel_label(ax, char('a' + (row - 1) * 2 + column - 1));
                    if row == 1
                        title(names{column}, 'FontSize', 8);
                        ylabel('Trend standard error (mm yr^{-1})');
                        legend(handles, {'Raw ALT-TG', 'Corrected ALT-TG'}, 'Location', 'northeast', ...
                            'Box', 'off', 'FontSize', 7);
                    else
                        yline(1, 'k--');
                        ylim([.2 1.4]);
                        xlabel('Window length from January 1993 (years)');
                        ylabel('Corrected / raw SE');
                    end
                end
            end

        otherwise
            error('Unknown figure number.');
        end
        figures(figIndex) = f;
    end
end


%% Figure helpers
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
