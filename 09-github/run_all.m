function run_all(makeFigures)
    % Rebuild from bundled monthly inputs; run_all(false) skips the figures.
    if nargin == 0
        makeFigures = true;
    end
    root = fileparts(mfilename('fullpath'));
    addpath(root, fullfile(root, 'functions'), fullfile(root, 'external', 'm_map'));
    tic
    a01_river_indices
    a02_prepare_sea_level
    a03_fit_models
    a04_results_and_tables
    if makeFigures
        main_figures;
        supplementary_figures;
    end
    fprintf('Completed in %.1f minutes.\n', toc / 60);
end
