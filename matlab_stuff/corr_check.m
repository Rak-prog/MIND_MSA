function out = corr_check(res, score_col, covar_cols, scorename, cfg, do_color, phenotype)
%CORR_CHECK Correlate clinical regressor with:
% (1) mean of NBS-significant edges
% (2) global mean of all edges
%
% score_col: column of res.GLM.X to use as score

if nargin < 2 || isempty(score_col), score_col = 4; end
if nargin < 6 || isempty(do_color), do_color = false; end
if nargin < 7
    phenotype = [];
end

% --- checks ---
if ~isfield(res,'GLM') || ~isfield(res.GLM,'y') || ~isfield(res.GLM,'X')
    error('res.GLM.y or res.GLM.X missing.');
end
if size(res.GLM.X,2) < score_col
    error('res.GLM.X has < %d columns.', score_col);
end
if ~isfield(res,'NBS') || ~isfield(res.NBS,'con_mat') || isempty(res.NBS.con_mat)
    error('res.NBS.con_mat missing/empty.');
end

score_raw = res.GLM.X(:, score_col);
COV = res.GLM.X(:, covar_cols);

% Residualize all edges
y_res = residualize_matrix(res.GLM.y, COV);
score_res = residualize_vector(score_raw, COV);

% Global mean from residualized similarities
global_mean = mean(y_res, 2, 'omitnan');

% --- map significant edges -> y columns ---
N = res.STATS.N;
ut = triu(true(N),1);
edge_index = zeros(N,N);
edge_index(ut) = 1:nnz(ut);
edge_index = edge_index + edge_index';

C = res.NBS.con_mat{1};
[iSig, jSig] = find(C);
sig_edge_indices = edge_index(sub2ind([N N], iSig, jSig));
sig_edge_indices = sig_edge_indices(sig_edge_indices > 0);

mean_sig = mean(y_res(:, sig_edge_indices), 2, 'omitnan');

% Axis limits
ax1_xlim = [min(mean_sig) max(mean_sig)]; %[min(mean_sig) max(mean_sig) min(score_res) max(score_res)];
ax1_ylim = [min(score_res) max(score_res)];
ax2 = [min(global_mean) max(global_mean) min(score_res) max(score_res)];

% --- plots ---
out = struct();
h = figure('Color',[1 1 1]);

%subplot(1,2,1);
out.sig = corr_panel_ci(mean_sig, score_res, ...
    'Mean MIND Residuals', [scorename , ' Residuals'], [ax1_xlim ax1_ylim], do_color, phenotype);
%subplot(1,2,2);
%out.glob = corr_panel_ci(global_mean, score_res, ...
%    'Residuals Mean (all edges)', ['Residuals ' scorename], ax2, do_color, phenotype);

saveas(h, fullfile(cfg.within_out, sprintf('%s_correlation_final.tiff', scorename)));
end


function s = corr_panel_ci(x, y, xlab, ylab, axis_lims, do_color, phenotype)
% Scatter + fitlm line + 95% CI band + Pearson/Spearman stats

if nargin < 6 || isempty(do_color), do_color = false; end
if nargin < 7, phenotype = []; end

mask = ~isnan(x) & ~isnan(y);

if do_color
    if isempty(phenotype)
        error('do_color is true, but phenotype was not provided.');
    end
    phenotype = string(phenotype);
    phenotype = phenotype(mask);
end

x = x(mask);
y = y(mask);
n = numel(x);

hold on; box off;

if do_color
    idxP = phenotype == "P";
    idxC = phenotype == "C";

    % plot all points first in gray
    scatter(x, y, 28, [0.7 0.7 0.7], 'filled');

    % overlay colored phenotype points
    scatter(x(idxP), y(idxP), 28, [0.85 0.33 0.10], 'filled');
    scatter(x(idxC), y(idxC), 28, [0 0.45 0.74], 'filled');

    legend({'All','P','C'}, 'Location', 'best', 'Box', 'off');
else
    scatter(x, y, 28, 'filled');
end

xlabel(xlab);
ylabel(ylab);

if n < 4
    title(sprintf('n=%d (too small)', n));
    s = struct('n',n,'rP',nan,'pP',nan,'rS',nan,'pS',nan,'r2',nan);
    return;
end

% Pearson + Spearman
[rP,pP] = corr(x, y, 'Type','Pearson',  'Rows','complete');
[rS,pS] = corr(x, y, 'Type','Spearman', 'Rows','complete');

% Linear model for line + CI band
lm = fitlm(x, y);

xx = linspace(min(x), max(x), 200)';
[yy, yCI] = predict(lm, xx, 'Prediction','curve', 'Alpha',0.05);

fill([xx; flipud(xx)], [yCI(:,1); flipud(yCI(:,2))], ...
     [0 0 0], 'FaceAlpha', 0.12, 'EdgeColor','none');

plot(xx, yy, 'k-', 'LineWidth', 1.6);

r2 = lm.Rsquared.Ordinary;
b1 = lm.Coefficients.Estimate(2);
p_b1 = lm.Coefficients.pValue(2);

%title(sprintf(['Pear r=%.2f | p=%.2g\n' ...
%               'Spear r=%.2f | p=%.2g\n' ...
%               'R^2=%.2f | p=%.2g'], ...
%               rP, pP, rS, pS, r2, p_b1));

text(axis_lims(2) - 0.35 * axis_lims(2),axis_lims(4) - 0.9*axis_lims(4),sprintf('R^2=%.2f', r2),'FontWeight','bold'); 

s = struct('n',n, ...
           'rP',rP,'pP',pP, ...
           'rS',rS,'pS',pS, ...
           'r2',r2,'slope',b1,'slope_p',p_b1);

axis([axis_lims]);
set(gca, 'FontWeight', 'Bold','FontSize',12,'FontName', 'Serif'  ) %,'XTick', [] , 'YTick', [])
axis(axis_lims)
xticks(round(linspace(axis_lims(1), axis_lims(2), 5),3))
%yticks(round(linspace(axis_lims(3), axis_lims(4), 5),3))
end