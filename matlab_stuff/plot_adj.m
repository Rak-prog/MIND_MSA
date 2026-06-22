function [degree, strength, effect_mat] = plot_adj(res, cfg, titname, do_within, g1_ind, g2_ind)

sig_edge = full(res.NBS.con_mat{1}); % matrix is sparse make it full
sig_edge = sig_edge + sig_edge';     % restore symmetry (nly upper triangle form before)
if ~isequal(sig_edge, sig_edge')
    error('Matrix is not symmetric')
end

% check diagonal to be rqual to zero
if any(sig_edge(eye(size(sig_edge))~=0))
    error("diagonal element must contain zero")
end
t_full   = full(res.NBS.test_stat);           % NxN test statistic matrix (assumed t), already symmetric
if ~isequal(t_full, t_full')
    error('Matrix is not symmetric')
end
% calculate degree and strength
t_mat    = t_full .* sig_edge; % keep only component edges

degree   = sum(sig_edge, 2); % binary degree in component 
strength = sum(t_mat, 2);    % weighted strength in component 

% calculate effect sizes
if do_within
    df = size(res.GLM.X,1) - rank(res.GLM.X);
    effect_mat = t_full ./ sqrt(t_full.^2 + df);

else
    if isempty(g1_ind) || isempty(g2_ind)
        error('Between-group analysis requires g1_ind and g2_ind.');
    end

    n1 = numel(g1_ind);
    n2 = numel(g2_ind);

    % d from t
    effect_mat = t_full .* sqrt(1/n1 + 1/n2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot adjacency matrices for final
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adj_plot = sig_edge;
%adj_plot(tril(true(size(adj_plot)),-1)) = NaN;
% hide zeros as white too
adj_plot(adj_plot == 0) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE FIGURE WITH SUBPLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Color',[1 1 1], 'Renderer','opengl');

% --- Prepare matrices ---
adj_plot = sig_edge;
adj_plot(adj_plot == 0) = NaN;

t_plot = sig_edge .* t_full;
eff_plot = sig_edge .* effect_mat;

% --- Create tiled layout (better than subplot) ---
t = tiledlayout(fig, 1, 3, 'Padding','compact', 'TileSpacing','compact');

%% =========================
% 1) Binary edges
%% =========================
ax1 = nexttile;
h1 = imagesc(ax1, adj_plot);
axis(ax1, 'square');
box(ax1, 'off');

colormap(ax1, winter);
clim(ax1, [0 1]);
cb1 = colorbar(ax1);

set(h1, 'AlphaData', ~isnan(adj_plot));
set(ax1, 'Color', 'w');
title(ax1, 'Binary edges');

apply_schaefer_grid_and_labels(ax1, cfg, false);

%% =========================
% 2) T-values
%% =========================
ax2 = nexttile;
h2 = imagesc(ax2, t_plot);
axis(ax2, 'square');
box(ax2, 'off');

colormap(ax2, parula);
clim(ax2, [str2double(cfg.TH) 4.5]);
cb2 = colorbar(ax2);

set(h2, 'AlphaData', ~isnan(adj_plot));
set(ax2, 'Color', 'w');
title(ax2, 't-values');

apply_schaefer_grid_and_labels(ax2, cfg, false);

%% =========================
% 3) Effect sizes
%% =========================
ax3 = nexttile;
h3 = imagesc(ax3, eff_plot);
axis(ax3, 'square');
box(ax3, 'off');

colormap(ax3, cool);
clim(ax3, [0.4 0.85]);
cb3 = colorbar(ax3);

set(h3, 'AlphaData', ~isnan(adj_plot));
set(ax3, 'Color', 'w');
title(ax3, 'Effect size');

apply_schaefer_grid_and_labels(ax3, cfg, false);

%% =========================
% FIGURE SIZE
%% =========================
set(fig, 'Units', 'centimeters');
set(fig, 'Position', [5 5 30 10]);  % wide figure

%% =========================
% SAVE
%% =========================
if do_within
    filename = fullfile(cfg.within_out, sprintf('%s_all.png', titname));
    filename_fig = fullfile(cfg.within_out, sprintf('%s_all.fig', titname));
else
    filename = fullfile(cfg.between_out, sprintf('%s_all.png', titname));
    filename_fig = fullfile(cfg.between_out, sprintf('%s_all.fig', titname));
end

exportgraphics(fig, filename, 'Resolution', 300);
savefig(fig, filename_fig);
close(fig);




% return nan for the yabplot toolbox 
degree(degree==0)=NaN;
strength(strength==0)=NaN;


end

function apply_schaefer_grid_and_labels(ax, cfg, showYLabels)

if ~isfield(cfg,'is_schaefer') || ~cfg.is_schaefer
    return;
end
if ~isfield(cfg,'net_bounds') || isempty(cfg.net_bounds)
    return;
end

hold(ax, 'on');
for b = cfg.net_bounds(:)'
    xline(ax, b + 0.5, 'k-', 'LineWidth', 0.8);
    yline(ax, b + 0.5, 'k-', 'LineWidth', 0.8);
end
hold(ax, 'off');

if ~isfield(cfg,'net_names') || ~isfield(cfg,'net_sizes')
    return;
end

net_names = string(cfg.net_names(:));
net_sizes = double(cfg.net_sizes(:));

starts  = cumsum([1; net_sizes(1:end-1)]);
ends    = cumsum(net_sizes);
centers = (starts + ends) / 2;

ax.YTick = centers;
ax.YTickLabel = net_names;
ax.XTick = [];
ax.XTickLabel = [];
%ax.XTickLabelRotation = 45;
ax.FontWeight = 'bold';
ax.FontSize = 11;

end