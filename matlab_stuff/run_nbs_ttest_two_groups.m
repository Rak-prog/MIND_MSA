function res = run_nbs_ttest_two_groups(conn2, demo2, maskA, maskB, g1, g2, cfg, use_protocol)

% maskA/maskB are logical vectors over ALL subjects (length N)
mask = maskA | maskB;

% subset data
Y = conn2.matrices(:,:,mask);     % nROI x nROI x Nsub
d = demo2(mask,:);

% build design + contrast (A > B)

[X, c] = build_design_ttest_covariates(d, maskA(mask), maskB(mask), cfg, use_protocol);

% output files
outdir = fullfile(cfg.path_out, sprintf("ttest_%s_gt_%s", g1, g2));
if ~exist(outdir,'dir'); mkdir(outdir); end

design_file   = fullfile(outdir, sprintf("ttest_%sgt%s.txt", g1, g2));
contrast_file = fullfile(outdir, sprintf("ttest_contrast_%sgt%s.txt", g1, g2));
mat_file      = fullfile(outdir, sprintf("ttest_contrast_%sgt%s.mat", g1, g2));

% write to disk (NBS expects files)
writematrix(X, design_file, 'Delimiter','\t');
writematrix(c, contrast_file, 'Delimiter','\t');

matrices = Y;
save(mat_file, 'matrices', '-v7');

% --- call NBS exactly like your block ---
UI = struct();
UI.method.ui   = 'Run NBS';
UI.test.ui     = 't-test';  % occhi oqua a seconda del test che fai 

UI.design.ui   = design_file;
UI.contrast.ui = contrast_file;
UI.thresh.ui   = cfg.TH;        % MUST be string
UI.matrices.ui = mat_file;

UI.node_coor.ui  = cfg.node_coor;
UI.node_label.ui = cfg.node_label;

UI.perms.ui    = cfg.NPERM;     % MUST be string
UI.alpha.ui    = cfg.ALPHA;     % MUST be string
UI.size.ui     = 'Extent';
UI.exchange.ui = '';

res = NBSrun_script(UI);

end
