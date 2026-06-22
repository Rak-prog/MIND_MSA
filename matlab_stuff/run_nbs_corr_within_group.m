function res = run_nbs_corr_within_group(conn2, demo2, maskG, varName, gName, cfg, signDir, use_protocol)
% signDir = +1 or -1
% Design: [1 age sex clinical]
% Contrast: [0 0 0 signDir]
%
% Requires: cfg.node_coor, cfg.node_label, cfg.path_out, cfg.TH/cfg.ALPHA/cfg.NPERM (strings)

if ~(signDir==1 || signDir==-1)
    error('signDir must be +1 or -1');
end

varName = string(varName);
gName   = string(gName);

% --- Build a valid mask: no NaNs in covariates/regressor ---
age0 = demo2.age_final;
sex0 = demo2.sex_final;
icv0 = demo2.ICV;
clin0 = demo2.(varName);
protocol0 = demo2.PROTOCOL;

maskValid = maskG & ~isnan(age0) & ~isnan(sex0) & ~isnan(icv0) & ~isnan(clin0) & ~isnan(protocol0);

n = sum(maskValid);
if n < 8
    error("Too few subjects for within-group NBS (%s, %s): n=%d", gName, varName, n);
end

% get only the demgraphic inside the group under examination 
age  = age0(maskValid);
sex  = sex0(maskValid);
icv  = icv0(maskValid);
clin = clin0(maskValid);
prot = protocol0(maskValid);

% demeaning only inside the group for the within group comparison
age  = age  - mean(age,  'omitnan');
icv  = icv  - mean(icv,  'omitnan');
clin = clin - mean(clin, 'omitnan');

% --- Subset matrices and build GLM y (upper triangle edges) ---
A = conn2.matrices(:,:,maskValid);  % nROI x nROI x n
N = size(A,1);

ind_upper = find(triu(ones(N,N),1));
Y = zeros(n, numel(ind_upper));
for s = 1:n
    M = A(:,:,s);
    Y(s,:) = M(ind_upper);
end

% --- Design matrix: [G1 age sex clinical] ---
if use_protocol == false
    X = [ones(n,1), age, sex, icv, clin];
elseif use_protocol == true
    X = [ones(n,1), age, sex, icv, prot, clin];
end

if any(vifs_mat(X) >= 5)
    error("there is a vif value in this dsign above 5 check please")
end

% --- Save files for NBS ---
tag = sprintf("%s_corr_%s_%s", gName, varName, ternary(signDir==1,"pos","neg"));

outdir = cfg.path_out;
if ~exist(outdir,'dir'); mkdir(outdir); end

matFile = fullfile(outdir, tag + ".mat");
designFile   = fullfile(outdir, tag + "_design.txt");
contrastFile = fullfile(outdir, tag + "_contrast.txt");

% NBS expects variable "matrices" inside .mat: (N x N x n)
matrices = A; %#ok<NASGU>
save(matFile, "matrices");

dlmwrite(designFile, X, 'delimiter','\t');
if use_protocol == false
    dlmwrite(contrastFile, [0 0 0 0 signDir], 'delimiter','\t');
elseif use_protocol == true 
    dlmwrite(contrastFile, [0 0 0 0 0 signDir], 'delimiter','\t'); % add also protocol here as covariate
end
% --- NBS UI ---
UI = struct();
UI.method.ui   = 'Run NBS';
UI.test.ui     = 't-test';   % t-test on the regressor via contrast
UI.design.ui   = char(designFile);
UI.contrast.ui = char(contrastFile);
UI.thresh.ui   = cfg.TH;
UI.matrices.ui = char(matFile);

UI.node_coor.ui  = cfg.node_coor;
UI.node_label.ui = cfg.node_label;

UI.perms.ui  = cfg.NPERM;
UI.alpha.ui  = cfg.ALPHA;
UI.size.ui   = 'Extent';
UI.exchange.ui = '';

res = NBSrun_script(UI);

% store some useful metadata
res.USER.maskValid = maskValid;
res.USER.gName = gName;
res.USER.varName = varName;
res.USER.signDir = signDir;

% save some subjects information for further visualization when g is MSA
if strcmp(gName, 'MSA')
    res.subject_mask = maskG;
    res.subject_idx = find(maskG);
    res.phenotype = demo2.Phenotype(maskG);
end
end

