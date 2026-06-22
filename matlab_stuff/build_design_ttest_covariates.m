function [X, c] = build_design_ttest_covariates(demo_sub, maskA_sub, maskB_sub, cfg, use_protocol)

% maskA_sub/maskB_sub are within the subset already, so they must sum to 1
assert(all(maskA_sub + maskB_sub == 1), "Subjects not exclusively in A or B.");

% group columns
gA = double(maskA_sub(:));
gB = double(maskB_sub(:));

% covariates (edit column names if needed)
age = demo_sub.age_final(:);
sex = code_sex_binary(demo_sub.sex_final);
icv = demo_sub.ICV;
protocol = demo_sub.PROTOCOL;

% demean covariates
age = age - mean(age,'omitnan');
icv = icv - mean(icv,'omitnan');

if use_protocol == false
    X = [gA, gB, age, sex, icv];
    X_vifs = [gA, age, sex, icv];
elseif use_protocol == true
    X = [gA, gB, age, sex, icv, protocol];  % addo protocol as well to se eif there is an effect of the protcol change on the results
    X_vifs = [gA, age, sex, icv, protocol];
end

% firts two columns are collinear G2 = 1-G1 this is needed to  do GLM ttest
% However, I drop G2 to caluclate vifs otherwise I get the warning on rank
% deficiency. In this case G1 and G2 are encoded as 0 and 1 by only the
% column 1
if any(vifs_mat(X_vifs) >= 5)
    error("there is a vif value in this dsign above 5 check please")
end

if use_protocol == false
    c = [1 -1 0 0 0];
elseif use_protocol == true
    c = [1 -1 0 0 0 0]; % added a zero for the protocol effect
end
end

