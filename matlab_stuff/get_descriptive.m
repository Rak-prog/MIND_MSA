function Summary = get_descriptive(demo, mask_tot)

% ---------- checks ----------
reqVars = {'sex_final','Phenotype','MCI_T0','age_final','Disease_Duration','UMSARSI','COGNITIVE','MOTOR'};
missingVars = reqVars(~ismember(reqVars, demo.Properties.VariableNames));
if ~isempty(missingVars)
    error('demo is missing required variables: %s', strjoin(missingVars, ', '));
end

baseGroups = {'HC','PD','PSP','MSA'};
for i = 1:numel(baseGroups)
    if ~isfield(mask_tot, baseGroups{i})
        error('mask_tot is missing field %s', baseGroups{i});
    end
end

% ---------- ensure formats ----------
if ~iscategorical(demo.sex_final)
    demo.sex_final = categorical(demo.sex_final);
end

phen = string(demo.Phenotype);

% ---------- build masks ----------
masks = struct();

% base groups
masks.HC  = mask_tot.HC;
masks.PD  = mask_tot.PD;
masks.PSP = mask_tot.PSP;
masks.MSA = mask_tot.MSA;

% optional PSP subgroup masks if present
optPSP = {'PSP_sub','PSP_R','PSP_cor'};
for i = 1:numel(optPSP)
    if isfield(mask_tot, optPSP{i})
        masks.(optPSP{i}) = mask_tot.(optPSP{i});
    end
end

% phenotype masks within MSA
masks.MSAC = masks.MSA & phen == "C";
masks.MSAP = masks.MSA & phen == "P";

% MCI masks in whole sample
masks.MCIyes = demo.MCI_T0 == 1;
masks.MCIno  = demo.MCI_T0 == 0;

% cross-combinations within MSA
masks.MSAC_MCIyes = masks.MSAC & demo.MCI_T0 == 1;
masks.MSAC_MCIno  = masks.MSAC & demo.MCI_T0 == 0;
masks.MSAP_MCIyes = masks.MSAP & demo.MCI_T0 == 1;
masks.MSAP_MCIno  = masks.MSAP & demo.MCI_T0 == 0;

% optional: MCI within other diagnostic groups too
masks.MSA_MCIyes = masks.MSA & demo.MCI_T0 == 1;
masks.MSA_MCIno  = masks.MSA & demo.MCI_T0 == 0;

masks.PD_MCIyes  = masks.PD & demo.MCI_T0 == 1;
masks.PD_MCIno   = masks.PD & demo.MCI_T0 == 0;

masks.PSP_MCIyes = masks.PSP & demo.MCI_T0 == 1;
masks.PSP_MCIno  = masks.PSP & demo.MCI_T0 == 0;

% ---------- group order ----------
groupOrder = { ...
    'HC','PD','PSP','MSA', ...
    'MSAC','MSAP', ...
    'MCIyes','MCIno', ...
    'MSA_MCIyes','MSA_MCIno', ...
    'MSAC_MCIyes','MSAC_MCIno', ...
    'MSAP_MCIyes','MSAP_MCIno', ...
    'PD_MCIyes','PD_MCIno', ...
    'PSP_MCIyes','PSP_MCIno', ...
    'PSP_sub','PSP_R','PSP_cor'};

% keep only those actually present
groups = groupOrder(ismember(groupOrder, fieldnames(masks)));
nGroups = numel(groups);

% ---------- sex categories ----------
sexCatsAll = categories(demo.sex_final);
nSexCats = numel(sexCatsAll);

% ---------- helper ----------
numDesc = @(x) struct( ...
    "Mean",   mean(x,'omitnan'), ...
    "SD",     std(x,'omitnan'), ...
    "Median", median(x,'omitnan'), ...
    "IQR",    iqr(x) );

% ---------- variables to summarize ----------
vars = ["age_final","Disease_Duration","UMSARSI","COGNITIVE","MOTOR"];
nVars = numel(vars);

% ---------- preallocate ----------
Group = strings(nGroups,1);
N = zeros(nGroups,1);

Mean   = nan(nGroups,nVars);
SD     = nan(nGroups,nVars);
Median = nan(nGroups,nVars);
IQRv   = nan(nGroups,nVars);

SexCount = zeros(nGroups,nSexCats);
SexPct   = nan(nGroups,nSexCats);

% ---------- loop ----------
for g = 1:nGroups
    groupName = groups{g};
    Group(g) = string(groupName);

    mask = masks.(groupName);
    demo_sub = demo(mask,:);
    N(g) = height(demo_sub);

    for v = 1:nVars
        x = demo_sub.(vars(v));
        d = numDesc(x);
        Mean(g,v)   = d.Mean;
        SD(g,v)     = d.SD;
        Median(g,v) = d.Median;
        IQRv(g,v)   = d.IQR;
    end

    if N(g) > 0
        [cnt, cats] = groupcounts(demo_sub.sex_final);
        for k = 1:numel(cats)
            idxSex = find(strcmp(sexCatsAll, char(cats(k))));
            if ~isempty(idxSex)
                SexCount(g,idxSex) = cnt(k);
            end
        end
        SexPct(g,:) = 100 * SexCount(g,:) / N(g);
    end
end

% ---------- output table ----------
Summary = table(Group, N);

for v = 1:nVars
    vn = vars(v);
    Summary.("Mean_" + vn) = Mean(:,v);
    Summary.("SD_"   + vn) = SD(:,v);
    % Summary.("Median_" + vn) = Median(:,v);
    % Summary.("IQR_"    + vn) = IQRv(:,v);
end

for s = 1:nSexCats
    catName = matlab.lang.makeValidName(sexCatsAll{s});
    Summary.("SexN_"   + catName) = SexCount(:,s);
    Summary.("SexPct_" + catName) = SexPct(:,s);
end

end