function T = get_demographic_table(demo, mask_tot)
%GET_DEMOGRAPHIC_TABLE  Demographic summary table.
%
%  Rows    : N | Age | Sex (F) | MCI+ | Disease Duration | UMSARS-I
%  Columns : HC | MSA | MSAC | MSAP | p(HC vs MSA) | p(MSAC vs MSAP)
%
%  Continuous : Shapiro-Wilk (fallback: Lilliefors). Normal -> mean±SD +
%               t-test; Non-normal -> median[IQR] + Mann-Whitney U.
%  Categorical (sex, MCI) : Fisher's exact test.
%  sex_final coding : 1 = Female, 0 = Male.

% ── masks ─────────────────────────────────────────────────────────────────
phen      = string(demo.Phenotype);
mMSA      = mask_tot.MSA;

subs.HC   = demo(mask_tot.HC,        :);
subs.MSA  = demo(mMSA,               :);
subs.MSAC = demo(mMSA & phen == "C", :);
subs.MSAP = demo(mMSA & phen == "P", :);

gNames = {'HC','MSA','MSAC','MSAP'};

% ═════════════════════════════ HELPERS ════════════════════════════════════

    function ok = isNormal(x)
        x = x(~isnan(x));
        % if numel(x) < 4, ok = false; return; end
        % try
        %     [h, ~] = swtest(x);   ok = ~h;
        % catch
            [h, ~] = lillietest(x); ok = ~h;
        %end
    end

    function s = fmtNum(x)
        x = x(~isnan(x));
        if isempty(x), s = "—"; return; end
        if isNormal(x)
            s = sprintf('%.1f ± %.1f', mean(x), std(x));
        else
            s = sprintf('%.1f [%.1f–%.1f]', median(x), quantile(x,.25), quantile(x,.75));
        end
    end

    function p = cmpNum(x1, x2)
        x1 = x1(~isnan(x1)); x2 = x2(~isnan(x2));
        if numel(x1) < 2 || numel(x2) < 2, p = NaN; return; end
        if isNormal(x1) && isNormal(x2)
            [~, p] = ttest2(x1, x2);
        else
            p = ranksum(x1, x2);
        end
    end

    function p = fisherExact(tbl)
        try
            [~, p] = fishertest(tbl);
        catch
            p = NaN;
            warning('fishertest unavailable — requires Statistics and ML Toolbox.');
        end
    end

    function s = fmtP(p)
        if isnan(p),     s = "—";
        elseif p < .001, s = "< 0.001";
        else,            s = sprintf('%.3f', p);
        end
    end

% ═══════════════════════ DESCRIPTIVE COLUMNS ══════════════════════════════

rowLabels = {'N'; 'Age (years)'; 'Sex – F (n, %)'; 'MCI+ (n, %)'; ...
             'Disease Duration (years)'; 'UMSARS-I'};

nRows = numel(rowLabels);
nG    = numel(gNames);
desc  = strings(nRows, nG);

for g = 1:nG
    sub = subs.(gNames{g});
    n   = height(sub);

    desc(1, g) = num2str(n);
    desc(2, g) = fmtNum(sub.age_final);

    % sex: 1=F, 0=M
    nF = sum(sub.sex_final == 1, 'omitnan');
    desc(3, g) = sprintf('%d (%.1f%%)', nF, 100 * nF / n);

    % MCI
    nMCI = sum(sub.MCI_T0 == 1, 'omitnan');
    desc(4, g) = sprintf('%d (%.1f%%)', nMCI, 100 * nMCI / n);

    desc(5, g) = fmtNum(sub.Disease_Duration);   % — for HC (all NaN)
    desc(6, g) = fmtNum(sub.UMSARSI);            % — for HC (all NaN)
end

% ═══════════════════════════ P-VALUE COLUMNS ══════════════════════════════

hc  = subs.HC;   msa = subs.MSA;
mac = subs.MSAC; map = subs.MSAP;

% 2x2 helpers for categorical variables (Fisher)
sexTbl = @(a,b) [sum(a.sex_final==1,'omitnan'), sum(a.sex_final==0,'omitnan'); ...
                 sum(b.sex_final==1,'omitnan'), sum(b.sex_final==0,'omitnan')];

mciTbl = @(a,b) [sum(a.MCI_T0==1,'omitnan'), sum(a.MCI_T0==0,'omitnan'); ...
                 sum(b.MCI_T0==1,'omitnan'), sum(b.MCI_T0==0,'omitnan')];

pHCvMSA  = [ "—"
             fmtP(cmpNum(hc.age_final,          msa.age_final))
             fmtP(fisherExact(sexTbl(hc,  msa)))
             fmtP(fisherExact(mciTbl(hc,  msa)))
             fmtP(cmpNum(hc.Disease_Duration,   msa.Disease_Duration))   % → — (HC=NaN)
             fmtP(cmpNum(hc.UMSARSI,            msa.UMSARSI)) ];         % → — (HC=NaN)

pMACvMAP = [ "—"
             fmtP(cmpNum(mac.age_final,         map.age_final))
             fmtP(fisherExact(sexTbl(mac, map)))
             fmtP(fisherExact(mciTbl(mac, map)))
             fmtP(cmpNum(mac.Disease_Duration,  map.Disease_Duration))
             fmtP(cmpNum(mac.UMSARSI,           map.UMSARSI)) ];

% ═══════════════════════ ASSEMBLE OUTPUT TABLE ═════════════════════════════

T = table(string(rowLabels), desc(:,1), desc(:,2), desc(:,3), desc(:,4), ...
          pHCvMSA, pMACvMAP, ...
          'VariableNames', {'Variable','HC','MSA','MSAC','MSAP', ...
                            'p_HC_vs_MSA','p_MSAC_vs_MSAP'});
end