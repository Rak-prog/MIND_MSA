function [TOT, HC_cells, PD_cells, PSP_cells, MSA_cells, HCyoung_cells] = loading_groups(ids, matrices, grouping)

HC_cells  = {};
PD_cells  = {};
PSP_cells = {};
MSA_cells = {};
HCyoung_cells = {};

for i = 1:numel(ids)

    % matrix (nROI x nROI)
    M = matrices(:,:,i);

    % IMPORTANT: take the i-th group (scalar)
    gi = grouping(i);

    % normalize group label robustly
    g2 = upper(strtrim(string(gi)));   
    switch g2
        case "HC"
            HC_cells{end+1} = M;
        case "PD"
            PD_cells{end+1} = M;
        case "PSP"
            PSP_cells{end+1} = M;
        case "MSA"
            MSA_cells{end+1} = M;
        case "HCYOUNG"
            HCyoung_cells{end+1} = M;
        otherwise
            warning('Unknown group "%s" for subject %s (skipping).', g2, string(ids(i)));
    end
end

% Convert to 3D (handle empty groups safely)
HC  = cat(3, HC_cells{:});
PD  = cat(3, PD_cells{:});
PSP = cat(3, PSP_cells{:});
MSA = cat(3, MSA_cells{:});
HCyoung = cat(3, HCyoung_cells{:});

% MUST BE IN THE SAME ORDER OF THE F TEST
TOT = cat(3, HC, PD, PSP, MSA, HCyoung);

end


