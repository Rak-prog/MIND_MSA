function perm = schaefer_perm_group_lr_by_network(roi_labels, mode)
% roi_labels: string/cellstr (Nx1)
% mode: "LH_then_RH" (default) or "interleave_LR"

if nargin < 2
    mode = "LH_then_RH";
end

roi_labels = string(roi_labels(:));

% Parse hemisphere + network + index from labels like:
% rh_7Networks_RH_DorsAttn_Post_4
hemi = strings(size(roi_labels));
net7 = strings(size(roi_labels));
idxN = zeros(size(roi_labels));

for i = 1:numel(roi_labels)
    s = roi_labels(i);

    % hemi from prefix
    if startsWith(s,"lh_"), hemi(i) = "LH";
    elseif startsWith(s,"rh_"), hemi(i) = "RH";
    else, hemi(i) = "UNK";
    end

    % Extract the part after "..._RH_" or "..._LH_"
    % Example tail: "DorsAttn_Post_4"
    tok = regexp(s, "7Networks_(LH|RH)_(.*)$", "tokens", "once");
    if isempty(tok)
        net7(i) = "UNK";
        idxN(i) = i;
        continue;
    end
    tail = string(tok{2});  % e.g., "DorsAttn_Post_4"

    % numeric index at end
    tokN = regexp(tail, "(\d+)$", "tokens", "once");
    if isempty(tokN)
        idxN(i) = i;
    else
        idxN(i) = str2double(tokN{1});
    end

    % remove trailing "_<num>"
    tailNoNum = regexprep(tail, "_\d+$", "");  % e.g., "DorsAttn_Post"

    % Map to 7-network bucket (ignore subnetwork suffixes)
    if startsWith(tailNoNum,"Vis")
        net7(i) = "Vis";
    elseif startsWith(tailNoNum,"SomMot")
        net7(i) = "SomMot";
    elseif startsWith(tailNoNum,"DorsAttn")
        net7(i) = "DorsAttn";
    elseif startsWith(tailNoNum,"SalVentAttn") || startsWith(tailNoNum,"VentAttn") || startsWith(tailNoNum,"Sal")
        net7(i) = "SalVentAttn";
    elseif startsWith(tailNoNum,"Limbic")
        net7(i) = "Limbic";
    elseif startsWith(tailNoNum,"Cont")
        net7(i) = "Cont";
    elseif startsWith(tailNoNum,"Default")
        net7(i) = "Default";
    else
        net7(i) = "UNK";
    end
end

% Define desired 7-network order
order = ["Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default","UNK"];

netRank = zeros(size(net7));
for k=1:numel(order)
    netRank(net7==order(k)) = k;
end
netRank(netRank==0) = numel(order);

% Within-network ordering:
% primarily by idxN, secondarily by original order (stable)
orig = (1:numel(roi_labels))';

switch string(mode)
    case "LH_then_RH"
        hemiRank = ones(size(hemi));
        hemiRank(hemi=="RH") = 2;   % LH=1, RH=2
    
        % network -> hemisphere -> index -> original
        [~, perm] = sortrows([netRank, hemiRank, idxN, orig], [1 2 3 4]);

    case "interleave_LR"
        % Interleave LH/RH for same idxN: (LH idx1, RH idx1, LH idx2, RH idx2)
        % Make hemi rank LH=1 RH=2 but put it after idx
        hemiRank = ones(size(hemi)); hemiRank(hemi=="RH") = 2;
        [~, perm] = sortrows([netRank, idxN, hemiRank, orig], [1 2 3 4]);

        % The sortrows already interleaves if both hemis share idxN.
        % This assumes idxN numbering matches across hemispheres (it usually does).
    otherwise
        error('Unknown mode: %s', mode);
end

end
