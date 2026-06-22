function conn = load_connectivity_from_all_minds_mat(mat_path)

S = load(mat_path);
assert(isfield(S, 'all_minds'), "Expected variable 'all_minds' in mat file.");

all_minds = S.all_minds;
ids = string(fieldnames(all_minds));

% infer ROI info from first subject
one = all_minds.(ids{1});
roi_labels = string(one.matrix_col_names);
nROI = size(one.matrix,1);

matrices = zeros(nROI, nROI, numel(ids));
groups   = strings(numel(ids),1);

for i = 1:numel(ids)
    subj = all_minds.(ids{i});
    matrices(:,:,i) = subj.matrix;
    if isfield(subj,'group')
        groups(i) = upper(strtrim(string(subj.group)));
    else
        groups(i) = "";
    end
end

conn = struct();
conn.matrices   = matrices;
conn.ids        = ids(:);
conn.group      = groups(:);
conn.roi_labels = roi_labels(:);

end
