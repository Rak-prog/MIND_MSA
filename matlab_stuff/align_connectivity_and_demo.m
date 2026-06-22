function [conn2, demo2, audit] = align_connectivity_and_demo(conn, demo, subcortex, cfg)

% --- validations ---
assert(istable(demo), "demo must be a table");
assert(istable(subcortex), "subcortex must be a table");
assert(any(strcmp(demo.Properties.VariableNames, cfg.demo_id_col)), ...
    "Demo ID column '%s' not found.", cfg.demo_id_col);

% allow different id col name for subcortex, otherwise reuse demo_id_col
if isfield(cfg, "subcortex_id_col") && ~isempty(cfg.subcortex_id_col)
    sub_id_col = cfg.subcortex_id_col;
else
    sub_id_col = cfg.demo_id_col;
end
assert(any(strcmp(subcortex.Properties.VariableNames, sub_id_col)), ...
    "Subcortex ID column '%s' not found.", sub_id_col);

% suffix for conflicting subcortex variable names
if isfield(cfg, "subcortex_suffix") && ~isempty(cfg.subcortex_suffix)
    sub_suffix = string(cfg.subcortex_suffix);
else
    sub_suffix = "_sub";
end

% --- pull and normalize IDs ---
demo_ids_raw = string(demo.(cfg.demo_id_col));
conn_ids_raw = string(conn.ids);
sub_ids_raw  = string(subcortex.(sub_id_col));

demo_ids = normalize_ids(demo_ids_raw, cfg);
conn_ids = normalize_ids(conn_ids_raw, cfg);
sub_ids  = normalize_ids(sub_ids_raw,  cfg);

% --- optional removals ---
remove_ids = string(cfg.remove_ids(:));
remove_ids = normalize_ids(remove_ids, cfg);

rm_conn = ismember(conn_ids, remove_ids);
rm_demo = ismember(demo_ids, remove_ids);
rm_sub  = ismember(sub_ids,  remove_ids);

conn_keep_mask = ~rm_conn;
demo_keep_mask = ~rm_demo;
sub_keep_mask  = ~rm_sub;

conn_ids2 = conn_ids(conn_keep_mask);
demo_ids2 = demo_ids(demo_keep_mask);
sub_ids2  = sub_ids(sub_keep_mask);

mat2 = conn.matrices(:,:,conn_keep_mask);
grp2 = conn.group(conn_keep_mask);
demo2 = demo(demo_keep_mask, :);
sub2  = subcortex(sub_keep_mask, :);

% --- intersect and reorder to common order across ALL three ---
% step 1: conn ∩ demo
[ids_cd, ia_conn, ia_demo] = intersect(conn_ids2, demo_ids2, 'stable');
mat2  = mat2(:,:,ia_conn);
grp2  = grp2(ia_conn);
demo2 = demo2(ia_demo, :);

% step 2: (conn∩demo) ∩ subcortex
[common_ids, ia_cd, ia_sub] = intersect(ids_cd, sub_ids2, 'stable');
mat2  = mat2(:,:,ia_cd);
grp2  = grp2(ia_cd);
demo2 = demo2(ia_cd, :);     % already ordered like ids_cd
sub2  = sub2(ia_sub, :);

% --- force ID columns to normalized common IDs ---
demo2.(cfg.demo_id_col) = common_ids;
sub2.(sub_id_col)       = common_ids;

% --- append subcortex columns into demo2 (excluding ID) ---
vars_to_add = setdiff(string(sub2.Properties.VariableNames), string(sub_id_col), 'stable');

% handle name collisions: if demo already has a var with same name, rename the incoming one
demo_vars = string(demo2.Properties.VariableNames);
new_names = vars_to_add;

collide = ismember(new_names, demo_vars);
new_names(collide) = new_names(collide) + sub_suffix;

% apply renaming only for the selected columns
sub_add = sub2(:, cellstr(vars_to_add));
sub_add.Properties.VariableNames = cellstr(new_names);

% horizontally concatenate (same row order guaranteed)
demo2 = [demo2 sub_add];

% --- assertions (fail fast) ---
assert(size(mat2,3) == height(demo2), ...
    "Mismatch: matrices N=%d, demo N=%d", size(mat2,3), height(demo2));
assert(all(common_ids == normalize_ids(string(demo2.(cfg.demo_id_col)), cfg)), ...
    "ID mismatch after alignment");

% --- audit report for debugging ---
missing_in_demo = setdiff(conn_ids2, demo_ids2, 'stable');
missing_in_conn = setdiff(demo_ids2, conn_ids2, 'stable');
missing_in_sub_from_cd = setdiff(ids_cd, sub_ids2, 'stable');
missing_in_cd_from_sub = setdiff(sub_ids2, ids_cd, 'stable');

audit = struct();
audit.report = table();
audit.report.common_N = numel(common_ids);
audit.report.removed_from_conn = sum(rm_conn);
audit.report.removed_from_demo = sum(rm_demo);
audit.report.removed_from_sub  = sum(rm_sub);

audit.report.missing_in_demo = {missing_in_demo};
audit.report.missing_in_conn = {missing_in_conn};
audit.report.missing_in_subcortex = {missing_in_sub_from_cd};
audit.report.missing_in_conn_or_demo = {missing_in_cd_from_sub};

audit.added_subcortex_columns = cellstr(new_names);
audit.colliding_subcortex_columns = cellstr(vars_to_add(collide));

% final output struct
conn2 = conn;
conn2.matrices = mat2;
conn2.ids      = common_ids;
conn2.group    = grp2;

end

