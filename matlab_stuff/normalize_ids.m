function ids_out = normalize_ids(ids_in, cfg)

ids_out = string(ids_in);
ids_out = strtrim(ids_out);

% unify case
ids_out = lower(ids_out);

% optionally strip prefix/suffix if needed for matching
if isfield(cfg, 'id_normalize_strip_prefix') && ~isempty(cfg.id_normalize_strip_prefix)
    p = lower(string(cfg.id_normalize_strip_prefix));
    ids_out = erase(ids_out, p);
end
if isfield(cfg, 'id_normalize_strip_suffix') && ~isempty(cfg.id_normalize_strip_suffix)
    s = lower(string(cfg.id_normalize_strip_suffix));
    ids_out = erase(ids_out, s);
end

% collapse multiple spaces (rare but can happen if IDs come messy)
ids_out = regexprep(ids_out, '\s+', '');

end
