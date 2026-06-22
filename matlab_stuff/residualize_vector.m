function r = residualize_vector(y, COV)
    X = [ones(size(COV,1),1) COV];
    mask = ~isnan(y) & all(~isnan(X),2);
    r = nan(size(y));
    if sum(mask) < 4, return; end
    b = X(mask,:) \ y(mask);
    r(mask) = y(mask) - X(mask,:) * b;
end