function Yres = residualize_matrix(Y, COV)
    X = [ones(size(COV,1),1) COV];
    mask = all(~isnan(X),2);
    Yres = nan(size(Y));
    if sum(mask) < 4, return; end
    B = X(mask,:) \ Y(mask,:);      % regress all edges at once
    Yhat = X(mask,:) * B;
    Yres(mask,:) = Y(mask,:) - Yhat;
end