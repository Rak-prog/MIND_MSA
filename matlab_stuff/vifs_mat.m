function v = vifs_mat(X)
    % Remove intercept if first column is all ones
    if ~isempty(X) && all(X(:,1) == 1)
        X = X(:,2:end);
    end

    p = size(X,2);
    v = NaN(p,1);

    for i = 1:p
        y = X(:,i);
        Xi = X(:, setdiff(1:p, i));

        % Keep only rows that are complete for y and all Xi
        rows = all(~isnan([y Xi]), 2);
        y = y(rows);
        Xi = Xi(rows,:);

        % Need enough samples
        if numel(y) < 3 || isempty(Xi)
            continue;
        end

        % Skip constant columns
        if var(y, 'omitnan') == 0
            continue;
        end

        % Add intercept for regression
        Xi = [ones(size(Xi,1),1) Xi];

        % Solve least squares (more stable than inv)
        b = Xi \ y;
        yhat = Xi * b;

        % R^2
        ss_res = sum((y - yhat).^2);
        ss_tot = sum((y - mean(y)).^2);

        if ss_tot == 0
            continue;
        end

        R2 = 1 - ss_res/ss_tot;

        % Guard against numerical issues
        if R2 >= 1
            v(i) = Inf;
        else
            v(i) = 1 / (1 - R2);
        end
    end
end
