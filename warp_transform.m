function [W,A,e] = warp_transform(p, q)
    N = size(p,1);
    px = repmat(p(:,1), 1, N);
    py = repmat(p(:,2), 1, N);
    pz = repmat(p(:,3), 1, N);
    K = sqrt((px - px').^2 + (py - py').^2 + (pz - pz').^2);

    P = [p, ones(N,1)];
    L = [K P; P' zeros(4,4)];
    D = [q - p; zeros(4,3)];
    warning off
    H = L \ D;
    warning on
    if any(isnan(H))
        H = pinv(L) * D;
    end
    W = H(1:N,:);
    A = H(N+1:end, :);
    e = sum(diag(W' * K * W));
end