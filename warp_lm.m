function rw = warp_lm(r, A, W, p)
    rw = r * A(1:3,1:3);
    rw = bst_bsxfun(@plus, rw, A(4,:));
    np = size(p,1);
    U = sqrt(bst_bsxfun(@minus, repmat(r(:,1),1,np), p(:,1)') .^ 2 + ...
             bst_bsxfun(@minus, repmat(r(:,2),1,np), p(:,2)') .^ 2 + ...
             bst_bsxfun(@minus, repmat(r(:,3),1,np), p(:,3)') .^ 2);
    rw = rw + U * W;   