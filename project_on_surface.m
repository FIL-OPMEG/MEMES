function [Q, dist] = project_on_surface(sSurf, P, center)
    % Center on (0,0,0)
    P = bst_bsxfun(@minus, P, center);
    sSurf.Vertices = bst_bsxfun(@minus, sSurf.Vertices, center);
    % Project points on surface
    Q = 0 .* P;
    for i = 1:length(P)
        proj = tess_ray_intersect(sSurf.Vertices, sSurf.Faces, [0 0 0], P(i,:))';
        if isempty(proj)
            Q(i,:) = P(i,:);
        elseif (size(proj,1) > 1)
            Q(i,:) = proj(1,:);
        else
            Q(i,:) = proj;
        end
    end
    % Compute the distance between each point and its projection (distance to the scalp surface)
    if (nargout >= 2)
        dist = sum(sqrt((P - Q).^2), 2);
    end
    % Restore center
    Q = bst_bsxfun(@plus, Q, center);
end