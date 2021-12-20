function [trans_fids] = warp_fid(headshape,fids_SPM_convert)
    % Function to warp the SPM canonical mesh based on fiducials supplied by
    % the headshape data
    elec2common  = ft_headcoordinates(headshape.fid.pos(1,:),...
        headshape.fid.pos(2,:), ...
        headshape.fid.pos(3,:));

    templ2common = ft_headcoordinates(fids_SPM_convert(1,:),...
        fids_SPM_convert(2,:), ...
        fids_SPM_convert(3,:));

    % compute the combined transform
    trans_fids       = inv(templ2common \ elec2common);

    end