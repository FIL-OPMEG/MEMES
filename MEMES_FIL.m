function MEMES_FIL(dir_name,headshape_downsampled,...
    path_to_MRI_library,method,varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRI Estimation for MEG Sourcespace (MEMES)
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
% - dir_name              = directory for saving
% - headshape_downsampled = headshape downsampled to 100-200 scalp points
% - path_to_MRI_library   = path to HCP MRI library
% - method                = method for creating pseudo head- and
%                           source-model: 'best' or 'average'
%
%%%%%%%%%%%%%%%%%%
% Variable Inputs:
%%%%%%%%%%%%%%%%%%
%
% - weight_face           = how much do you want to weight towards the
%                           facial information (1 = no weighting;
%                           10 = very high weighting. RS recommends
%                           weight_face = 3;
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
% - shape                   = headshape and fiducial information
% - headshape_downsampled   = headshape downsampled to 100 points
% - trans_matrix            = transformation matrix applied to headmodel
%                           and sourcemodel
% - sourcemodel3d           = sourcemodel warped to MNI space
% - headmodel               = singleshell headmodel (10000 vertices)
%
%%%%%%%%%%%%%%%%%%%%%
% Other Information:
%%%%%%%%%%%%%%%%%%%%%
%
% Example function call:
% MEMES_FIL(dir_name,grad_trans,headshape_downsampled,...
%    path_to_MRI_library,method,[0.98:1.02],8)

% I have introduced a variable scaling parameter for the MRIs to
% help with coregistration. For example to apply -2% to +2% scaling to
% every MRI specify: scaling = [0.98:0.01:1.2].
%
% However NOTE: the more scaling factors you apply the longer it will take
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['\nThis is MEMES for data collected from FIL OPM Lab...\n\nMake sure you have asked Robert'...
    'for an MRI library\n\n']);
warning('on')

%% Check inputs
disp('Performing input check');
% If Path to MRI library doesn't end with / or \ throw up and error
if ismember(path_to_MRI_library(end),['/','\']) == 0
    error('!!! Path to MRI library must end with / or \ !!!');
end


% If variable inputs are empty use defaults
if isempty(varargin)
    weight_face         = [];
else
    weight_face         = varargin{1};
end

% Convert headshape_downsampled to mm if required
if headshape_downsampled.unit ~= 'mm'
    headshape_downsampled = ft_convert_units(headshape_downsampled,'mm');
end

%% Extract subject names from your MRI library

try
    cd(path_to_MRI_library);
    % Get a list of all files and folders in this folder.
    files = dir(path_to_MRI_library);
    files(1:2) = [];
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    
    % Now these names to a variable called subject
    subject = [];
    
    for sub = 1 : length(subFolders)
        subject{sub} = subFolders(sub).name;
    end
    
    fprintf('%d subjects found in the MRI library: from %s to %s\n',...
        length(subject),subject{1}, subject{end});
    
catch
    warning('Something is wrong with your MRI library... Check the path!\n');
end

% Now try to load relevent information from the first subject
fprintf('Now checking the MRI library is organised correctly...\n');
try
    load([path_to_MRI_library subject{1} '/mesh.mat']);
    load([path_to_MRI_library subject{1} '/mri_realigned.mat']);
    clear mesh mri_realigned
    fprintf('...Subject %s is organised correctly!\n',subject{1});
    
catch
    warning('Your MRI library is not organised correctly');
    disp('Each folder should contain: mesh.mat, mri_realigned.mat');
end

%% CD to the right place

% CD to right place
cd(dir_name); fprintf('\nCDd to the right place\n');

%% Set up varaibles outside the pseudo MRI loop

% Error term variable - MEMES will crash here if your MRI library path is
% wrong..
error_term = zeros(1,length(subject));
% Variable to hold the transformation matrices
trans_matrix_library = [];

%% Facial Weighting?
% Weight towards facial information, if specified
if ~isempty(weight_face)
    % Find facial points
    count_facialpoints = find(headshape_downsampled.pos(:,3)<30 &...
        headshape_downsampled.pos(:,1)>70);
    
    % Create an array
    w = ones(size(headshape_downsampled.pos,1),1).* (1/weight_face);
    % Replace facial points with 1
    w(count_facialpoints) = 1;
    weights = @(x)assignweights(x,w);
    fprintf('Applying Weighting of %.2f \n',weight_face);
end

%% Subject Loop
% For each subject...
for m = 1:length(subject)
    
   disp(m);
    
    % Load the mesh
    sss = load(fullfile(path_to_MRI_library, subject{m},'mesh.mat'));
    mesh = sss.mesh;
    load([path_to_MRI_library subject{m} '/fids_SPM_convert.mat']);
    
    % Perform initial realign based on FIDS
    trans_fids = warp_fid(headshape_downsampled,fids_SPM_convert);
    
    mesh2 = ft_transform_geometry(trans_fids,mesh);
    
    %     figure; ft_plot_mesh(mesh2,'facealpha',0.4);
    %     ft_plot_headshape(headshape_downsampled); view([90 0]); camlight;
    %     ft_plot_mesh(fids_SPM_convert,'vertexcolor','g','vertexsize',40);
    
    numiter = 30;
    
    % Perform ICP
    % If we are applying weighting...
    if ~isempty(weight_face)
        [R, t, err, ~, ~] = icp(mesh2.pos', ...
            headshape_downsampled.pos', numiter, 'Minimize', 'plane',...
            'Extrapolation', true,'Weight', weights,'WorstRejection', 0.05);
        % If not applying weighting...
    else
        [R, t, err, ~, ~] = icp(mesh2.pos', ...
            headshape_downsampled.pos', numiter, 'Minimize', 'plane',...
            'Extrapolation', true,'WorstRejection', 0.1);
    end
    
    % Update error term
    error_term(m) = err(end);
    
    % Combine trans matrices
    trans_matrix_temp = inv([real(R) real(t);0 0 0 1]);
    trans_matrix_library{m}  = trans_fids*trans_matrix_temp;
    
    % Clear mesh for next loop
    clear mesh mesh2 fids_SPM_convert
end

fprintf(' Finished the iterations\n');

%% Make pretty figure
fprintf('\n Finding good, OK and bad examples\n');

error_term_sorted = sort(error_term, 'ascend');
middle_num = length(error_term_sorted)/2;
winners = find(ismember(error_term,error_term_sorted(1:3)));
middles = find(ismember(error_term,error_term_sorted(middle_num-1:middle_num+1)));
losers = find(ismember(error_term,error_term_sorted(end-2:end)));

concat = [winners middles losers];

% Create figure to summarise the losers,middles and winners
figure;
for i = 1:9
    load([path_to_MRI_library subject{(concat(i))} '/mesh.mat'])
    mesh_spare = mesh;
    mesh_spare.pos = ft_warp_apply(trans_matrix_library{(concat(i))}, mesh_spare.pos);
    
    subplot(3,3,i)
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; view([-270,-10]);
    if ismember(i,1:3)
        title(sprintf('BEST: %.3f', error_term((concat(i)))));
    elseif ismember(i,4:6)
        title(sprintf('MIDDLE: %.3f', error_term((concat(i)))));
    elseif ismember(i,7:9)
        title(sprintf('WORST: %.3f', error_term((concat(i)))));
    end
    
    ft_plot_headshape(headshape_downsampled);
    
    clear mesh mesh_spare
    
    if i == 9
        print('best_middle_worst_examples','-dpng','-r100');
    end
end


switch method
    case 'average'
        fprintf('USE WITH CAUTION - Still testing \n');
        
        % Average over how many? N=20 the best?
        average_over_n = 20;
        % Variable to hold average sourcemodel .pos
        average_sourcemodel_all = [];
        % Variable to hold average sourcemodel .pos
        average_headmodel_all = [];
        
        for rep = 1:average_over_n
            % Find the number of the nth MRI
            winner_rep = find(ismember(error_term,error_term_sorted(rep)));
            
            % Update the user
            fprintf('Loaded MRI %d of %d : %s ... Scaling factor: %.2f\n',...
                rep,average_over_n,subject{winner_rep},...
                scaling_factor_all(winner_rep));
            
            % Get the transformation matrix of the winner
            trans_matrix = trans_matrix_library{winner_rep};
            
            %% Get mesh
            % Get facial mesh of 1st winner
            if rep == 1
                load([path_to_MRI_library subject{winner_rep} '/mesh.mat'])
                mesh.pos = ft_warp_apply([scaling_factor_all(winner_rep) 0 0 0;0 ...
                    scaling_factor_all(winner_rep) 0 0; 0 0 scaling_factor_all(winner_rep) 0;...
                    0 0 0 1],mesh.pos);
                mesh.pos = ft_warp_apply(trans_matrix, mesh.pos);
                mesh_spare = mesh;
            end
            
            clear mesh
            
            %% Create Headmodel (in mm)
            load([path_to_MRI_library subject{winner_rep} '/headmodel.mat']);
            
            % Scale
            headmodel.bnd.pos = ft_warp_apply([scaling_factor_all(winner_rep) 0 0 0;0 ...
                scaling_factor_all(winner_rep) 0 0; 0 0 scaling_factor_all(winner_rep) 0; 0 0 0 1],...
                headmodel.bnd.pos);
            
            % Transform (MESH --> coreg via ICP adjustment)
            headmodel.bnd.pos = ft_warp_apply(trans_matrix,headmodel.bnd.pos);
            
            % Add the pos field to the array outside the loop
            average_headmodel_all(rep,:,:) = headmodel.bnd.pos(:,:);
            
            % Reserve the first headmodel for later
            if rep == 1
                headmodel_for_outside_loop = headmodel;
            end
            
            clear headmodel
            
            %% Create Sourcemodel (in mm)
            
            % Load specified sized sourcemodel
            load([path_to_MRI_library ...
                subject{winner_rep} '/sourcemodel3d_' num2str(sourcemodel_size) 'mm.mat']);
            
            % Scale
            sourcemodel3d.pos = ft_warp_apply([scaling_factor_all(winner_rep)...
                0 0 0;0 scaling_factor_all(winner_rep) 0 0; 0 0 ...
                scaling_factor_all(winner_rep) 0; 0 0 0 1],sourcemodel3d.pos);
            
            % Transform (MESH --> coreg via ICP adjustment)
            sourcemodel3d.pos = ft_warp_apply(trans_matrix,sourcemodel3d.pos);
            
            average_sourcemodel_all(rep,:,:) = sourcemodel3d.pos;
            
            % Reserve the first headmodel for later
            if rep == 1
                sourcemodel_for_outside_loop = sourcemodel3d;
            end
            
            clear trans_matrix sourcemodel3d winner_rep
            
        end
        
        % Average Headmodel
        fprintf('Averaging Headmodel\n');
        headmodel = headmodel_for_outside_loop;
        headmodel.bnd.pos = squeeze(mean(average_headmodel_all,1));
        
        % Average Sourcemodel
        fprintf('Averaging Sourcemodel\n');
        sourcemodel3d = sourcemodel_for_outside_loop;
        sourcemodel3d.pos = squeeze(mean(average_sourcemodel_all,1));
        
        % Create figure to check headodel and sourcemodel match
        figure;
        ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
        alpha 0.4; camlight;
        ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:),'vertexsize',5);
        view([0 0]);
        
        view_angle = [0 90 180 270];
        
        % Create figure to show final coregiration (with mesh of 1st place
        % MRI)
        figure; hold on;
        for rep = 1:4
            subplot(2,2,rep);
            ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.6; camlight;
            ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:),'vertexsize',3);
            ft_plot_headshape(headshape_downsampled) %plot headshape
            view([view_angle(rep),0]);
            ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.5);
            camlight; lighting phong; material dull;
        end
        
        print('coregistration_volumetric_quality_check','-dpng','-r100');
        
        %% SAVE
        fprintf('\nSaving the necessary data\n');
        
        save headmodel headmodel
        %save trans_matrix trans_matrix
        save sourcemodel3d sourcemodel3d
        %save mri_realigned_MEMES mri_realigned_MEMES
        
        fprintf('\nCOMPLETED - check the output for quality control\n');
        
        
    case 'best'
        
        % Find the MRI with the lowest ICP error between Polhemus points
        % and 3D scalp mesh
        winner = find(error_term == min(min(error_term)));
        fprintf('\nThe winning MRI is number %d of %d : %s\n',winner,length(subject),subject{winner});
        
        % Get the transformation matrix of the winner
        trans_matrix = trans_matrix_library{winner};
        
        % Get facial mesh of winner
        load([path_to_MRI_library subject{winner} '/mesh.mat'])
        mesh_spare = mesh;
        mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);
        
        % Get MRI of winning subject
        fprintf('Transforming the MRI\n');
        load([path_to_MRI_library subject{winner} '/mri_realigned.mat'],'mri_realigned');
        mri_realigned_MEMES = ft_transform_geometry(trans_matrix,...
            mri_realigned);
        
        % Create figure to show final coregiration
        figure; hold on;
        view_angle = [0 180 90 -90];
        for rep = 1:4
            subplot(2,2,rep);
            ft_plot_headshape(headshape_downsampled) %plot headshape
            view([view_angle(rep),0]);
            ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.5);
            camlight; lighting phong; material dull;
        end
        
        print('coregistration_volumetric_quality_check','-dpng','-r100');

        %% SAVE
        fprintf('\nSaving the necessary data\n');
        
        save mri_realigned_MEMES mri_realigned_MEMES
        
        % Export MRI
        fprintf('\nExporting MRI...\n');
        cfg             = [];
        cfg.parameter   = 'anatomy';
        cfg.filename    = 'mri_realigned_MEMES';
        cfg.filetype    = 'nifti';
        cfg.spmversion  = 'spm12';
        cfg.datatype    = 'double';
        ft_volumewrite(cfg,mri_realigned_MEMES);
        
        fprintf('\nCOMPLETED - check the output for quality control\n');
        
    otherwise
        fprintf('Something went wrong - did you specify *average* or *best*')
end
end