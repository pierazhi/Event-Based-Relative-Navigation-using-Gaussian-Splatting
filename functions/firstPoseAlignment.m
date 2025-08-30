% ========================================================================
% Trajectory Evaluation and Pose Error Analysis Script
% Author: Piercarlo Fontana
%
% Description:
% This script loads ground truth and estimated trajectories (in KITTI pose
% format) and performs a complete evaluation of translational and rotational
% accuracy. It includes:
%   - Visualization of raw trajectories before and after alignment
%   - Umeyama alignment (with optional scale) for position correction
%   - First-pose orientation alignment for rotational consistency
%   - Per-axis translation error statistics (RMSE, mean, max)
%   - Per-axis orientation error evaluation in yaw, pitch, and roll
%   - Global angular drift computation across all frames
%   - Automated generation of LaTeX tables with translation and rotation
%     metrics for quick integration into reports or theses
%   - Optional saving of all generated figures to a specified directory
%
% The script provides both numerical metrics (printed to console and LaTeX)
% and plots for visual analysis, supporting fair comparison between estimated
% and reference trajectories in 3D space.
%

%
% Functions:
%   - load_kitti_poses: loads poses from a KITTI-style text file
%   - umeyama_alignment: performs least-squares similarity transformation
%   - saveAllFiguresToPath: exports all open figures as .eps/.png to a directory
%
% ========================================================================

clearvars; clc; close all;

% === GLOBAL FONT SETTINGS ===
AXES_FONT_NAME   = 'Helvetica';
AXES_FONT_SIZE   = 15; % ticks
LABEL_FONT_SIZE  = 15; % xlabel/ylabel/zlabel
LEGEND_FONT_SIZE = 10; % legends

% Helper to apply fonts to current axes
applyFontSettings = @(ax) set(ax, 'FontName', AXES_FONT_NAME, 'FontSize', AXES_FONT_SIZE, 'LineWidth', 1);
applyLabelSettings = @(h) set(h, 'FontName', AXES_FONT_NAME, 'FontSize', LABEL_FONT_SIZE);
applyLegendSettings = @(h) set(h, 'FontName', AXES_FONT_NAME, 'FontSize', LEGEND_FONT_SIZE);

% === Settings ===
blender_fix = false;
purple = [173/255, 7/255, 249/255];

% === Load poses ===
gt = load_kitti_poses("gt_office3.txt");
est = load_kitti_poses("est_office3.txt");
N = size(est, 3);

% === Apply Blender convention fix ===
if blender_fix
    R_flipZ = diag([1 1 -1]);
    for i = 1:N
        gt(1:3,1:3,i)  = gt(1:3,1:3,i) * R_flipZ;
        est(1:3,1:3,i) = est(1:3,1:3,i) * R_flipZ;
    end
end

gt_xyz = squeeze(gt(1:3, 4, :));
est_xyz = squeeze(est(1:3, 4, :));

figure('Name', 'Trajectory Before Alignment');
plot3(gt_xyz(1,:), gt_xyz(2, :), gt_xyz(3,:), 'k-', 'LineWidth', 2); hold on;
plot3(est_xyz(1,:), est_xyz(2,:), est_xyz(3,:), 'Color', purple, 'LineStyle', '--', 'LineWidth', 2);
hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]');
hL = legend('Ground Truth', 'Estimated', 'Location', 'best');
grid on; axis equal;
applyFontSettings(gca); applyLabelSettings([hX hY hZ]);  applyLegendSettings(hL);

% === Align only positions ===

[R_u, t_u, s_u] = umeyama_alignment(est_xyz, gt_xyz, true);
est_traj_aligned = zeros(4, 4, N);
for i = 1:N
    T = est(:,:,i);
    est_traj_aligned(1:3,1:3,i) = s_u * R_u * T(1:3,1:3);
    est_traj_aligned(1:3,4,i)   = s_u * R_u * T(1:3,4) + t_u;
    est_traj_aligned(4,4,i)     = 1;
end

% === First-pose orientation alignment ===
T_align = gt(:,:,1) / est(:,:,1);
est_orient_aligned = zeros(4,4,N);
for i = 1:N
    est_orient_aligned(:,:,i) = T_align * est(:,:,i);
end

% === Scene scale & quiver size ===
scene_min = min(gt_xyz, [], 2);
scene_max = max(gt_xyz, [], 2);
scene_diag = norm(scene_max - scene_min);
quiver_scale = 0.05 * scene_diag;

% === Translation error ===
gt_xyz_vec  = gt_xyz';
est_xyz_vec = squeeze(est_traj_aligned(1:3, 4, :))';
trans_error = sqrt(sum((gt_xyz_vec - est_xyz_vec).^2, 2));

% === Plot trajectories ===

figure('Name', 'Trajectory After Alignment');
plot3(gt_xyz_vec(:,1), gt_xyz_vec(:,2), gt_xyz_vec(:,3), 'k-', 'LineWidth', 2); hold on;
plot3(est_xyz_vec(:,1), est_xyz_vec(:,2), est_xyz_vec(:,3), 'Color', purple, 'LineStyle', '--', 'LineWidth', 2);
hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]');
hL = legend('Ground Truth', 'Estimated', 'Location', 'best');
grid on; axis equal;
applyFontSettings(gca); applyLabelSettings([hX hY hZ]);  applyLegendSettings(hL);

% === Translation error plot ===
% figure('Name', 'Error');
% plot(trans_error * 100, 'LineWidth', 2, 'color', purple);
% hX = xlabel('Frame Index'); hY = ylabel('Translation Error [cm]');
% grid on;
% applyFontSettings(gca); applyLabelSettings([hX hY]); 

% === Stats ===
rmse = sqrt(mean(trans_error.^2));
fprintf('Translation RMSE: %.6f cm\n', rmse*100);
fprintf('Max Translation Error: %.6f cm\n', max(trans_error)*100);
fprintf('Mean Translation Error: %.6f cm\n', mean(trans_error)*100);

diff_xyz_cm = (gt_xyz_vec - est_xyz_vec) * 100;
rmse_xyz = sqrt(mean(diff_xyz_cm.^2, 1));
max_xyz  = max(abs(diff_xyz_cm), [], 1);
labels = {'X', 'Y', 'Z'};
for i = 1:3
    fprintf('\n[%s axis]\n', labels{i});
    fprintf('  RMSE:   %.3f cm\n', rmse_xyz(i));
    fprintf('  Max:    %.3f cm\n', max_xyz(i));
end

% === Per-axis evolution plots ===
figure('Name', 'Per_Axis');
for j = 1:3
    subplot(3,1,j);
    plot(gt_xyz_vec(:,j), 'k-', 'LineWidth', 2); hold on;
    plot(est_xyz_vec(:,j), 'Color', purple, 'LineWidth', 2, 'LineStyle', '--');
    hX = xlabel('Frame Index'); hY = ylabel([labels{j} ' [m]']);
    hL = legend('Ground Truth', 'Estimated', 'Location', 'best');
    grid on;
    applyFontSettings(gca); applyLabelSettings([hX hY]);  applyLegendSettings(hL);
end

% === Orientation axes plot ===
% figure('Name', 'Orientation');
% hold on; axis equal; grid on;
% hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]');
% plot3(gt_xyz_vec(:,1), gt_xyz_vec(:,2), gt_xyz_vec(:,3), 'k-', 'LineWidth', 2);
% plot3(est_xyz_vec(:,1), est_xyz_vec(:,2), est_xyz_vec(:,3), 'Color', purple, 'LineStyle', '--', 'LineWidth', 2);
% axis_colors = {'r', 'g', 'b'};
% step = max(1, round(N / 10));
% for i = 1:step:N
%     origin_gt = gt(1:3,4,i);
%     origin_est = est_traj_aligned(1:3,4,i);
%     for a = 1:3
%         dir_gt = quiver_scale * (gt(1:3,a,i) / norm(gt(1:3,a,i)));
%         dir_est = quiver_scale * (est_orient_aligned(1:3,a,i) / norm(est_orient_aligned(1:3,a,i)));
%         quiver3(origin_gt(1), origin_gt(2), origin_gt(3), dir_gt(1), dir_gt(2), dir_gt(3), ...
%             'Color', axis_colors{a}, 'LineWidth', 1.2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
%         quiver3(origin_est(1), origin_est(2), origin_est(3), dir_est(1), dir_est(2), dir_est(3), ...
%             'Color', axis_colors{a}, 'LineStyle', '--', 'LineWidth', 1.2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
%     end
% end
% hL = legend('GT', 'Estimated', 'Location', 'best');
% applyFontSettings(gca); applyLabelSettings([hX hY hZ]);  applyLegendSettings(hL);

% === Orientation errors ===
rot_err_rpy = zeros(N,3);
rpy_gt = zeros(N,3); rpy_est = zeros(N,3);
for i = 1:N
    R_gt = gt(1:3,1:3,i);
    R_est = est_orient_aligned(1:3,1:3,i);
    rpy_gt(i,:)  = rotm2eul(R_gt, 'ZYX');
    rpy_est(i,:) = rotm2eul(R_est, 'ZYX');
end
rpy_gt = rad2deg(unwrap(rpy_gt));
rpy_est = rad2deg(unwrap(rpy_est));
rot_err_rpy = wrapTo180(rpy_est - rpy_gt);

labels = {'Yaw', 'Pitch', 'Roll'};
figure('Name', 'RPY');
for j = 1:3
    subplot(3,1,j);
    plot(rpy_gt(:,j), 'k-', 'LineWidth', 2); hold on;
    plot(rpy_est(:,j), 'Color', purple, 'LineStyle', '--', 'LineWidth', 2);
    hX = xlabel('Frame Index'); hY = ylabel([labels{j} ' [°]']);
    if j == 1
        hL = legend('Ground Truth', 'Estimated', 'Location', 'Best');
    end
    grid on;
    applyFontSettings(gca); applyLabelSettings([hX hY]);  applyLegendSettings(hL);
end

% === Orientation error stats ===
rpy_rmse = sqrt(mean(rot_err_rpy.^2, 1));
rpy_mean = mean(abs(rot_err_rpy), 1); 
rpy_max  = max(abs(rot_err_rpy), [], 1);
fprintf('\n📐 Rotation Errors per Axis [°]:\n');
for i = 1:3
    fprintf('[%s]\n', labels{i});
    fprintf('  RMSE: %.3f°\n', rpy_rmse(i));
    fprintf('  Mean: %.3f°\n', rpy_mean(i));
    fprintf('  Max:  %.3f°\n\n', rpy_max(i));
end

% === Angular difference plot ===
ang_diff = zeros(N, 1);
for i = 1:N
    R_gt = gt(1:3,1:3,i);
    R_est = est_orient_aligned(1:3,1:3,i);
    R_err = R_gt' * R_est;
    theta_rad = acos( max(-1, min(1, (trace(R_err) - 1) / 2)) );
    ang_diff(i) = rad2deg(theta_rad);
end

figure('Name', 'Drift');
plot(ang_diff, 'LineWidth', 2, 'Color', purple);
hX = xlabel('Frame Index'); hY = ylabel('Angular Error [°]');
grid on;
applyFontSettings(gca); applyLabelSettings([hX hY]); 

fprintf('🔁 Angular Difference (all axes combined):\n');
fprintf('  RMSE: %.3f°\n', sqrt(mean(ang_diff.^2)));
fprintf('  Drift @ End Frame %.3f°\n', max(ang_diff));
fprintf('  Mean: %.3f°\n', mean(ang_diff));

%%

latex = sprintf([
    '\\begin{minipage}{.5\\textwidth}\n' ...
    '   \\begin{table}[H]\n' ...
    '       \\centering\n' ...
    '           \\begin{tabular}{cc}\n' ...
    '               \\toprule\n' ...
    '               \\textbf{Translation Metric} & \\textbf{RMSE} \\\\\n' ...
    '               \\midrule\n' ...
    '                       \\textbf{X} & %.3f~cm \\\\\n' ...
    '                       \\textbf{Y} & %.3f~cm \\\\\n' ...
    '                       \\textbf{Z} & %.3f~cm \\\\\n' ...
    '               \\midrule\n' ...
    '                       \\textbf{Total} & %.3f~cm \\\\\n' ...
    '               \\bottomrule\n' ...
    '           \\end{tabular}\n' ...
    '           \\caption{Translation Metrics}\n' ...
    '           \\label{tab:pose_eval_translation}\n' ...
    '   \\end{table}\n' ...
    '\\end{minipage}\n' ...
    '\\hfill\n' ...
    '\\begin{minipage}{.5\\textwidth}\n' ...
    '   \\begin{table}[H]\n' ...
    '       \\centering\n' ...
    '           \\begin{tabular}{cc}\n' ...
    '               \\toprule\n' ...
    '               \\textbf{Rotation Metric} & \\textbf{RMSE} \\\\\n' ...
    '               \\midrule\n' ...
    '                   \\textbf{Yaw} & %.3f$^{\\circ}$  \\\\\n' ...
    '                   \\textbf{Pitch} & %.3f$^{\\circ}$  \\\\\n' ...
    '                   \\textbf{Roll} & %.3f$^{\\circ}$  \\\\\n' ...
    '               \\midrule\n' ...
    '                   \\textbf{Drift} & %.3f$^{\\circ}$  \\\\\n' ...
    '               \\bottomrule\n' ...
    '           \\end{tabular}\n' ...
    '           \\caption{Rotation Metrics}\n' ...
    '           \\label{tab:pose_eval_rotation}\n' ...
    '   \\end{table}\n' ...
    '\\end{minipage}\n' ...
    ], ...
    rmse_xyz(1), rmse_xyz(2), rmse_xyz(3), rmse*100, ...
    rpy_rmse(1), rpy_rmse(2), rpy_rmse(3), sqrt(mean(ang_diff.^2)) );

% Print once
fprintf('\n%s\n', latex);


% Optional: save to a .tex file with both tables
fid = fopen('pose_metrics_tables.tex','w');
if fid ~= -1
    fprintf(fid, '%s\n%s', latex);
    fclose(fid);
end

%%

saveAllFiguresToPath('/Users/pierazhi/Documents/ImagesTesi/')

%% Function

function poses = load_kitti_poses(filename)
    data = readmatrix(filename);
    N = size(data, 1);
    poses = zeros(4, 4, N);
    for i = 1:N
        T = reshape(data(i, :), [4, 3])';
        T = [T; 0 0 0 1];
        poses(:, :, i) = T;
    end
end

function [R, t, c] = umeyama_alignment(x, y, with_scale)
%UMEYAMA_ALIGNMENT Least-squares similarity transformation between two point sets
%   [R, t, c] = umeyama_alignment(x, y, with_scale)
%   x and y: m x n matrices (m = dimensions, n = points)
%   with_scale: true/false, whether to estimate scale
%
%   Returns:
%     R: rotation matrix (m x m)
%     t: translation vector (m x 1)
%     c: scale factor (scalar)

    if nargin < 3
        with_scale = false;
    end
    
    if ~isequal(size(x), size(y))
        error('Input matrices x and y must be of the same size');
    end
    
    [m, n] = size(x); % m = dimensions, n = number of points
    
    % Compute means
    mean_x = mean(x, 2); % column vector
    mean_y = mean(y, 2);
    
    % Variance of x
    sigma_x = (1/n) * sum(vecnorm(x - mean_x, 2, 1).^2);
    
    % Covariance matrix
    cov_xy = zeros(m, m);
    for i = 1:n
        cov_xy = cov_xy + (y(:,i) - mean_y) * (x(:,i) - mean_x)';
    end
    cov_xy = cov_xy / n;
    
    % SVD
    [U, D, V] = svd(cov_xy);
    
    % Handle degenerate case
    % if rank(D) < m - 1
    %     error('Degenerate covariance rank, Umeyama alignment is not possible');
    % end
    
    % Construct S matrix
    S = eye(m);
    if det(U) * det(V) < 0
        S(end, end) = -1;
    end
    
    % Rotation
    R = U * S * V';
    
    % Scale
    if with_scale
        c = trace(D * S) / sigma_x;
    else
        c = 1.0;
    end
    
    % Translation
    t = mean_y - c * R * mean_x;

end

function saveAllFiguresToPath(userDefinedPath)
% SAVEALLFIGURESTOPATH Saves all open figures to a user-defined path using the figure name as filename.
% Maximizes each figure before saving.
%
%   saveAllFiguresToPath(userDefinedPath) saves all open MATLAB figures to the
%   specified userDefinedPath. Filenames are based on the figure's Name property.
%   Spaces are replaced with underscores. Falls back to 'Figure_N' if no name is found.

% --- Input validation
if nargin == 0
    error('User-defined path is required.');
end

if ~(ischar(userDefinedPath) || (isstring(userDefinedPath) && isscalar(userDefinedPath)))
    error('User-defined path must be a non-empty string or char array.');
end

userDefinedPath = char(userDefinedPath); % Ensure char format

if ~exist(userDefinedPath, 'dir')
    mkdir(userDefinedPath);
    disp(['Created directory: ' userDefinedPath]);
end

% --- Get all open figures
openFigures = findall(0, 'Type', 'figure');

if isempty(openFigures)
    disp('No open figures found.');
    return
end

% --- Save each figure
for i = 1:numel(openFigures)
    currentFigure = openFigures(i);

    % Maximize figure window
    %set(currentFigure, 'WindowState', 'maximized');
    drawnow; % Ensure resizing is applied

    % Get figure name
    figName = get(currentFigure, 'Name');
    if ischar(figName) || isstring(figName)
        figName = strtrim(char(figName));
    else
        figName = '';
    end

    % Default fallback if no name is set
    if isempty(figName)
        fileName = sprintf('Figure_%d', i);
    else
        % Replace spaces with underscores and strip illegal characters
        fileName = regexprep(figName, '\s+', '_');
        fileName = regexprep(fileName, '[^\w]', '_'); % Keep only letters, numbers, and underscores
    end

    % Construct full path
    fullFilePath = fullfile(userDefinedPath, [fileName '.eps']);

    % Save as high-quality EPS
    try
        exportgraphics(currentFigure, fullFilePath, 'Resolution', 300);
        disp(['Figure ' num2str(i) ' saved to: ' fullFilePath]);
    catch ME
        disp(['Error saving figure ' num2str(i) ': ' ME.message]);
    end
end
end
