% ========================================================================
% Ground-Truth Analytics and Estimated Trajectory Evaluation
% Author: Piercarlo Fontana
%
% Description:
% This script loads a ground-truth sequence from a MAT file (poses, timestamps,
% velocities, positions) and computes the nominal sampling period and mean frame
% rate, then evaluates an estimated trajectory against ground truth using a
% similarity transformation for positional alignment and a first-pose rotation
% alignment for attitude comparison. The pipeline reads KITTI-format text files
% for both GT and Estimated poses, applies Umeyama alignment on positions to
% remove global offset, rotation, and scale, and then computes translation RMSE,
% per-axis RMSE and maxima, and produces diagnostic 3D trajectory plots and
% per-axis evolution figures. For orientation, it aligns the estimated pose
% sequence to the ground-truth first pose, unwraps ZYX Euler angles over time,
% reports yaw/pitch/roll errors (RMSE, mean, max), and computes the geodesic
% angular difference per frame to summarize overall rotational discrepancy.
% The script also derives velocity from the aligned estimated trajectory using
% timestamp differences, compares it to ground-truth velocities on the common
% window, builds a per-frame velocity error norm, normalizes it by altitude,
% and reports the resulting scalar score E(i) as a compact proxy for dynamic
% consistency. Finally, it measures total path length for both GT and Estimated
% trajectories and expresses translation RMSE as a percentage of GT path length
% for scale-aware reporting. Publication-ready LaTeX tables are generated for
% translation and rotation metrics, and figures can be batch-exported to disk.
%
% Inputs:
% - 0000_groundtruth.mat must define: poses (Nx4x4), timestamps (Nx1),
%   velocities ((N-1)x3), positions (Nx3). Timestamps are strictly
%   increasing; their mean spacing defines the nominal frequency.
% - GT and Estimated text files are in KITTI pose format (Nx12 row-wise
%   SE(3) flattened; loader reshapes to 4x4 per frame). start_idx must be set
%   to align the time window between GT and Estimated; end_idx is derived.
%
% Outputs:
% - Console metrics for translation (total and per-axis), orientation, 
%   velocity score E(i), path lengths, and RMSE normalized by path length.
% - Figures named consistently for trajectories, per-axis evolution, RPY, and
%   velocity comparison; saveAllFiguresToPath can export PNGs to a target folder.
% - Two LaTeX tables (booktabs) for direct inclusion in papers or theses.
%
% Implementation notes:
% - Position-only Umeyama alignment estimates rotation, translation, and
%   optional scale; orientation comparison uses first-pose alignment to
%   isolate relative drift. Angle unwrapping is applied before RPY error
%   computation and errors are wrapped back to [-180, 180] degrees.
% - Velocity is computed via finite differences of aligned positions using
%   timestamp deltas; indexing is matched to ground-truth arrays to avoid
%   off-by-one issues when forming norms and normalization by altitude.
% ========================================================================


clc; close all; clearvars;

% === Load Ground Truth ===
load('0000_groundtruth.mat');  % Contains: poses (Nx4x4), timestamps (Nx1), velocities ((N-1)x3), positions (Nx3)

timestamps = timestamps(:);    % Ensure column vector
N = size(poses, 1);
purple = [173/255, 7/255, 249/255];

% === GLOBAL FONT SETTINGS ===
AXES_FONT_NAME   = 'Helvetica';
AXES_FONT_SIZE   = 12; % ticks
LABEL_FONT_SIZE  = 14; % axis labels
LEGEND_FONT_SIZE = 10; % legends

% Helpers to apply fonts
applyFontSettings   = @(ax) set(ax, 'FontName', AXES_FONT_NAME, 'FontSize', AXES_FONT_SIZE, 'LineWidth', 1);
applyLabelSettings  = @(h)  set(h,  'FontName', AXES_FONT_NAME, 'FontSize', LABEL_FONT_SIZE);
applyLegendSettings = @(h)  set(h,  'FontName', AXES_FONT_NAME, 'FontSize', LEGEND_FONT_SIZE);

% === Compute overall frame frequency ===
dt_all = diff(timestamps);      % all time intervals
mean_dt = mean(dt_all);         % average time between timestamps
freq_hz = 1 / mean_dt;          % frequency in Hz

fprintf('🕒 Mean time step: %.6f s (%.2f Hz)\n', mean_dt, freq_hz);

% === Define Parameters ===
T_total = timestamps(end);  % Optional, used only if timestamps aren't regular

% 
% % === Ground Truth Analysis ===
% % 1. Plot 3D trajectory
% figure;
% p1 = plot3(positions(:,1), positions(:,2), positions(:,3), 'k-', 'LineWidth', 1.5); hold on;
% p2 = plot3(positions(start_idx:end_idx,1), positions(start_idx:end_idx,2), positions(start_idx:end_idx,3), 'r-', 'LineWidth', 2.5);
% hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]');
% hL = legend('Full GT', 'Compared Segment', 'Location', 'best'); grid on;
% applyFontSettings(gca); applyLabelSettings([hX hY hZ]);  applyLegendSettings(hL);
% 
% % 2. Plot velocity
% figure;
% subplot(3,1,1);
% plot(1:N-1, velocities(:,1), 'k'); hold on;
% plot(start_idx:end_idx-1, velocities(start_idx:end_idx-1,1), 'r', 'LineWidth', 2);
% hY = ylabel('V_x [m/s]'); grid on; applyFontSettings(gca); applyLabelSettings(hY);
% 
% subplot(3,1,2);
% plot(1:N-1, velocities(:,2), 'k'); hold on;
% plot(start_idx:end_idx-1, velocities(start_idx:end_idx-1,2), 'r', 'LineWidth', 2);
% hY = ylabel('V_y [m/s]'); grid on; applyFontSettings(gca); applyLabelSettings(hY);
% 
% subplot(3,1,3);
% plot(1:N-1, velocities(:,3), 'k'); hold on;
% plot(start_idx:end_idx-1, velocities(start_idx:end_idx-1,3), 'r', 'LineWidth', 2);
% hX = xlabel('Frame Index'); hY = ylabel('V_z [m/s]');
% applyFontSettings(gca); applyLabelSettings([hX hY]); 
% 
% % 3. Altitude
% altitude = abs(positions(:,3));
% figure;
% plot(1:N, altitude, 'k-'); hold on;
% plot(start_idx:end_idx, altitude(start_idx:end_idx), 'r', 'LineWidth', 2);
% hX = xlabel('Frame Index'); hY = ylabel('|Z| [m]'); grid on;
% applyFontSettings(gca); applyLabelSettings([hX hY]); 
% 
% % 4. Angular velocity (kept commented in your code)
% % figure;
% % plot(1:N-1, angular_velocity, 'LineWidth', 1.5);
% % xlabel('Frame Index');
% % ylabel('Angular Velocity [rad/s]');
% % legend('\omega_x', '\omega_y', '\omega_z');
% % title('Angular Velocity Over Time'); grid on;
% 
% % 5. Euler angles
% euler_angles = zeros(N, 3); angular_velocity = zeros(N-1, 3);
% for i = 1:N
%     R = squeeze(poses(i,1:3,1:3));
%     euler_angles(i,:) = rotm2eul(R, 'XYZ');
% end
% for i = 1:N-1
%     R0 = squeeze(poses(i,1:3,1:3)); R1 = squeeze(poses(i+1,1:3,1:3));
%     R_rel = R0' * R1; axang = rotm2axang(R_rel);
%     angular_velocity(i,:) = axang(1:3) * axang(4) / (timestamps(i+1) - timestamps(i));
% end
% 
% figure;
% plot(1:N, rad2deg(euler_angles), 'k-'); hold on;
% plot(start_idx:end_idx, rad2deg(euler_angles(start_idx:end_idx,:)), 'LineWidth', 2);
% hX = xlabel('Frame Index'); hY = ylabel('Euler Angles [deg]');
% hL = legend('Full Roll', 'Full Pitch', 'Full Yaw', 'Seg Roll', 'Seg Pitch', 'Seg Yaw', 'Location', 'best');
% grid on;
% applyFontSettings(gca); applyLabelSettings([hX hY]); applyLegendSettings(hL); 

% === Estimated Trajectory Comparison ===
close all; 

% === Settings ===
blender_fix = false;


% === Load poses ===
gt = load_kitti_poses("gt_ELOPE_0000.txt");
est = load_kitti_poses("est_ELOPE_0000.txt");
N = size(est, 3);
start_idx = 0;             % <-- You must set this
end_idx = start_idx + N;

% === Apply Blender convention fix (flip camera Z) ===
if blender_fix
    R_flipZ = diag([1 1 -1]);
    for i = 1:N
        gt(1:3,1:3,i)  = gt(1:3,1:3,i) * R_flipZ;
        est(1:3,1:3,i) = est(1:3,1:3,i) * R_flipZ;
    end
end

% === Align only positions using Umeyama ===
gt_xyz = squeeze(gt(1:3, 4, :));
est_xyz = squeeze(est(1:3, 4, :));
[R_u, t_u, s_u] = umeyama_alignment(est_xyz, gt_xyz, true);
est_traj_aligned = zeros(4, 4, N);
for i = 1:N
    T = est(:,:,i);
    est_traj_aligned(1:3,1:3,i) = s_u * R_u * T(1:3,1:3);
    est_traj_aligned(1:3,4,i)   = s_u * R_u * T(1:3,4) + t_u;
    est_traj_aligned(4,4,i)     = 1;
end

% === First-pose orientation alignment ===
T_align = gt(:,:,1) * inv(est(:,:,1));  % only rotation correction
est_orient_aligned = zeros(4,4,N);
for i = 1:N
    est_orient_aligned(:,:,i) = T_align * est(:,:,i);
end

% === Compute scene scale and set quiver size ===
scene_min = min(gt_xyz, [], 2);
scene_max = max(gt_xyz, [], 2);
scene_diag = norm(scene_max - scene_min);
quiver_scale = 0.05 * scene_diag;

% === Compute translation error ===
gt_xyz_vec  = gt_xyz';  % Nx3
est_xyz_vec = squeeze(est_traj_aligned(1:3, 4, :))';
trans_error = sqrt(sum((gt_xyz_vec - est_xyz_vec).^2, 2));

% === Plot trajectories ===
figure('Name', 'Trajectory')
plot3(gt_xyz_vec(:,1), gt_xyz_vec(:,2), gt_xyz_vec(:,3), 'k-', 'LineWidth', 2); hold on;
plot3(est_xyz_vec(:,1), est_xyz_vec(:,2), est_xyz_vec(:,3), 'Color', purple, 'Linestyle', '--', 'LineWidth', 2);
hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]');
hL = legend('Ground Truth', 'Estimated', 'Location', 'best'); grid on;
applyFontSettings(gca); applyLabelSettings([hX hY hZ]);  applyLegendSettings(hL);

% === RMSE & per-axis stats ===
rmse = sqrt(mean(trans_error.^2));
fprintf('Translation RMSE: %.6f m\n', rmse);
fprintf('Max Translation Error: %.6f m\n', max(trans_error));
fprintf('Mean Translation Error: %.6f m\n', mean(trans_error));

diff_xyz = (gt_xyz_vec - est_xyz_vec);
rmse_xyz = sqrt(mean(diff_xyz.^2, 1));
max_xyz  = max(abs(diff_xyz), [], 1);
labels = {'X', 'Y', 'Z'};
for i = 1:3
    fprintf('\n[%s axis]\n', labels{i});
    fprintf('  RMSE:   %.3f m\n', rmse_xyz(i));
    fprintf('  Max:    %.3f m\n', max_xyz(i));
end

% === Per-axis evolution plots ===
figure('Name','Per_axis');
for j = 1:3
    subplot(3,1,j);
    plot(abs(gt_xyz_vec(:,j)), 'k-', 'LineWidth', 2); hold on;
    plot(abs(est_xyz_vec(:,j)), 'Color', purple, 'Linestyle', '--', 'LineWidth', 2);
    hX = xlabel('Frame Index'); hY = ylabel([labels{j} ' [m]']);
    hL = legend('Ground Truth', 'Estimated', 'Location', 'best'); legend boxoff; grid on;
    applyFontSettings(gca); applyLabelSettings([hX hY]);  applyLegendSettings(hL);
end

% === Orientation Axis Plot ===
% figure; hold on; axis equal; grid on;
% hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]');
% plot3(gt_xyz_vec(:,1), gt_xyz_vec(:,2), gt_xyz_vec(:,3), 'k-', 'LineWidth', 2);
% plot3(est_xyz_vec(:,1), est_xyz_vec(:,2), est_xyz_vec(:,3), 'Color', purple, 'Linestyle', '--', 'LineWidth', 2);
% axis_colors = {'r', 'g', 'b'};
% step = max(1, round(N / 10));
% for i = 1:step:N
%     origin_gt = gt(1:3,4,i);
%     origin_est = est_traj_aligned(1:3,4,i);
%     for a = 1:3
%         % Normalize and scale axis arrows
%         dir_gt = quiver_scale * (gt(1:3,a,i) / norm(gt(1:3,a,i)));
%         dir_est = quiver_scale * (est_orient_aligned(1:3,a,i) / norm(est_orient_aligned(1:3,a,i)));
%         quiver3(origin_gt(1), origin_gt(2), origin_gt(3), ...
%                 dir_gt(1), dir_gt(2), dir_gt(3), ...
%                 'Color', axis_colors{a}, 'LineWidth', 1.2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
%         quiver3(origin_est(1), origin_est(2), origin_est(3), ...
%                 dir_est(1), dir_est(2), dir_est(3), ...
%                 'Color', axis_colors{a}, 'LineStyle', '--', ...
%                 'LineWidth', 1.2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
%     end
% end
% hL = legend('GT', 'Estimated','Location', 'best');
% applyFontSettings(gca); applyLabelSettings([hX hY hZ]);  applyLegendSettings(hL);

% === Orientation errors ===
rot_err_rpy   = zeros(N,3);
rpy_gt = zeros(N,3);
rpy_est = zeros(N,3);
for i = 1:N
    R_gt = gt(1:3,1:3,i);
    R_est = est_orient_aligned(1:3,1:3,i);
    rpy_gt(i,:)  = rotm2eul(R_gt, 'ZYX');
    rpy_est(i,:) = rotm2eul(R_est, 'ZYX');
end

% === Unwrap angles along time axis ===
rpy_gt = rad2deg(unwrap(rpy_gt));   % [N x 3], unwrapped and converted to deg
rpy_est = rad2deg(unwrap(rpy_est)); % same

% === Compute RPY error (angle difference) ===
rot_err_rpy = wrapTo180(rpy_est - rpy_gt);  % wrap to [-180,180] deg

% === RPY Evolution Plot ===
labels = {'Yaw [°]', 'Pitch [°]', 'Roll [°]'};
figure('Name', 'RPY');
for j = 1:3
    subplot(3,1,j);
    plot(rpy_gt(:,j), 'k-', 'LineWidth', 2); hold on;
    plot(rpy_est(:,j), 'Color', purple, 'Linestyle', '--', 'LineWidth', 2);
    hY = ylabel([labels{j}]); hX = xlabel('Frame Index');
    hL = legend('Ground Truth', 'Estimated', 'Location', 'best');
    grid on; legend boxoff;
    applyFontSettings(gca); applyLabelSettings([hX hY]);  applyLegendSettings(hL);
end

labels = {'Yaw (Z)', 'Pitch (Y)', 'Roll (X)'};
rpy_rmse = sqrt(mean(rot_err_rpy.^2, 1));
rpy_mean = mean(abs(rot_err_rpy), 1);  % abs in case bias flips
rpy_max  = max(abs(rot_err_rpy), [], 1);

fprintf('\n📐 Rotation Errors per Axis (deg):\n');
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
    theta_rad = acos(max(-1, min(1, (trace(R_err) - 1) / 2))); % Ensure value is within valid range for acos
    ang_diff(i) = rad2deg(theta_rad);
end

% figure('Name', 'Drift');
% plot(ang_diff, 'LineWidth', 2, 'Color', purple);
% hX = xlabel('Frame Index'); hY = ylabel('Angular Error [°]');
% grid on;
% applyFontSettings(gca); applyLabelSettings([hX hY]); 

fprintf('🔁 Angular Difference (all axes combined):\n');
fprintf('  RMSE: %.3f°\n', sqrt(mean(ang_diff.^2)));
fprintf('  Mean: %.3f°\n', mean(ang_diff));
fprintf('  Max:  %.3f°\n\n', max(ang_diff));

dt = diff(timestamps(start_idx:end_idx-1));   % [N-1 x 1]
est_vel = diff(est_xyz_vec) ./ dt;            % dt is [N-1 x 1]
gt_vel_used = velocities(start_idx : end_idx - 2, :);
z_gt = abs(positions(start_idx:end_idx-2, 3));     % also fix this to match
diff_vel = est_vel - gt_vel_used;
rms_vel_err = sqrt(sum(diff_vel.^2, 2));           % norm of velocity error per frame
norm_vel_err = rms_vel_err ./ z_gt;                % normalized by altitude
E_i = mean(norm_vel_err);                          % final score
fprintf('📉 Velocity Score E(i): %.4f\n', E_i);

figure('Name','Velocity');
vlabels = {'V_x', 'V_y', 'V_z'};
for d = 1:3
    subplot(3,1,d);
    plot(gt_vel_used(:,d), 'k-', 'LineWidth', 2); hold on;
    plot(est_vel(:,d), 'Color', purple, 'Linestyle', '--', 'LineWidth', 2);
    hX = xlabel('Frame Index'); hY = ylabel([vlabels{d} ' [m/s]']);
    hL = legend('Ground Truth', 'Estimated','Location', 'best');
    grid on; legend boxoff;
    applyFontSettings(gca); applyLabelSettings([hX hY]); applyLegendSettings(hL); 
end

% === Compute total path length for GT and Estimated ===
gt_deltas  = diff(gt_xyz_vec);                        % [N-1 x 3]
est_deltas = diff(est_xyz_vec);                       % [N-1 x 3]

gt_path_length  = sum(sqrt(sum(gt_deltas.^2, 2)));    % scalar
est_path_length = sum(sqrt(sum(est_deltas.^2, 2)));   % scalar

fprintf('GT path length:  %.3f m\n', gt_path_length);
fprintf('Est path length: %.3f m\n', est_path_length);

% === RMSE normalized by GT path length (in %) ===
rmse_norm_pct = (rmse / gt_path_length) * 100;
fprintf('Normalized RMSE w.r.t. path length: %.2f %%\n', rmse_norm_pct);


%%

% === Translation RMSE table ===
latex_tbl_trans = sprintf([ ...
    '\\begin{table}[H]\n' ...
    '    \\centering\n' ...
    '        \\begin{tabular}{cc}\n' ... % <-- both columns centered
    '            \\toprule\n' ...
    '            \\textbf{Translation Metric} & \\textbf{RMSE [cm]} \\\\\n' ...
    '            \\midrule\n' ...
    '            \\textbf{Total} & %.3f~cm \\\\\n' ...
    '            \\textbf{X} & %.3f~cm \\\\\n' ...
    '            \\textbf{Y} & %.3f~cm \\\\\n' ...
    '            \\textbf{Z} & %.3f~cm \\\\\n' ...
    '            \\bottomrule\n' ...
    '        \\end{tabular}\n' ...
    '        \\caption{Translation RMSE Metrics}\n' ...
    '        \\label{tab:pose_eval_translation}\n' ...
    '\\end{table}\n'], ...
    rmse*100, rmse_xyz(1), rmse_xyz(2), rmse_xyz(3));

% === Rotation RMSE table ===
latex_tbl_rot = sprintf([ ...
    '\\begin{table}[H]\n' ...
    '    \\centering\n' ...
    '        \\begin{tabular}{cc}\n' ... % <-- both columns centered
    '            \\toprule\n' ...
    '            \\textbf{Rotation Metric} & \\textbf{RMSE [$^{\\circ}$]} \\\\\n' ...
    '            \\midrule\n' ...
    '            \\textbf{Yaw} & %.3f$^{\\circ}$ \\\\\n' ...
    '            \\textbf{Pitch} & %.3f$^{\\circ}$ \\\\\n' ...
    '            \\textbf{Roll} & %.3f$^{\\circ}$ \\\\\n' ...
    '            \\textbf{Drift} & %.3f$^{\\circ}$ \\\\\n' ...
    '            \\bottomrule\n' ...
    '        \\end{tabular}\n' ...
    '        \\caption{Rotation RMSE Metrics}\n' ...
    '        \\label{tab:pose_eval_rotation}\n' ...
    '\\end{table}\n'], ...
    rpy_rmse(1), rpy_rmse(2), rpy_rmse(3), sqrt(mean(ang_diff.^2)));

% Print both tables in order
fprintf('\n%s\n%s\n', latex_tbl_trans, latex_tbl_rot);

% Optional: save to a .tex file with both tables
fid = fopen('pose_metrics_tables.tex','w');
if fid ~= -1
    fprintf(fid, '%s\n%s', latex_tbl_trans, latex_tbl_rot);
    fclose(fid);
end


%%

saveAllFiguresToPath('Z:\Datasets\Images\ELOPE\0016')

%% === Umeyama Alignment Function ===

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

function saveAllFiguresToPath(userDefinedPath)
% SAVEALLFIGURESTOPATH Saves all open figures to a user-defined path using the plot title as filename.
%
%   saveAllFiguresToPath(userDefinedPath) saves all open MATLAB figures to the
%   specified userDefinedPath. Filenames are based on the axes title (from title()).
%   Spaces are replaced with underscores. Falls back to 'Figure_N' if no title is found.

% Input Validation
if nargin == 0
    error('User-defined path is required.');
end

if ~ischar(userDefinedPath) || isempty(userDefinedPath)
    error('User-defined path must be a non-empty string.');
end

if ~exist(userDefinedPath, 'dir')
    mkdir(userDefinedPath);
    disp(['Created directory: ' userDefinedPath]);
end

% List all open figures
openFigures = findall(0, 'Type', 'figure');

if isempty(openFigures)
    disp('No open figures found.');
    return
end

% Save figures
for i = 1:numel(openFigures)
    currentFigure = openFigures(i);
    axesHandles = findall(currentFigure, 'Type', 'axes');
    
    % Default fallback name
    fileName = sprintf('Figure_%d', i);
    
    % Try to extract the title from the first axes
    if ~isempty(axesHandles)
        titleText = get(get(axesHandles(1), 'Title'), 'String');
        
        if ischar(titleText) && ~isempty(strtrim(titleText))
            % Replace spaces with underscores and strip illegal characters
            fileName = regexprep(strtrim(titleText), '\s+', '_');
            fileName = regexprep(fileName, '[^\w]', '_');  % Replace non-word chars
        end
    end

    % Construct full path
    fullFilePath = fullfile(userDefinedPath, [fileName '.png']);

    % Save the figure
    try
        saveas(currentFigure, fullFilePath);
        disp(['Figure ' num2str(i) ' saved to: ' fullFilePath]);
    catch ME
        disp(['Error saving figure ' num2str(i) ': ' ME.message]);
    end
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
    if rank(D) < m - 1
        error('Degenerate covariance rank, Umeyama alignment is not possible');
    end
    
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