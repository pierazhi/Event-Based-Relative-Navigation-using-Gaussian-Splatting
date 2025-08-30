% ========================================================================
% Partial-Trajectory Alignment and Comparative Pose Evaluation
% Author: Piercarlo Fontana
%
% Description:
% This script compares a pipeline's estimated trajectory against ground
% truth and a COLMAP baseline over a guaranteed-overlap partial window.
% It trims sequences to a common index range, aligns each candidate to the
% ground-truth segment using Umeyama similarity (rotation, translation,
% and optional scale), and produces metrics and figures that actually let
% you see who's lying. After alignment, it extracts positions, computes
% framewise translation errors, filters COLMAP outliers with a configurable
% threshold, and reports both per-axis and total RMSE. For attitude, it
% performs a first-pose orientation alignment to isolate relative drift,
% then evaluates yaw–pitch–roll RMSE and end-to-end rotational drift using
% the geodesic angle. It generates publication-ready LaTeX tables to drop
% straight into a report and offers optional batch saving of figures.
%
% Inputs:
%   - Ground truth poses in KITTI format:            "kitti.txt"
%   - COLMAP trajectory in KITTI format:             "colmap_trajectory_kitti.txt"
%   - Pipeline estimate (partial) in KITTI format:   "pose_file_est_BA_f073_whole.txt"
%   All inputs are Nx12 row-wise SE(3) matrices; the loader reshapes them
%   to 4x4 per-frame homogeneous transforms.
%
% What it does, in order:
%   1) Chooses a partial window starting at 'partial_start_idx' and clamps
%      all sequences to the overlapping region so metrics aren't polluted
%      by missing frames.
%   2) Aligns COLMAP→GT and Pipeline→GT on the partial segment via Umeyama
%      to neutralize global offsets and scale mismatches.
%   3) Plots aligned 3D trajectories and per-axis position traces for a
%      quick sanity check.
%   4) Computes translation errors, applies an inlier mask to COLMAP based
%      on a distance threshold in meters, and prints thresholded RMSE.
%   5) Performs first-pose orientation alignment and reports RPY RMSE plus
%      end-to-end drift in degrees for both methods.
%   6) Assembles side-by-side LaTeX tables (booktabs) comparing Pipeline
%      and COLMAP on translation and rotation, clearly stating the inlier
%      policy used for COLMAP.
%
% Key parameters you may want to touch:
%   - partial_start_idx: starting frame (1-based) for the comparison window.
%   - thr_m: COLMAP inlier threshold in meters for translation RMSE.
%
% Outputs:
%   - Console printouts for RMSE per axis and total (cm), RPY RMSE (deg),
%     and rotational drift (deg).
%   - Figures named 'Aligned', 'Error', 'Per_Axis_Fixed', 'Error_Fixed',
%     and 'RPY' for trajectory and error visualization.
%   - A LaTeX block with two tables ready to paste into your paper or thesis.
%   - Optional PNG exports of all open figures if you call saveAllFiguresToPath.
%
% Caveats:
%   - If your inputs don't overlap, the guard throws an error. Fix your
%     indices or your data.
%   - COLMAP RMSE is reported on inliers only, by design; the threshold
%     communicates that policy in the LaTeX caption.
%
% ========================================================================

clearvars; clc; close all;

% === Settings ===
partial_start_idx = 1;    % starting frame for the partial window (1-based)
purple = [173/255, 7/255, 249/255];

% === GLOBAL FONT SETTINGS ===
AXES_FONT_NAME   = 'Helvetica';
AXES_FONT_SIZE   = 14; % ticks
LABEL_FONT_SIZE  = 14; % xlabel/ylabel/zlabel
LEGEND_FONT_SIZE = 12; % legends

% Helpers to apply fonts to current axes and text objects
applyFontSettings   = @(ax) set(ax, 'FontName', AXES_FONT_NAME, 'FontSize', AXES_FONT_SIZE, 'LineWidth', 1);
applyLabelSettings  = @(h)  set(h,  'FontName', AXES_FONT_NAME, 'FontSize', LABEL_FONT_SIZE);
applyLegendSettings = @(h)  set(h,  'FontName', AXES_FONT_NAME, 'FontSize', LEGEND_FONT_SIZE);

% === Load poses (use your own paths as needed) ===
gt_full     = load_kitti_poses("kitti.txt");     
colmap_full = load_kitti_poses("colmap_trajectory_kitti.txt");
est_partial = load_kitti_poses("pose_file_est_BA_f073_whole.txt");

% === Sizes and partial window ===
N_gt  = size(gt_full,     3);
N_col = size(colmap_full, 3);
N_est = size(est_partial, 3);

partial_end_idx = partial_start_idx + N_est - 1;
partial_end_idx = min([partial_end_idx, N_gt, N_col]);      % enforce overlap
N_part = partial_end_idx - partial_start_idx + 1;

% Guard if indices collapse
if N_part <= 0
    error('Partial window is empty. Check partial_start_idx and input sizes.');
end

% === Extract partial segments that are guaranteed to overlap ===
gt_part     = gt_full(:,:,partial_start_idx:partial_end_idx);
colmap_part = colmap_full(:,:,partial_start_idx:partial_end_idx);
est_partial = est_partial(:,:,1:N_part);  % trim in case it was longer

% === Raw (unaligned) positions for quick visualization (PARTIAL ONLY) ===
gt_xyz_part_raw     = squeeze(gt_full(1:3, 4, partial_start_idx:partial_end_idx))';
colmap_xyz_part_raw = squeeze(colmap_full(1:3, 4, partial_start_idx:partial_end_idx))';
est_xyz_part_raw    = squeeze(est_partial(1:3, 4, :))';

% figure('Name','Not_Aligned');
% plot3(gt_xyz_part_raw(:,1), gt_xyz_part_raw(:,2), gt_xyz_part_raw(:,3), 'k-', 'LineWidth', 2); hold on;
% plot3(est_xyz_part_raw(:,1), est_xyz_part_raw(:,2), est_xyz_part_raw(:,3), '--', 'Color', purple, 'LineWidth', 1.5);
% plot3(colmap_xyz_part_raw(:,1), colmap_xyz_part_raw(:,2), colmap_xyz_part_raw(:,3), 'c--', 'LineWidth', 1.5);
% hX = xlabel('X'); hY = ylabel('Y'); hZ = zlabel('Z'); grid on; axis equal;
% hL = legend('GT', 'Pipeline', 'Colmap', 'Location', 'Best');
% applyFontSettings(gca); applyLabelSettings([hX hY hZ]); applyLegendSettings(hL);

% === Align COLMAP (partial) -> GT (partial) ===
gt_xyz_part_3xN     = squeeze(gt_part(1:3, 4, :));        % 3 x N_part
colmap_xyz_part_3xN = squeeze(colmap_part(1:3, 4, :));    % 3 x N_part
[R_c, t_c, s_c] = umeyama_alignment(colmap_xyz_part_3xN, gt_xyz_part_3xN, true);

colmap_aligned_part = zeros(4,4,N_part);
for i = 1:N_part
    T = colmap_part(:,:,i);
    colmap_aligned_part(1:3,1:3,i) = s_c * R_c * T(1:3,1:3);
    colmap_aligned_part(1:3,4,i)   = s_c * R_c * T(1:3,4) + t_c;
    colmap_aligned_part(4,4,i)     = 1;
end

% === Align Pipeline (partial) -> GT (partial) ===
est_xyz_part_3xN = squeeze(est_partial(1:3, 4, :));       % 3 x N_part
[R_e, t_e, s_e] = umeyama_alignment(est_xyz_part_3xN, gt_xyz_part_3xN, true);

est_aligned_part = zeros(4,4,N_part);
for i = 1:N_part
    T = est_partial(:,:,i);
    est_aligned_part(1:3,1:3,i) = s_e * R_e * T(1:3,1:3);
    est_aligned_part(1:3,4,i)   = s_e * R_e * T(1:3,4) + t_e;
    est_aligned_part(4,4,i)     = 1;
end

% === Extract positions for plots (PARTIAL ONLY) ===
gt_xyz_part_vec       = squeeze(gt_part(1:3, 4, :))';                 
colmap_xyz_part_vec   = squeeze(colmap_aligned_part(1:3, 4, :))';     
est_xyz_part_vec      = squeeze(est_aligned_part(1:3, 4, :))';        

% === Plot (after alignment) ===
figure('Name','Aligned');
plot3(gt_xyz_part_vec(:,1), gt_xyz_part_vec(:,2), gt_xyz_part_vec(:,3), 'k-', 'LineWidth', 2); hold on;
plot3(colmap_xyz_part_vec(:,1), colmap_xyz_part_vec(:,2), colmap_xyz_part_vec(:,3), 'c--', 'LineWidth', 1.5);
plot3(est_xyz_part_vec(:,1), est_xyz_part_vec(:,2), est_xyz_part_vec(:,3), '-', 'Color', purple, 'LineWidth', 2);
hX = xlabel('X [m]'); hY = ylabel('Y [m]'); hZ = zlabel('Z [m]'); grid on; axis equal;
hL = legend('GT', 'Colmap', 'Pipeline', 'Location', 'Best');
applyFontSettings(gca); applyLabelSettings([hX hY hZ]); applyLegendSettings(hL); 

% === XYZ position evolution ===
xyz_labels = {'X [m]', 'Y [m]', 'Z [m]'};
% figure('Name', 'Per_Axis');
% for j = 1:3
%     subplot(3,1,j);
%     plot(partial_start_idx:partial_end_idx, gt_xyz_part_vec(:,j), 'k-', 'LineWidth', 2); hold on;
%     plot(partial_start_idx:partial_end_idx, est_xyz_part_vec(:,j), '--', 'Color', purple, 'LineWidth', 2);
%     plot(partial_start_idx:partial_end_idx, colmap_xyz_part_vec(:,j), 'c-.', 'LineWidth', 2);
%     hY = ylabel(xyz_labels{j});
% 
%     if j == 3
%         hX = xlabel('Frame Index');
%     else
%         hX = []; % nothing to set on rows 1-2
%     end
%     hL = legend('GT', 'Pipeline', 'COLMAP', 'Location', 'Best'); grid on;  
%     applyFontSettings(gca);
%     if ~isempty(hX), applyLabelSettings(hX); end
%     applyLabelSettings(hY);
%     if exist('hT','var') && ~isempty(hT) && j == 1,  end
%     applyLegendSettings(hL);
% end

% === Translation errors on partial window ===
err_colmap_part = sqrt(sum((gt_xyz_part_vec - colmap_xyz_part_vec).^2, 2));
err_est_part    = sqrt(sum((gt_xyz_part_vec - est_xyz_part_vec).^2,  2));

% === Remove COLMAP points with large deviation from GT ===
thr_m = 1;  % threshold in meters (e.g., 0.10 m = 10 cm)
inliers_colmap = err_colmap_part <= thr_m;

% Filter RMSE calculation
rmse_colmap_thr = sqrt(mean(err_colmap_part(inliers_colmap).^2)) * 100; % cm

fprintf('\n📏 Partial Translation RMSEs (threshold filtered):\n');
fprintf('  COLMAP (≤ %.0f cm): %.2f cm  [kept %d/%d frames]\n', ...
        thr_m*100, rmse_colmap_thr, nnz(inliers_colmap), numel(inliers_colmap));
fprintf('  Pipeline:           %.2f cm\n', sqrt(mean(err_est_part.^2)) * 100);

% Optional: plot with inliers only
figure('Name','Error');
idx = partial_start_idx:partial_end_idx;
plot(idx(inliers_colmap), err_colmap_part(inliers_colmap)*100, 'c-', 'LineWidth', 2); hold on;
plot(idx(~inliers_colmap), err_colmap_part(~inliers_colmap)*100, 'r.', 'MarkerSize', 12); % excluded
plot(idx, err_est_part*100, '-', 'Color', purple, 'LineWidth', 2);
hX = xlabel('Frame Index'); hY = ylabel('Translation Error [cm]'); grid on;
hL = legend('Colmap (inliers)', 'Colmap (excluded)', 'Pipeline', 'Location', 'Best');
applyFontSettings(gca); applyLabelSettings([hX hY]); applyLegendSettings(hL); 

figure('Name', 'Per_Axis_Fixed');
for j = 1:3
    subplot(3,1,j);
    plot(idx, gt_xyz_part_vec(:,j), 'k-', 'LineWidth', 2); hold on;
    plot(idx, est_xyz_part_vec(:,j), '--', 'Color', purple, 'LineWidth', 2);
    plot(idx(inliers_colmap), colmap_xyz_part_vec(inliers_colmap,j), 'c-.', 'LineWidth', 2);
    hY = ylabel(xyz_labels{j});
   
    if j == 3
        hX = xlabel('Frame Index');
    else
        hX = [];
    end
    hL = legend('GT', 'Pipeline', 'COLMAP', 'Location', 'Best'); grid on;
    applyFontSettings(gca);
    if ~isempty(hX), applyLabelSettings(hX); end
    applyLabelSettings(hY);
    if exist('hT','var') && ~isempty(hT) && j == 1 end
    applyLegendSettings(hL);
end

figure('Name','Error_Fixed');
plot(idx(inliers_colmap), err_colmap_part(inliers_colmap)*100, 'c--', 'LineWidth', 2); hold on;
plot(idx, err_est_part*100, '-', 'Color', purple, 'LineWidth', 2);
hX = xlabel('Frame Index'); hY = ylabel('Translation Error [cm]'); grid on;
hL = legend('Colmap', 'Pipeline', 'Location', 'best');
applyFontSettings(gca); applyLabelSettings([hX hY]); applyLegendSettings(hL); 

% === First-pose orientation alignment (rotation-only) for RPY plots — PARTIAL ===
R_gt0     = gt_part(1:3,1:3,1);
R_est0    = est_partial(1:3,1:3,1);
R_colmap0 = colmap_part(1:3,1:3,1);

R_align_est    = R_gt0 * R_est0';
R_align_colmap = R_gt0 * R_colmap0';

N_part = size(gt_part, 3);
rpy_gt     = zeros(N_part, 3);
rpy_est    = zeros(N_part, 3);
rpy_colmap = zeros(N_part, 3);

for i = 1:N_part
    Rg = gt_part(1:3,1:3,i);
    Re = R_align_est    * est_partial(1:3,1:3,i);
    Rc = R_align_colmap * colmap_part(1:3,1:3,i);
    % ZYX: yaw-pitch-roll
    rpy_gt(i,:)     = rad2deg(rotm2eul(Rg, 'ZYX'));
    rpy_est(i,:)    = rad2deg(rotm2eul(Re, 'ZYX'));
    rpy_colmap(i,:) = rad2deg(rotm2eul(Rc, 'ZYX'));
end

% unwrap each column independently, then wrap to [-180, 180] if desired
rpy_gt     = unwrap(rpy_gt);
rpy_est    = unwrap(rpy_est);
rpy_colmap = unwrap(rpy_colmap);

rpy_labels = {'Yaw', 'Pitch', 'Roll'};
figure('Name', 'RPY');
for j = 1:3
    subplot(3,1,j);
    plot(partial_start_idx:partial_end_idx, rpy_gt(:,j), 'k-', 'LineWidth', 2); hold on;
    plot(partial_start_idx:partial_end_idx, rpy_est(:,j), '--', 'Color', purple, 'LineWidth', 2);
    plot(partial_start_idx:partial_end_idx, rpy_colmap(:,j), 'c-.', 'LineWidth', 2);
    hY = ylabel([rpy_labels{j} ' [°]']);
 
    if j == 3
        hX = xlabel('Frame Index');
    else
        hX = [];
    end
    hL = legend('GT', 'Pipeline', 'Colmap', 'Location', 'Best'); grid on;
    applyFontSettings(gca);
    if ~isempty(hX), applyLabelSettings(hX); end
    applyLabelSettings(hY);
    if exist('hT','var') && ~isempty(hT) && j == 1,  end
    applyLegendSettings(hL);
end

% === METRICS + LATEX TABLES (Pipeline vs COLMAP) ===
% Helpers
rmse = @(x) sqrt(mean(x.^2,'omitnan'));
ang_from_rotm = @(R) rad2deg(acos(max(-1,min(1,(trace(R)-1)/2))));  % numerically safe

% --- Translation RMSEs (cm) ---
% Errors (m) on the overlapping, aligned partial window
e_est_xyz    = est_xyz_part_vec    - gt_xyz_part_vec;      % N x 3
e_colmap_xyz = colmap_xyz_part_vec - gt_xyz_part_vec;      % N x 3

% Per-axis RMSE (cm). For COLMAP, respect your inlier mask.
rmse_est_xyz_cm    = 100 * sqrt(mean(e_est_xyz.^2,           1, 'omitnan'));        % [X Y Z]
rmse_colmap_xyz_cm = 100 * sqrt(mean(e_colmap_xyz(inliers_colmap,:).^2, 1, 'omitnan'));

% Total translation RMSE (cm) using the 3D norm per frame
err_est_norm_cm    = 100 * sqrt(sum(e_est_xyz.^2,    2, 'omitnan'));
err_colmap_norm_cm = 100 * sqrt(sum(e_colmap_xyz.^2, 2, 'omitnan'));
rmse_est_total_cm    = rmse(err_est_norm_cm);
rmse_colmap_total_cm = rmse(err_colmap_norm_cm(inliers_colmap));

% --- Rotation RMSEs (deg) on ZYX angles (already unwrapped by you) ---
rpy_rmse_est_deg    = sqrt(mean((rpy_est    - rpy_gt).^2,    1, 'omitnan'));  % [Yaw Pitch Roll]
rpy_rmse_colmap_deg = sqrt(mean((rpy_colmap - rpy_gt).^2,    1, 'omitnan'));

% --- Orientation drift (deg) from first to last frame (after your first-pose alignment) ---
Rg_last = gt_part(1:3,1:3,end);
Re_last = (R_align_est    * est_partial(1:3,1:3,end));
Rc_last = (R_align_colmap * colmap_part(1:3,1:3,end));
drift_est_deg    = ang_from_rotm(Re_last * Rg_last');   % geodesic angle between end attitudes
drift_colmap_deg = ang_from_rotm(Rc_last * Rg_last');

%% === LaTeX: side-by-side tables (requires \usepackage{booktabs}) ===
latex_block = sprintf([ ...
'\\begin{minipage}{.5\\textwidth}\n' ...
'  \\begin{table}[H]\n' ...
'    \\centering\n' ...
'    \\begin{tabular}{lcc}\n' ...
'      \\toprule\n' ...
'      \\textbf{Translation Metric} & \\textbf{Pipeline [cm]} & \\textbf{COLMAP [cm]} \\\\\n' ...
'      \\midrule\n' ...
'      \\textbf{X}      & %.3f & %.3f \\\\\n' ...
'      \\textbf{Y}      & %.3f & %.3f \\\\\n' ...
'      \\textbf{Z}      & %.3f & %.3f \\\\\n' ...
'      \\midrule\n' ...
'      \\textbf{Total}  & %.3f & %.3f \\\\\n' ...
'      \\bottomrule\n' ...
'    \\end{tabular}\n' ...
'    \\caption{Orbit Ambient: Translation RMSE (COLMAP computed on inliers $\\leq$ %.0f cm)}\n' ...
'    \\label{tab:pose_eval_translation}\n' ...
'  \\end{table}\n' ...
'\\end{minipage}\n' ...
'\\hfill\n' ...
'\\begin{minipage}{.5\\textwidth}\n' ...
'  \\begin{table}[H]\n' ...
'    \\centering\n' ...
'    \\begin{tabular}{lcc}\n' ...
'      \\toprule\n' ...
'      \\textbf{Rotation Metric} & \\textbf{Pipeline [$^{\\circ}$]} & \\textbf{COLMAP [$^{\\circ}$]} \\\\\n' ...
'      \\midrule\n' ...
'      \\textbf{Yaw}   & %.3f & %.3f \\\\\n' ...
'      \\textbf{Pitch} & %.3f & %.3f \\\\\n' ...
'      \\textbf{Roll}  & %.3f & %.3f \\\\\n' ...
'      \\midrule\n' ...
'      \\textbf{Drift} & %.3f & %.3f \\\\\n' ...
'      \\bottomrule\n' ...
'    \\end{tabular}\n' ...
'    \\caption{Orbit Ambient: Rotation Metrics (RMSE per axis and end-to-end drift)}\n' ...
'    \\label{tab:pose_eval_rotation}\n' ...
'  \\end{table}\n' ...
'\\end{minipage}\n' ...
], ...
rmse_est_xyz_cm(1), rmse_colmap_xyz_cm(1), ...
rmse_est_xyz_cm(2), rmse_colmap_xyz_cm(2), ...
rmse_est_xyz_cm(3), rmse_colmap_xyz_cm(3), ...
rmse_est_total_cm,  rmse_colmap_total_cm,  thr_m*100, ...
rpy_rmse_est_deg(1),    rpy_rmse_colmap_deg(1), ...
rpy_rmse_est_deg(2),    rpy_rmse_colmap_deg(2), ...
rpy_rmse_est_deg(3),    rpy_rmse_colmap_deg(3), ...
drift_est_deg,          drift_colmap_deg);

fprintf('\n%s\n', latex_block);


%%

%saveAllFiguresToPath('C:\Users\Pierazhi\Documents\MATLAB\Tesi\ELOPE_images')

%% Function

function applyFontStyle(axHandle, axesFontName, axesFontSize, textFontName, textFontSize)
    set(axHandle, 'FontName', axesFontName, 'FontSize', axesFontSize, 'LineWidth', 1);
    txtObjects = findall(ancestor(axHandle, 'figure'), 'Type', 'text');
    set(txtObjects, 'FontName', textFontName, 'FontSize', textFontSize);
    lgd = findobj(ancestor(axHandle, 'figure'), 'Type', 'Legend');
    set(lgd, 'FontName', textFontName, 'FontSize', textFontSize);
end

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

function saveAllFiguresToPath(userDefinedPath)
% SAVEALLFIGURESTOPATH Saves all open figures to a user-defined path using the figure name as filename.
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

