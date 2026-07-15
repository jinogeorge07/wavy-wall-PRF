% Input Output Post Processing (Wavy Wall) — HillP-style formatting
clc; clear; close all;

%% ------------------------------------------------------------------------
% Paths
% E:\dedalus\Wavy Wall\flexible signed distance function\optimum_vpm_symmetric_wavywall\epsilon 0.25\dedalus_22859770_Re500\snapshots_channel\mean_v_Re500.00_c1.00

folderpath = './'

addpath('/gpfs/homefs1/jig23007/ray_tracing_matlab/cbrewer2-master');

%% ------------------------------------------------------------------------
% Wavy wall + domain parameters  (MATCH PYTHON: y in [-Ly/2, Ly/2])
Ny = 148; Nx = 108; Nz = 2; % Nz to Nx for streamwise terms
N = Ny*Nx;

epsilon = 0.12;

Lx = 2.4;
Ly = 2.0 + 1.1*2*epsilon;      % total height (python: [-Ly/2, Ly/2])
ymin = -Ly/2;
ymax =  Ly/2;

y0 = 1.0;
A1 = epsilon;
A2 = epsilon;
Re = 500;
%Re = 235;

quiver = "OFF";

%% ------------------------------------------------------------------------
% Load results (adjust filename as you use)
% Read the 3 files directly and stack the already-saved fields
% fname1 = 'stability_results_Re500_168x120_kz11_c1.mat';
% fname2 = 'stability_results_Re500_168x120_kz21_c1.mat';
% fname3 = 'stability_results_Re500_168x120_kz33_c1.mat';


fnames = {
        
        'stability_results_Re500_kz5_148x108.mat'
        'stability_results_Re500_kz20_148x108.mat'
        'stability_results_Re500_kz32_148x108.mat'

}

% S1 = load(fname1);
% S2 = load(fname2);
% S3 = load(fname3);


S = cell(size(fnames));   % storage for loaded structs
kz_idx = zeros(size(fnames));

for i = 1:numel(fnames)
    fname = fnames{i};
    S{i} = load(fname);

    % extract kz index
    tokens = regexp(fname, 'kz(\d+)', 'tokens');
    kz_idx(i) = str2double(tokens{1});

    % extract the actual kz value from the loaded struct
    kz_val(i) = S{i}.kz_list(kz_idx(i));
end

%% take x, y, kz info from one file
x  = S{1}.x;
y  = S{1}.y;

c_index = 1
% kz_index = [kz kz_list(21) kz_list(33)];
nc = numel(kz_val);

% pick only the desired c-slice from each file
U_hat  = cat(3, S{1}.U_hat(:,:,1),  S{2}.U_hat(:,:,1),  S{3}.U_hat(:,:,1));
V_hat  = cat(3, S{1}.V_hat(:,:,1),  S{2}.V_hat(:,:,1),  S{3}.V_hat(:,:,1));
W_hat  = cat(3, S{1}.W_hat(:,:,1),  S{2}.W_hat(:,:,1),  S{3}.W_hat(:,:,1));

U2_hat = cat(3, S{1}.U2_hat(:,:,1), S{2}.U2_hat(:,:,1), S{3}.U2_hat(:,:,1));
V2_hat = cat(3, S{1}.V2_hat(:,:,1), S{2}.V2_hat(:,:,1), S{3}.V2_hat(:,:,1));
W2_hat = cat(3, S{1}.W2_hat(:,:,1), S{2}.W2_hat(:,:,1), S{3}.W2_hat(:,:,1));

c_number = 12; 
% Fine grid mean flow etc.
load data_x.mat
load data_y.mat
load u_mean_zt.mat
load x_solid.mat
load y_solid.mat

%mask_smooth_wavywall_Re500_148x108_epsilon_0.25
filename_mask = ['mask_smooth_wavywall_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '_epsilon_' num2str(epsilon) '.mat'];
load(filename_mask)   % expects mask_smooth (same name as in hillp code)

% ---- Solid mask (robust) ----
% If your mask uses 100 in solid, this still works.
solid_mask = (mask_smooth == Re);

% ---- Interpolate mean flow onto coarse grid (x,y from stability file) ----
[X_old,Y_old] = meshgrid(data_x, data_y);
[Xc, Yc]      = meshgrid(x, y);

Umean_coarse  = interp2(X_old, Y_old, u_mean_zt, Xc, Yc, 'linear');
clear data_x data_y u_mean_zt X_old Y_old Xc Yc

%% ------------------------------------------------------------------------
% Wavy-wall curves on the SAME x-grid as contour plots
x_wall  = x(:).';   % row vector
y1_wall = -y0 - A1*sin(2*pi/Lx * x_wall);   % bottom wall
y2_wall =  y0 + A2*sin(2*pi/Lx * x_wall);   % top wall

%% ------------------------------------------------------------------------
% Make subfolder for snapshots
snapdir = fullfile(folderpath, sprintf('snapshots_wavywall_Re%d_c%d_%dx%d',Re, c_index, Nx,Ny));
if ~exist(snapdir,'dir'), mkdir(snapdir); end

hasBWR = exist('bluewhitered','file');

% ---- Layout + tick style (from hillp script) ----
BIG_TICKS = 40;                     % axis tick fontsize
NXT       = 4;                      % number of major x ticks
NYT       = 4;                      % number of major y ticks

FIG_POS   = [100 100 1350 1000];    % figure size in pixels
AX_POS    = [0.18 0.26 0.55 0.62];  % axes position [left bottom width height]
CB_POS    = [0.78 0.16 0.05 0.78];  % colorbar position

idx_list = 1:nc;

    Umax_common  = 0; U2max_common = 0;
    Vmax_common  = 0; V2max_common = 0;
    Wmax_common  = 0; W2max_common = 0;

    for ii = idx_list
        tmp = U_hat(:,:,ii);   tmp(solid_mask) = NaN;
        Umax_common = max(Umax_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = U2_hat(:,:,ii);  tmp(solid_mask) = NaN;
        U2max_common = max(U2max_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = V_hat(:,:,ii);   tmp(solid_mask) = NaN;
        Vmax_common = max(Vmax_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = V2_hat(:,:,ii);  tmp(solid_mask) = NaN;
        V2max_common = max(V2max_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = W_hat(:,:,ii);   tmp(solid_mask) = NaN;
        Wmax_common = max(Wmax_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = W2_hat(:,:,ii);  tmp(solid_mask) = NaN;
        W2max_common = max(W2max_common, max(abs(tmp(:)), [], 'omitnan'));
    end

    
%% ------------------------------------------------------------------------
% Main loop (matches hillp plotting style)
for i = 1:nc % min(12, size(U_hat,3))

    % ===== U: response =====
    %fU1 = figure('Visible','off','Position',FIG_POS);
    fU1 = figure('Visible','off','Position',[100 100 1000 800]);
    U_plot = U_hat(:,:,i);
    U_plot(solid_mask) = NaN;

    contourf(x, y, U_plot, 40, 'LineWidth', 0.5);
    set(gcf,'Renderer','painters');
    set(gca,'Color','white');
    shading flat
    hold on
    %caxis([-max(abs(U_plot(:))) max(abs(U_plot(:)))]);
    %clim = max(abs(U_plot(:)), [], 'omitnan');
    %caxis([-clim clim]);
    caxis([-Umax_common Umax_common]);

    % Fill solid regions (below y1 and above y2) in the CURRENT axes
    fill_patch_wavywall(x_wall, y1_wall, y2_wall);

    % Draw walls
    plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
    plot(x_wall, y2_wall, 'k', 'LineWidth', 2);

    % Apply hillp-style ticks/positions/colorbar
    wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, Umax_common);
    %title(sprintf('U Response mode, c=%.2f, kz=%.2f', c_list(i), kz), 'FontSize', 24);
    daspect([1 1 1])
    exportgraphics(fU1, fullfile(snapdir, sprintf('U_response_Re%d_kz%g_c%d.png',Re, kz_val(i), c_index)), 'Resolution', 200);
    close(fU1);

    % ===== X: forcing (U2_hat) =====
    %fU3 = figure('Visible','off','Position',FIG_POS);
    fU3 = figure('Visible','off','Position',[100 100 1000 800]);

    U2_plot = U2_hat(:,:,i);
    U2_plot(solid_mask) = NaN;

    contourf(x, y, U2_plot, 40, 'LineWidth', 0.5);
    set(gcf,'Renderer','painters');
    set(gca,'Color','white');
    shading flat
    hold on
    caxis([-U2max_common U2max_common]);

    fill_patch_wavywall(x_wall, y1_wall, y2_wall);
    plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
    plot(x_wall, y2_wall, 'k', 'LineWidth', 2);

    wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, U2max_common);
    daspect([1 1 1])
    exportgraphics(fU3, fullfile(snapdir, sprintf('X_forcing_Re%d_kz%g_c%d.png',Re, kz_val(i), c_index)), 'Resolution', 200);
    close(fU3);

    % ===== V: response =====
    %fV1 = figure('Visible','off','Position',FIG_POS);
    fV1 = figure('Visible','off','Position',[100 100 1000 800]);

    V_plot = V_hat(:,:,i);
    V_plot(solid_mask) = NaN;
    contourf(x, y, V_plot, 40, 'LineWidth', 0.5);
    set(gcf,'Renderer','painters');
    set(gca,'Color','white');
    shading flat
    caxis([-Vmax_common Vmax_common]);

    fill_patch_wavywall(x_wall, y1_wall, y2_wall);
    plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
    plot(x_wall, y2_wall, 'k', 'LineWidth', 2);

    wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, Vmax_common);

    %title(sprintf('V Response mode, c=%.2f, kz=%.2f', c_list(i), kz), 'FontSize', 24);
    daspect([1 1 1])
    exportgraphics(fV1, fullfile(snapdir, sprintf('V_response_Re%d_kz%g_c%d.png',Re, kz_val(i), c_index)), 'Resolution', 200);
    close(fV1);

    % ===== Y: forcing (V2_hat) =====
    %fV3 = figure('Visible','off','Position',FIG_POS);
    fV3 = figure('Visible','off','Position',[100 100 1000 800]);

    V2_plot = V2_hat(:,:,i);
    V2_plot(solid_mask) = NaN;

    contourf(x, y, V2_plot, 40, 'LineWidth', 0.5);
    set(gcf,'Renderer','painters');
    set(gca,'Color','white');
    shading flat
    caxis([-V2max_common V2max_common]);

    fill_patch_wavywall(x_wall, y1_wall, y2_wall);
    plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
    plot(x_wall, y2_wall, 'k', 'LineWidth', 2);

    wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, V2max_common);

    %title(sprintf('Y Forcing mode, c=%.2f, kz=%.2f', c_list(i), kz), 'FontSize', 24);
    daspect([1 1 1])
    exportgraphics(fV3, fullfile(snapdir, sprintf('Y_forcing_Re%d_kz%g_c%d.png',Re, kz_val(i), c_index)), 'Resolution', 200);
    close(fV3);

    % ===== W: response =====
    %fW1 = figure('Visible','off','Position',FIG_POS);
    fW1 = figure('Visible','off','Position',[100 100 1000 800]);

    W_plot = W_hat(:,:,i);
    W_plot(solid_mask) = NaN;

    contourf(x, y, W_plot, 40, 'LineWidth', 0.5);
    set(gcf,'Renderer','painters');
    set(gca,'Color','white');
    shading flat
    caxis([-Wmax_common Wmax_common]);

    % Fill solid regions
    fill_patch_wavywall(x_wall, y1_wall, y2_wall);

    % Draw walls
    plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
    plot(x_wall, y2_wall, 'k', 'LineWidth', 2);

    wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, Wmax_common);

    %title(sprintf('W Response mode, c=%.2f, kz=%.2f', c_list(i), kz), 'FontSize', 24);
    daspect([1 1 1])
    exportgraphics(fW1, fullfile(snapdir, sprintf('W_response_Re%d_kz%g_c%d.png', Re, kz_val(i), c_index)), 'Resolution', 300);
    close(fW1);


    % ===== Z: forcing (W2_hat) =====
    %fW3 = figure('Visible','off','Position',FIG_POS);
    fW3 = figure('Visible','off','Position',[100 100 1000 800]);

    W2_plot = W2_hat(:,:,i);
    W2_plot(solid_mask) = NaN;

    contourf(x, y, W2_plot, 40, 'LineWidth', 0.5);
    set(gcf,'Renderer','painters');
    set(gca,'Color','white');
    shading flat
    caxis([-W2max_common W2max_common]);

    fill_patch_wavywall(x_wall, y1_wall, y2_wall);
    plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
    plot(x_wall, y2_wall, 'k', 'LineWidth', 2);

    wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, W2max_common);

    %title(sprintf('W Forcing mode, c=%.2f, kz=%.2f', c_list(i), kz), 'FontSize', 24);
    daspect([1 1 1])
    exportgraphics(fW3, fullfile(snapdir, sprintf('Z_forcing_Re%d_kz%d_c%d.png',Re, kz_val(i), c_index)), 'Resolution', 200);
    close(fW3);

end

%% ------------------------------------------------------------------------
% Local helpers (no separate .m files needed)

function fill_patch_wavywall(x_wall, y1_wall, y2_wall)
% Fill below bottom wall and above top wall using current axis limits.
ax = gca;
hold(ax,'on');

ymin = ax.YLim(1);
ymax = ax.YLim(2);

xcol = x_wall(:);
y1   = y1_wall(:);
y2   = y2_wall(:);

solid_color = [0.30 0.30 0.30];
solid_alpha = 1.0;

% Bottom region: [ymin, y1(x)]
xb = [xcol; flipud(xcol)];
yb = [y1;   ymin*ones(size(xcol))];
patch(xb, yb, solid_color, 'EdgeColor','none', 'FaceAlpha',solid_alpha, 'Parent',ax);

% Top region: [y2(x), ymax]
xt = [xcol; flipud(xcol)];
yt = [ymax*ones(size(xcol)); flipud(y2)];
patch(xt, yt, solid_color, 'EdgeColor','none', 'FaceAlpha',solid_alpha, 'Parent',ax);
end

function wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, clim_val)
% HillP-style axis/tick/colorbar formatting, but for wavy wall.
ax = gca;

xlabel('x','FontSize',32);
ylabel('y','FontSize',32);

% Force colorbar limits symmetric about zero - THIS MUST BE BEFORE COLORBAR
caxis(ax, [-clim_val clim_val]);

% Now create the colorbar (it will inherit the symmetric limits)
c = colorbar;
c.Ticks = linspace(-clim_val, clim_val, 7);
c.TickLabels = strip_zeros(compose('%.3f', c.Ticks));
c.FontSize = BIG_TICKS;

% IMPORTANT: Ensure the colormap is applied to this figure
colormap(ax, bluewhitered);  % Re-apply colormap to ensure it's using bluewhitered

% Axis styling
ax.FontSize = BIG_TICKS;
ax.LineWidth = 1.6;
ax.TickLength = [0.02 0.02];
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;

% Ticks (from beginning to end)
xmin = min(x(:));
ax.XLim = [xmin Lx];
ax.XTick = linspace(xmin, Lx, NXT);
ax.XTickLabel = strip_zeros(compose('%.2f', ax.XTick));

ymin = min(y(:));
ymax = max(y(:));
ax.YLim = [ymin ymax];
ax.YTick = linspace(ymin, ymax, NYT);
ax.YTickLabel = strip_zeros(compose('%.2f', ax.YTick));

daspect([1 1 1]);

% Positioning
set(ax, 'Units','normalized', 'Position', AX_POS);
set(c, 'Units','normalized', 'Position', CB_POS);
end


function out = strip_zeros(lbls)
% Remove trailing zeros so labels have <= 2 decimals but look clean (1 not 1.00).
out = cellstr(lbls);
for k = 1:numel(out)
    s = out{k};
    s = regexprep(s,'(\.\d*?)0+$','$1');  % drop trailing zeros
    s = regexprep(s,'\.$','');            % drop trailing dot
    out{k} = s;
end
end
