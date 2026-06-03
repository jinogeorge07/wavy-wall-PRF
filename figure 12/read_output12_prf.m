% Input Output Post Processing for wavy wall
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
folderpath = './'

%% Periodic hill and Domain Properties
Ny = 140; Nx = 100; Nz = 2; % Nz to Nx for streamwise terms
N = Ny*Nx;

epsilon = 0.25;
h = 1;
y0 = h;
A1 = epsilon;
A2 = epsilon;
%Lx = 0.6 * pi;
Lx = 2.4;
Ly = 2.0 + 1.1*2*epsilon;      % total height (python: [-Ly/2, Ly/2])
ymin = -Ly/2;
ymax =  Ly/2;
Re = 235
x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only

% figure;
% plot(x_val,y1_vals)
% hold on
% quiver = "OFF";
kz_idx = [5,20,32];
% load('stability_results_Re100_120x96_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
% load('stability_results_Re100_120x96_kz20.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
% load('stability_results_Re100_120x96_kz32.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega


show_velocity = "true";
c_number = 12; % 96 % MULTIPLE OF NUMBER OF CORES
% ---- choose the kz index you ran (adjust if needed) ----
% <--- set to the i you used in the run

load data_x.mat
load data_y.mat
load u_mean_zt.mat               % fine-grid mean U from Dedalus
load x_solid.mat
load y_solid.mat
filename_mask = ['mask_smooth_wavywall_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '_epsilon_' num2str(epsilon) '.mat'];
load(filename_mask)
rwb1 = bluewhitered
%load mask_smooth_hillperiodic_Re100_120x96.mat

solid_mask = (mask_smooth == Re);   % <<< IMPORTANT: > 0, not == 0

%guards for walls & colormap
haveWalls = exist('x_val','var') && exist('y1_vals','var') && exist('y2_vals','var');
hasBWR   = exist('bluewhitered','file');

% % ---- layout + tick style (same idea as working case) ----
% BIG_TICKS = 40;                         % axis tick fontsize
% NXT       = 4;                          % number of major x ticks
% NYT       = 4;                          % number of major y ticks
% 
% FIG_POS   = [100 100 1350 1000];        % figure size in pixels
% AX_POS  = [0.12 0.28 0.78 0.58];
% CB_POS = [0.18 0.20 0.64 0.04];

% ---- layout + tick style (same idea as working case) ----
BIG_TICKS = 40;                     % axis tick fontsize
NXT       = 4;                      % number of major x ticks
NYT       = 4;                      % number of major y ticks

FIG_POS   = [100 100 1350 1000];    % figure size in pixels
AX_POS    = [0.18 0.26 0.55 0.62];  % axes position [left bottom width height]
CB_POS    = [0.78 0.16 0.05 0.78];  % colorbar position

kz_idx_list = [5,20,32]
i_fixed = 1;

Umax_common = 0;
Vmax_common = 0;
Wmax_common = 0;

%stability_results_Re235_kz10_140x100
for kz_idx = kz_idx_list
    fname = ['stability_results_Re' num2str(Re) '_kz' num2str(kz_idx) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
    load(fname);

    tmp = U2_hat(:,:,i_fixed);
    tmp(solid_mask) = NaN;
    Umax_common = max(Umax_common, max(abs(tmp(:)), [], 'omitnan'));

    tmp = V2_hat(:,:,i_fixed);
    tmp(solid_mask) = NaN;
    Vmax_common = max(Vmax_common, max(abs(tmp(:)), [], 'omitnan'));

    tmp = W2_hat(:,:,i_fixed);
    tmp(solid_mask) = NaN;
    Wmax_common = max(Wmax_common, max(abs(tmp(:)), [], 'omitnan'));
end


%% ------------------------------------------------------------------------
% Wavy-wall curves on the SAME x-grid as contour plots
x_wall  = x(:).';   % row vector
y1_wall = -y0 - A1*sin(2*pi/Lx * x_wall);   % bottom wall
y2_wall =  y0 + A2*sin(2*pi/Lx * x_wall);   % top wall
%% ------------------------------------------------------------------------

for kz_idx = [5,20,32]

    fname = ['stability_results_Re' num2str(Re) '_kz' num2str(kz_idx) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
    load(fname);

    kz = kz_list(kz_idx);

    [X_old,Y_old] = meshgrid(data_x, data_y);
    [Xc, Yc]      = meshgrid(x, y);  % 'x','y' came from stability_results.mat

    Umean_coarse  = interp2(X_old, Y_old, u_mean_zt, Xc, Yc, 'linear');  % note transpose
    [dUdx_coarse, dUdy_coarse] = gradient(Umean_coarse, x, y);

    % (optional) clear big vars:
    %clear u_mean_zt X_old Y_old Xc Yc

    % grid for quiver (same grid you used for pcolor/contour)
    [Xg, Yg] = meshgrid(x, y);

    % thin the arrows a bit so the plot stays readable
    qstep = max(1, round(min(numel(x), numel(y))/40));   % ~25–35 arrows each way
    Xs = Xg(1:qstep:end, 1:qstep:end);
    Ys = Yg(1:qstep:end, 1:qstep:end);

    for i = [1]

        crit_contour = Umean_coarse - c_list(i);   % Ny x Nx

        %     % ---------- U: contour forcing mode ----------
        fU3 = figure('Visible','on','Position',[100 100 1000 800]);
        %contourf(x, y, U2_hat(:,:,i), 30, 'LineWidth', 1/2);
        % Copy U and blank out the solid region
        U2_plot = U2_hat(:,:,i);
        U2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, U2_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        %caxis([-max(abs(U2_plot(:))) max(abs(U2_plot(:)))]);  % ensures 0 is centered
        caxis([-Umax_common Umax_common]);
        hold on

        writematrix(U2_plot, fullfile(sprintf('X_forcing_c%02d_kz%g.csv', i, kz)));
        % fill solid patch
        fill_patch_wavywall(x_wall, y1_wall, y2_wall);
        plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
        plot(x_wall, y2_wall, 'k', 'LineWidth', 2);
        hold on

        wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, Umax_common);
        daspect([1 1 1]);
        exportgraphics(fU3, sprintf('X_forcing_c%02d_kz%g.png', i, kz), ...
        'Resolution', 200, 'BackgroundColor', 'white');
        close(fU3);

        %% ---------- V: contour forcing mode----------
        fV3 = figure('Visible','on','Position',[100 100 1000 800]);
        %contourf(x, y, V2_hat(:,:,i), 30, 'LineWidth', 1/2);

        %% Copy U and blank out the solid region
        V2_plot = V2_hat(:,:,i);
        V2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, V2_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        %caxis([-max(abs(V2_plot(:))) max(abs(V2_plot(:)))]);  % ensures 0 is centered
        caxis([-Vmax_common Vmax_common]);
        hold on

        writematrix(V2_plot, fullfile(sprintf('Y_forcing_c%02d_kz%g.csv', i, kz)));
        % fill solid patch
        fill_patch_wavywall(x_wall, y1_wall, y2_wall);
        plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
        plot(x_wall, y2_wall, 'k', 'LineWidth', 2);
        hold on

        wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, Vmax_common);
        daspect([1 1 1]);
        exportgraphics(fV3, sprintf('Y_forcing_c%02d_kz%g.png', i, kz), ...
        'Resolution', 200, 'BackgroundColor', 'white');
        close(fV3);

        %     % ---------- Z: contour Forcing mode ----------
            fW3 = figure('Visible','off','Position',[100 100 1000 800]);
            %contourf(x, y, W2_hat(:,:,i), 30, 'LineWidth', 1/2);
            %% Copy U and blank out the solid region
            W2_plot = W2_hat(:,:,i);
            W2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
            contourf(x, y, W2_plot, 40, 'LineWidth', 1/2);
            shading interp
            colormap(bluewhitered);                     % ensures the center color is white
            %caxis([-max(abs(W2_plot(:))) max(abs(W2_plot(:)))]);  % ensures 0 is centered
            caxis([-Wmax_common Wmax_common]);
            hold on
        
            writematrix(W2_plot, fullfile(sprintf('Z_forcing_c%02d_kz%g.csv', i, kz)));
            % fill solid patch
            fill_patch_wavywall(x_wall, y1_wall, y2_wall);
            plot(x_wall, y1_wall, 'k', 'LineWidth', 2);
            plot(x_wall, y2_wall, 'k', 'LineWidth', 2);
            hold on

            wavywall_tick_function(x, y, Lx, BIG_TICKS, NXT, NYT, AX_POS, CB_POS, Wmax_common);
            daspect([1 1 1]);
            exportgraphics(fW3, sprintf('Z_forcing_c%02d_kz%g.png', i, kz), ...
            'Resolution', 200, 'BackgroundColor', 'white');
            %%saveas(fW3, fullfile(sprintf('Z_forcing_c%02d_kz%g.png', i, kz)));
            close(fW3);

    end

end

% sigma_list = result_sigma(:,kz_idx)
% filename = ['sigma_list_Re' num2str(Re) '_kz' num2str(kz_idx) '.mat'];
% save(filename, 'sigma_list');

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

% Force colorbar limits symmetric about zero

c = colorbar;
c.Ticks      = linspace(-clim_val, clim_val, 7);
%c.Ticks      = linspace(c.Limits(1), c.Limits(2), 6);
c.TickLabels = strip_zeros(compose('%.3f', c.Ticks));
c.FontSize   = BIG_TICKS;

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
ax.XLim  = [xmin Lx];
ax.XTick = linspace(xmin, Lx, NXT);
ax.XTickLabel = strip_zeros(compose('%.2f', ax.XTick));

ymin = min(y(:));
ymax = max(y(:));
ax.YLim  = [ymin ymax];
ax.YTick = linspace(ymin, ymax, NYT);
ax.YTickLabel = strip_zeros(compose('%.2f', ax.YTick));

daspect([1 1 1]);
% Positioning
set(ax, 'Units','normalized', 'Position', AX_POS);
set(c,  'Units','normalized', 'Position', CB_POS);

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
