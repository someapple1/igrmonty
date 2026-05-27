clear all;
clf reset;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool("Threads");
end

archive_dir = fullfile('E:\Raytracing\igrmonty', [char([24402 26723]), '(1)']);
old_dir = pwd;
cleanup_obj = onCleanup(@() cd(old_dir));
cd(archive_dir);

% USER DIAGNOSTIC PATCH: merge a wider viewing-angle range to reduce
% Monte Carlo noise in the high-frequency scattered tail.
theta_center_deg = 45;
theta_half_width_deg = 10;
theta_min = deg2rad(theta_center_deg - theta_half_width_deg);
theta_max = deg2rad(theta_center_deg + theta_half_width_deg);


for kk = 0:300
output = load(['SgrA_forward_data/output', num2str(kk), '.mat']).output;


obs_theta = output.theta_c;
nu_left = output.nu_left;
obs_nu_c = output.nu_c;


% USER DIAGNOSTIC PATCH: use overlap-weighted dcos(theta) averaging over
% [theta_min, theta_max], instead of a narrow fixed number of bins.
theta_left = obs_theta(1:end-1);
theta_right = obs_theta(2:end);
overlap_left = max(theta_left, theta_min);
overlap_right = min(theta_right, theta_max);
theta_mask = overlap_right > overlap_left;
quanzhong = abs(cos(overlap_left(theta_mask)) ...
    - cos(overlap_right(theta_mask)));
quanzhong = quanzhong/sum(quanzhong);
quanzhong = reshape(quanzhong, 1, []);

p_nu_L_nu_all(kk+1,:) = sum(output.res(:, theta_mask).*quanzhong, 2);

 p_nu_L_nu_0(kk+1,:) = sum(output.res0(:, theta_mask).*quanzhong, 2);
 p_nu_L_nu_1(kk+1,:) = sum(output.res1(:, theta_mask).*quanzhong, 2);

 p_nu_L_nu_2(kk+1,:) = sum(output.res2(:, theta_mask).*quanzhong, 2);
 p_nu_L_nu_3(kk+1,:) = sum(output.res3(:, theta_mask).*quanzhong, 2);

end
nu_Lnu_all = mean(p_nu_L_nu_all);
nu_Lnu_0 = mean(p_nu_L_nu_0);

nu_Lnu_1 = mean(p_nu_L_nu_1);

nu_Lnu_2 = mean(p_nu_L_nu_2);
nu_Lnu_3 = mean(p_nu_L_nu_3);



%% spectrum.h5, same processing style as tools/plspec.py
h5_path = fullfile('E:\Raytracing\igrmonty\igrmonty', 'spectrum.h5');
ME = 9.1093897e-28;
CL = 2.99792458e10;
HPL = 6.6260755e-27;
LSUN = 3.827e33;

h5_lnu = h5read(h5_path, '/output/lnu');
h5_nu = (10.^h5_lnu) * ME * CL * CL / HPL;
h5_nuLnu = h5read(h5_path, '/output/nuLnu') * LSUN;

% In Python/h5py this dataset is [component, frequency, theta].
% MATLAB h5read returns it as [theta, frequency, component].
% Match the original forward-data angular average with the same cos(theta)
% weighting over the selected theta interval.
% USER DIAGNOSTIC PATCH: use the same widened theta interval for h5.
h5_theta_min = theta_min;
h5_theta_max = theta_max;
h5_nth = size(h5_nuLnu, 1);
h5_dtheta = (pi/2)/h5_nth;
h5_nuLnu_theta = zeros(size(h5_nuLnu, 2), size(h5_nuLnu, 3));
h5_weight_sum = 0;

for h5_j = 1:h5_nth
    bin_left = (h5_j - 1)*h5_dtheta;
    bin_right = h5_j*h5_dtheta;
    overlap_left = max(bin_left, h5_theta_min);
    overlap_right = min(bin_right, h5_theta_max);

    if overlap_right > overlap_left
        h5_weight = abs(cos(overlap_left) - cos(overlap_right));
        h5_weight_sum = h5_weight_sum + h5_weight;
        h5_nuLnu_theta = h5_nuLnu_theta + h5_weight*squeeze(h5_nuLnu(h5_j,:,:));
    end
end

h5_nuLnu_theta = h5_nuLnu_theta/h5_weight_sum; % [frequency, component]
h5_nuLnu_total = sum(h5_nuLnu_theta, 2);
h5_nuLnu_0 = h5_nuLnu_theta(:,1) + h5_nuLnu_theta(:,5);
h5_nuLnu_1 = h5_nuLnu_theta(:,2) + h5_nuLnu_theta(:,6);
h5_nuLnu_2 = h5_nuLnu_theta(:,3) + h5_nuLnu_theta(:,7);
h5_nuLnu_3 = h5_nuLnu_theta(:,4) + h5_nuLnu_theta(:,8);


lws = 1.5;
h5_total_col = [0.00,0.20,0.80];
h5_col_0 = [0.00,0.62,0.45];
h5_col_1 = [0.84,0.37,0.00];
h5_col_2 = [0.55,0.16,0.70];
h5_col_3 = [0.93,0.69,0.13];


%% PRL-style figure parameters
alp = 0.8;
fig_width_cm  = 17.8*alp;
fig_height_cm = 12.0*alp;

%% Create figure
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0, 0, fig_width_cm, fig_height_cm]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [fig_width_cm, fig_height_cm]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0, 0, fig_width_cm, fig_height_cm]);
set(gcf, 'Renderer', 'painters');
hold on

stairs(nu_left, nu_Lnu_all, ...
    'LineWidth', lws, 'Color', [0,0,0])

stairs(h5_nu, h5_nuLnu_total, '--', ...
    'LineWidth', 2, 'Color', h5_total_col)




hold on 

stairs(nu_left, nu_Lnu_0, 'LineWidth', lws)

stairs(nu_left, nu_Lnu_1, 'LineWidth', lws)

stairs(nu_left, nu_Lnu_2, 'LineWidth', lws)

stairs(nu_left, nu_Lnu_3, 'LineWidth', lws)

stairs(h5_nu, h5_nuLnu_0, ':', 'LineWidth', 1.4, 'Color', h5_col_0)
stairs(h5_nu, h5_nuLnu_1, ':', 'LineWidth', 1.4, 'Color', h5_col_1)
stairs(h5_nu, h5_nuLnu_2, ':', 'LineWidth', 1.4, 'Color', h5_col_2)
stairs(h5_nu, h5_nuLnu_3, ':', 'LineWidth', 1.4, 'Color', h5_col_3)


stairs(nu_left, nu_Lnu_all+eps, ...
    'LineWidth', lws, 'Color', [0,0,0])
% 
% loglog(obs_nu_c,output(:,4))
% 
% 
% loglog(obs_nu_c,output(:,5))

lg = legend('forward', 'spectrum.h5 total', 'n=0', ...
    'n=1', 'n=2', 'n=3', 'h5 n=0', 'h5 n=1', ...
    'h5 n=2', 'h5 n=3', 'Location', 'northwest');

lg.NumColumns = 3;





xlim([1e10,1e23])
ylim([1e30,1e40])

ax = gca;

ax.XScale = "log";
ax.YScale = 'log';

xlabel(ax, '$\nu$', 'Interpreter', 'latex')
ylabel(ax, '$\nu L_\nu$', 'Interpreter', 'latex')




cd(old_dir);
exportgraphics(gcf, fullfile(old_dir, 'SgrA_fdf_with_spectrum_h5.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(old_dir, 'SgrA_fdf_with_spectrum_h5.png'), 'Resolution', 300);
