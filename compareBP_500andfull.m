% Load both files into separate structures
data_full = load("H:\PhD work\bp_cross\NewZealand\result\2025recaculate_wholedepth.mat");
data_500m = load("H:\PhD work\bp_cross\NewZealand\result\2025recaculate_500depth.mat");

% Define target locations as [lat, lon]
target_coords = [
    -39.869, 178.110;
    -39.869, 178.383;
    -40.063, 178.566
];

% Preallocate output arrays for each dataset
nLocs = size(target_coords,1);
nTimes = size(data_full.eff_bpga,2);  % assume same length in both datasets

% Create arrays to hold nearest time series from both depths
bpga_full = zeros(nLocs, nTimes);
bpssh_full = zeros(nLocs, nTimes);
bpt_full = zeros(nLocs, nTimes);

bpga_500 = zeros(nLocs, nTimes);
bpssh_500 = zeros(nLocs, nTimes);
bpt_500 = zeros(nLocs, nTimes);

% Loop through each target coordinate
for k = 1:nLocs
    target_lat = target_coords(k, 1);
    target_lon = target_coords(k, 2);
    
    % === Find nearest index in full-depth data ===
    dist_full = sqrt((data_full.eff_sta(:,1) - target_lon).^2 + ...
                     (data_full.eff_sta(:,2) - target_lat).^2);
    [~, idx_full] = min(dist_full);

    % === Find nearest index in 500m-depth data ===
    dist_500 = sqrt((data_500m.eff_sta(:,1) - target_lon).^2 + ...
                    (data_500m.eff_sta(:,2) - target_lat).^2);
    [~, idx_500] = min(dist_500);
    
    fprintf('Target #%d:\n', k);
    fprintf('  Full-depth nearest = (%.3f, %.3f)\n', ...
            data_full.eff_sta(idx_full,2), data_full.eff_sta(idx_full,1));
    fprintf('  500m-depth nearest = (%.3f, %.3f)\n', ...
            data_500m.eff_sta(idx_500,2), data_500m.eff_sta(idx_500,1));

    % Extract time series
    bpga_full(k,:)  = data_full.eff_bpga(idx_full, :);
    bpssh_full(k,:) = data_full.eff_bpssh(idx_full, :);
    bpt_full(k,:)   = data_full.eff_bpt(idx_full, :);

    bpga_500(k,:)  = data_500m.eff_bpga(idx_500, :);
    bpssh_500(k,:) = data_500m.eff_bpssh(idx_500, :);
    bpt_500(k,:)   = data_500m.eff_bpt(idx_500, :);
    idx_full
    idx_500
end

%% === Example: plot comparison for target #1 ===
% Define exact location from the target_coords array
k =3;
target_lat = target_coords(k, 1);
target_lon = target_coords(k, 2);

figure;

subplot(2,1,1);
plot(data_full.Timedate,bpga_full(k,:), 'b', 'DisplayName','Full-depth','LineWidth',1); hold on;
plot(data_full.Timedate,bpga_500(k,:), 'r', 'DisplayName','0-500m','LineWidth',1);
plot(data_full.Timedate,bpssh_full(k,:), 'k','DisplayName','SSH pressure','LineWidth',1);
title(sprintf('Steric Pressure Comparison at (Lat: %.3f, Lon: %.3f)', target_lat, target_lon));
legend('Location','southwest','FontSize',10); ylabel('Pressure(hPa)');
ax=gca;
ax.FontSize = 16; 

subplot(2,1,2);
plot(data_full.Timedate,bpt_full(k,:), 'b','LineWidth',1,'DisplayName','Full-depth'); hold on;
plot(data_full.Timedate,bpt_500(k,:), 'r','LineWidth',1, 'DisplayName','0-500m');
title(sprintf('Bottom Pressure Comparison at (Lat: %.3f, Lon: %.3f)', target_lat, target_lon));
ylabel('Pressure(hPa)');
xlabel('Time');
ax=gca;
ax.FontSize = 16; 
legend('Location','southwest','FontSize',10);