addpath("H:\PhD work\Miguel\seawater_ver3_3.1")
addpath(genpath("H:\PhD work\natsortfiles")) % If you use natsortfiles (optional)

% === Step 1: Get file lists (excluding subfolders like 'data2') ===
sal_path = "H:\PhD work\bp_cross\NewZealand\data\salinity\";
temp_path = "H:\PhD work\bp_cross\NewZealand\data\temperature\";

C = dir(fullfile(sal_path, '*.nc'));
D = dir(fullfile(temp_path, '*.nc'));

% Remove files in subdirectories
C = C(~[C.isdir]);
D = D(~[D.isdir]);

% Natural sort by filenames (ensures 1.nc, 2.nc, ..., 12.nc)
C = natsortfiles(C);
D = natsortfiles(D);

% === Step 2: Load coordinates ===
Lat1 = ncread(fullfile(C(3).folder, C(3).name), 'lat');
Lon1 = ncread(fullfile(C(3).folder, C(3).name), 'lon');
Depth1 = ncread(fullfile(C(3).folder, C(3).name), 'depth');
ssh1 = ncread("H:\PhD work\bp_cross\NewZealand\data\ssh\1.nc", 'surf_el');

m1 = length(Lon1);
n1 = length(Lat1);
bp1_N = nan(size(ssh1));
bpga = nan(size(ssh1));

rho = 1040; % (kg/m3)
g = 9.8;    % (m/s2)
P1 = rho * Depth1 * g / 1e4; % pressure in hPa

l1_N = zeros(length(C),1);
l2_N = zeros(length(D),1);
nz = length(Depth1);  % total vertical layers

% === Preallocate 4D arrays ===
bpga = nan(m1, n1, nz, sum(l1_N));
bp1_N = nan(m1, n1, nz, sum(l1_N));  % bpt

% Reset counter
num1 = 1;

for mo1 = 1:length(C)
    sal_file = fullfile(C(mo1).folder, C(mo1).name);
    temp_file = fullfile(D(mo1).folder, D(mo1).name);

    fprintf('Processing: %s + %s\n', C(mo1).name, D(mo1).name);

    S1 = ncread(sal_file, 'salinity');
    T1 = ncread(temp_file, 'water_temp');

    l1_N(mo1) = size(S1,4);
    l2_N(mo1) = size(T1,4);

    for yr1 = 1:size(S1,4)
        for x1 = 1:m1
            for y1 = 1:n1
                for z1 = 2:nz
                    % Check for valid salinity profile up to current depth
                    if ~isnan(S1(x1,y1,z1,yr1)) && all(~isnan(S1(x1,y1,1:z1,yr1)))
                        sal_prof = squeeze(S1(x1,y1,1:z1,yr1));
                        temp_prof = squeeze(T1(x1,y1,1:z1,yr1));
                        ga = sw_gpan(sal_prof, temp_prof, P1(1:z1));  % geopotential anomaly profile

                        % Store bpga and bpt at all depths 1:z1
                        bpga(x1, y1, 1:z1, num1) = rho * ga;
                        bp1_N(x1, y1, 1:z1, num1) = rho * g * ssh1(x1, y1, num1) - bpga(x1, y1, 1:z1, num1);
                    end
                end
            end
        end
        num1 = num1 + 1;
    end

    fprintf('Completed file %2d (%s): %d time steps added, num1 now = %d\n', ...
            mo1, C(mo1).name, l1_N(mo1), num1 - 1);
end
%%

bp2_N = nan(size(bp1_N)); 
for x1=1:m1
    for y1=1:n1
        for z1 = 1:nz
        bp2_N(x1,y1,z1,:)=bp1_N(x1,y1,z1,1:end)-nanmean(bp1_N(x1,y1,z1,1:end));
        end
    end
end

bpga2 = nan(size(bpga)); 
for x1=1:m1
    for y1=1:n1
        for z1 = 1:nz
        bpga2(x1,y1,z1,:)=bpga(x1,y1,z1,1:end)-nanmean(bpga(x1,y1,z1,1:end));
        end
    end
end

bpssh = rho*g*ssh1; 
bpssh2 = nan(size(bpssh)); 
for x1=1:m1
    for y1=1:n1
        bpssh2(x1,y1,:)=bpssh(x1,y1,1:end)-nanmean(bpssh(x1,y1,1:end));
    end
end
%%
% Define target coordinates [lat, lon]
target_coords = [
    -39.869, 178.110;
    %-39.869, 178.383;
    %-40.063, 178.566
];

% Ensure Lat1 and Lon1 are 2D grid arrays (as from ROMS/HYCOM)
% If they are vectors, convert to meshgrid
if isvector(Lat1) && isvector(Lon1)
    [Lon1_grid, Lat1_grid] = meshgrid(Lon1, Lat1);
else
    Lon1_grid = Lon1;
    Lat1_grid = Lat1;
end

% Initialize output
nearest_indices = zeros(size(target_coords,1), 2);  % [ilat, ilon]

for i = 1:size(target_coords, 1)
    target_lat = target_coords(i, 1);
    target_lon = target_coords(i, 2);

    % Compute distance matrix (Euclidean distance in lat-lon space)
    dist = sqrt((Lat1_grid - target_lat).^2 + (Lon1_grid - target_lon).^2);

    % Find index of the minimum distance
    [~, min_idx] = min(dist(:));
    [ilat, ilon] = ind2sub(size(dist), min_idx);

    % Save indices
    nearest_indices(i,:) = [ilat, ilon];

    fprintf('Target #%d: Nearest grid point at (%.3f, %.3f), indices = (%d, %d)\n', ...
        i, Lat1_grid(ilat, ilon), Lon1_grid(ilat, ilon), ilat, ilon);
end
%%
% Define a specific location (adjust indices as needed)

% Extract profiles at this location: [depth x time]
profile = squeeze(bpga2(ilon, ilat, :, :)) / 100;   % steric pressure (hPa)
bpt_profile = squeeze(bp2_N(ilon, ilat, :, :)) / 100; % total pressure anomaly (hPa)

%% mask out the 0 data
valid_depth_idx = any(bpt_profile ~= 0, 2);  % size = [depth x 1]
bpt_profile_masked = bpt_profile(valid_depth_idx, :);
profile_masked     = profile(valid_depth_idx, :);
depth_masked       = Depth1(valid_depth_idx);
fprintf('Depth range: %.1f m to %.1f m\n', min(depth_masked), max(depth_masked));
%%
figure;

% Panel 1
subplot(3,1,1)
imagesc(Timedate, depth_masked, profile_masked);
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth (m)');
title(sprintf('Cumulative Steric Pressure (hPa) profile at (Lat: %.3f, Lon: %.3f)',target_lat,target_lon));
%colorbar;
ax=gca;
ax.FontSize = 14; 

% Panel 2
subplot(3,1,2)
imagesc(Timedate, depth_masked, bpt_profile_masked);
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth (m)');
title(sprintf('Bottom Pressure (hPa) Profile at (Lat: %.3f, Lon: %.3f)',target_lat,target_lon));
%colorbar;
ax=gca;
ax.FontSize = 14; 

% Panel 3
subplot(3,1,3)
plot(Timedate, squeeze(bpt_profile_masked(end,:)), 'b','DisplayName','Full-depth');
hold on
plot(Timedate, squeeze(bpt_profile_masked(34,:)), 'r', 'DisplayName','0–1250 m');
xlabel('Time'); ylabel('Pressure (hPa)');
title(sprintf('Bottom Pressure Time Series at (Lat: %.3f, Lon: %.3f)', target_lat, target_lon));

ax=gca;
ax.FontSize = 14; 
legend('Location', 'northeast','FontSize',11);

%% plot std for each layer's differences
nz = length(depth_masked);
nt = size(profile_masked, 2);

% Preallocate
std_pressure_diff = nan(nz - 1, 1);

% Loop over layer pairs
for i = 1:(nz - 1)
    % Pressure difference between adjacent depth layers
    layer_pressure_diff = profile_masked(i + 1, :) - profile_masked(i, :);  % [1 × nt]

    % Standard deviation over time
    std_pressure_diff(i) = std(layer_pressure_diff, 'omitnan');
    
end

% Find index of max std
[max_std, idx_max] = max(std_pressure_diff);

% Depths involved
depth_top = depth_masked(idx_max);       % top of the layer
depth_bottom = depth_masked(idx_max + 1);% bottom of the layer
depth_center = (depth_top + depth_bottom) / 2;

% === Plot ===
figure;
plot(std_pressure_diff, depth_masked(2:end), '-o', 'LineWidth', 2);
set(gca, 'YDir', 'reverse');
xlabel('Std Dev of Layer Steric Pressure Difference (hPa)', 'FontSize', 14);
ylabel('Depth (m)', 'FontSize', 14);
title(sprintf('\\sigma(\\Delta P) over Time between Adjacent Layers at (Lat: %.3f, Lon: %.3f)', ...
    target_lat, target_lon), 'FontSize', 16);
grid on;
ax = gca;
ax.FontSize = 13;

% Mark max std point
hold on;
plot(max_std, depth_bottom, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% Add annotation
text(max_std, depth_bottom, ...
    sprintf('  Max σ = %.2f hPa\n  Layer: %.1f–%.1f m', max_std, depth_top, depth_bottom), ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'red');

%% plot normalized difference: delta P/ delta Z

%% plot std for each layer's differences (normalized by layer thickness)
nz = length(depth_masked);
nt = size(profile_masked, 2);

% Preallocate
std_pressure_diff = nan(nz - 1, 1);
layer_thickness = nan(nz - 1, 1);

% Loop over layer pairs
for i = 1:(nz - 1)
    % Compute layer thickness between adjacent levels
    layer_thickness(i) = depth_masked(i + 1) - depth_masked(i);

    % Pressure difference between adjacent depth layers
    layer_pressure_diff = profile_masked(i + 1, :) - profile_masked(i, :);  % [1 × nt]

    % Normalize by layer thickness (in meters)
    normalized_diff = layer_pressure_diff / layer_thickness(i);  % hPa/m

    % Standard deviation over time
    std_pressure_diff(i) = std(normalized_diff, 'omitnan');
end

% === Plot ===
depth_mid = (depth_masked(1:end-1) + depth_masked(2:end)) / 2;


%% Compute normalized steric pressure differences
nz = length(depth_masked);       % number of vertical levels
nt = size(profile_masked, 2);    % number of time steps

% Preallocate output matrix: (nz-1) x nt
normalized_diff = nan(nz - 1, nt);
depth_mid = (depth_masked(1:end-1) + depth_masked(2:end)) / 2;

for i = 1:(nz - 1)
    % Calculate pressure difference between adjacent layers
    layer_pressure_diff = profile_masked(i + 1, :) - profile_masked(i, :);  % [1 x nt]

    % Layer thickness
    delta_z = depth_masked(i + 1) - depth_masked(i);

    % Normalize by thickness
    normalized_diff(i, :) = layer_pressure_diff / delta_z;  % [1 x nt], units: hPa/m
end

%% Plot
figure;
subplot(211)
imagesc(Timedate, depth_mid, normalized_diff);
set(gca, 'YDir', 'reverse');  % So that depth increases downward
xlabel('Time');
ylabel('Depth (m)');
title(sprintf('Normalized Steric Pressure Gradient (\\DeltaP/\\Deltaz) at (Lat: %.3f, Lon: %.3f)', ...
    target_lat, target_lon));
colorbar;
cb = colorbar;
cb.Label.String = 'hPa/m';
ax = gca;
ax.FontSize = 13;

subplot(212)
plot(std_pressure_diff, depth_mid, '-o', 'LineWidth', 2);
set(gca, 'YDir', 'reverse');
xlabel('Std Dev of Normalized Steric Pressure Difference (hPa/m)', 'FontSize', 14);
ylabel('Midpoint Depth (m)', 'FontSize', 14);
title(sprintf('\\sigma(\\Delta P/\\Delta z) over Time at (Lat: %.3f, Lon: %.3f)', ...
    target_lat, target_lon), 'FontSize', 16);
grid on;
ax = gca;
ax.FontSize = 13;
