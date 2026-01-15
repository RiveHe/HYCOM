addpath("H:\PhD work\Miguel\seawater_ver3_3.1")
addpath(genpath("H:\PhD work\natsortfiles"))

% === Step 1: Get file lists ===
sal_path = "H:\PhD work\bp_cross\NewZealand\data\salinity\";
temp_path = "H:\PhD work\bp_cross\NewZealand\data\temperature\";

C = dir(fullfile(sal_path, '*.nc'));
D = dir(fullfile(temp_path, '*.nc'));
C = C(~[C.isdir]);
D = D(~[D.isdir]);

C = natsortfiles(C);
D = natsortfiles(D);

% === Step 2: Load coordinates and depth (trim to ≤ 500 m) ===
Lat1 = ncread(fullfile(C(3).folder, C(3).name), 'lat');
Lon1 = ncread(fullfile(C(3).folder, C(3).name), 'lon');
Depth1_full = ncread(fullfile(C(3).folder, C(3).name), 'depth');

depthMask = Depth1_full <= 500;
Depth1 = Depth1_full(depthMask);
P1 = 1040 * 9.8 * Depth1 / 1e4;  % hPa

ssh1 = ncread("H:\PhD work\bp_cross\NewZealand\data\ssh\1.nc", 'surf_el');
m1 = length(Lon1);
n1 = length(Lat1);
bp1_N = nan(size(ssh1));
bpga = nan(size(ssh1));


rho = 1040;
g = 9.8;

num1 = 1;
l1_N = zeros(length(C),1);
l2_N = zeros(length(D),1);

% === Step 3: Loop and compute using top 0–500 m ===
for mo1 = 1:length(C)
    sal_file = fullfile(C(mo1).folder, C(mo1).name);
    temp_file = fullfile(D(mo1).folder, D(mo1).name);

    fprintf('Processing: %s + %s\n', C(mo1).name, D(mo1).name);

    S_full = ncread(sal_file, 'salinity');
    T_full = ncread(temp_file, 'water_temp');

    S1 = S_full(:,:,depthMask,:);
    T1 = T_full(:,:,depthMask,:);

    l1_N(mo1) = size(S1,4);
    l2_N(mo1) = size(T1,4);

    for yr1 = 1:size(S1,4)
        for x1 = 1:m1
            for y1 = 1:n1
                if ~isnan(S1(x1,y1,end,yr1))
                    ga = sw_gpan(squeeze(S1(x1,y1,:,yr1)), squeeze(T1(x1,y1,:,yr1)), P1);
                    bp1_N(x1,y1,num1) = rho*g*(ssh1(x1,y1,num1) - ga(end)/g);
                    bpga(x1,y1,num1)   = rho*g*ga(end)/g;
                else
                    for depth1 = 2:size(S1,3)
                        if isnan(S1(x1,y1,depth1,yr1)) && ~isnan(S1(x1,y1,depth1-1,yr1))
                            ga = sw_gpan(squeeze(S1(x1,y1,1:(depth1-1),yr1)), ...
                                         squeeze(T1(x1,y1,1:(depth1-1),yr1)), P1(1:(depth1-1)));
                            bp1_N(x1,y1,num1) = rho*g*(ssh1(x1,y1,num1) - ga(end)/g);
                            bpga(x1,y1,num1)  = rho*g*ga(end)/g;
                          
                        end
                    end
                end
            end
        end
        num1 = num1 + 1;
    end

    fprintf('Completed file %2d (%s): %d time steps added, num1 now = %d\n', ...
            mo1, C(mo1).name, l1_N(mo1), num1 - 1);
end

% === Step 4: Total time steps ===
fprintf('Total salinity time steps: %d\n', sum(l1_N));
fprintf('Total temperature time steps: %d\n', sum(l2_N));

% === Step 5: Check time order ===
expected_timesteps = sum(l1_N);
actual_timesteps = size(bp1_N, 3);
if expected_timesteps == actual_timesteps
    fprintf('✅ Time steps match: %d time steps accumulated in order.\n', actual_timesteps);
else
    error('❌ Mismatch in expected (%d) and actual (%d) time steps.', expected_timesteps, actual_timesteps);
end

% === Step 6: Print file-by-file breakdown ===
fprintf('\nTime steps per file (should be in order):\n');
cumulative = 0;
for i = 1:length(C)
    cumulative = cumulative + l1_N(i);
    fprintf('File %2d (%s): %d time steps, cumulative end index = %d\n', ...
            i, C(i).name, l1_N(i), cumulative);
end

return


%%
bp2_N = nan(size(bp1_N)); 
for x1=1:m1
    for y1=1:n1
        bp2_N(x1,y1,:)=bp1_N(x1,y1,1:end)-nanmean(bp1_N(x1,y1,1:end));
    end
end

bpssh = rho*g*ssh1; 
bpssh2 = nan(size(bpssh)); 
for x1=1:m1
    for y1=1:n1
        bpssh2(x1,y1,:)=bpssh(x1,y1,1:end)-nanmean(bpssh(x1,y1,1:end));
    end
end

bpga2 = nan(size(bpga)); 
for x1=1:m1
    for y1=1:n1
        bpga2(x1,y1,:)=bpga(x1,y1,1:end)-nanmean(bpga(x1,y1,1:end));
    end
end
%%
Lat1 = Lat1';
desiredSize = [38, 76]; %lon,lat

% Replicate the columns to create Lon with the desired size
Lon1 = repmat(Lon1, 1, desiredSize(2));
Lat1 = repmat(Lat1, desiredSize(1), 1);
%%
yr = 100;
eff_num = sum(sum(~isnan(bp2_N(:,:,yr))));
eff_sta = zeros(eff_num,2);
eff_bpga = zeros(eff_num,size(bpga2,3));
eff_bpssh = zeros(eff_num,size(bpssh2,3));
eff_bpt = zeros(eff_num,size(bp2_N,3));
num=1;

for is1 = 1:size(Lon1,1)
    for is2 = 1:size(Lon1,2)
        if ~isnan(bp2_N(is1,is2,yr))
            eff_sta(num,:) = [Lon1(is1,is2), Lat1(is1,is2)];
            eff_bpga(num,:) = squeeze(bpga2(is1,is2,:)) / 10^2;
            eff_bpssh(num,:) = squeeze(bpssh2(is1,is2,:)) / 10^2;
            eff_bpt(num,:) = squeeze(bp2_N(is1,is2,:)) / 10^2;
            num = num + 1;
        end
    end
end
%%
eff_bpga_nearest = [];
eff_bpssh_nearest = [];
eff_bpt_nearest = [];
% Target location
target_lat = -39.869;  % Latitude of the target location
target_lon = 178.110; % Longitude of the target location

% eff_sta contains the coordinates (longitude, latitude) of the effective stations

% Calculate the distance between each point in eff_sta and the target location
distances = sqrt((eff_sta(:,1) - target_lon).^2 + (eff_sta(:,2) - target_lat).^2);

% Find the index of the minimum distance
[~, nearest_idx] = min(distances);

% Nearest point's coordinates
nearest_lat = eff_sta(nearest_idx, 2);
nearest_lon = eff_sta(nearest_idx, 1);

% Display the nearest point's coordinates
disp(['Nearest Latitude: ', num2str(nearest_lat)]);
disp(['Nearest Longitude: ', num2str(nearest_lon)]);

% Extract the corresponding time series data using nearest_idx
eff_bpga_nearest = eff_bpga(nearest_idx, :);
eff_bpssh_nearest = eff_bpssh(nearest_idx, :);
eff_bpt_nearest = eff_bpt(nearest_idx, :);

%%
referenceDate = datetime('2000-01-01 00:00:00'); % Reference date for time conversion
timeData = ncread("H:\PhD work\bp_cross\NewZealand\data\ssh\1.nc", 'time');
Timedate = referenceDate + hours(timeData);
%%
figure;
plot(Timedate, eff_bpt_nearest / 100, 'r', 'LineWidth', 1);
legend('HYCOM','prs\_bsp4')
xlabel('Time');
ylabel('Effective Bottom Pressure (dbar)');
title(['Bottom Pressure at Lat: ', num2str(nearest_lat), ', Lon: ', num2str(nearest_lon-360)]);
grid on;
%%
function S = sort_nat(files)
    [~, idx] = sort_nat_internal(files);
    S = files(idx);
end

function [sorted, idx] = sort_nat_internal(files)
    [~, idx] = sort(lower(files));
    sorted = files(idx);
end
