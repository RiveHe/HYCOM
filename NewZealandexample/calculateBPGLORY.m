ssh = ncread("H:\PhD work\GLORY\temp\mercatorglorys12v1_gl12_mean_201407.nc",'zos');
sal = ncread("H:\PhD work\GLORY\temp\mercatorglorys12v1_gl12_mean_201407.nc",'so');
tem = ncread("H:\PhD work\GLORY\temp\mercatorglorys12v1_gl12_mean_201407.nc",'thetao');
%%
depth = ncread("H:\PhD work\GLORY\temp\mercatorglorys12v1_gl12_mean_201407.nc",'depth');
time = ncread("H:\PhD work\GLORY\temp\mercatorglorys12v1_gl12_mean_201407.nc",'time');
referencetime = datetime('1950-01-01 00:00:00'); 
Timedate = hours(time) + referencetime; 
%%
addpath("H:\PhD work\Miguel\seawater_ver3_3.1")
addpath(genpath("H:\PhD work\natsortfiles"))

% === Step 1: Setup paths and file list ===
data_path = "H:\PhD work\GLORY\";
files = dir(fullfile(data_path, '*.nc'));
files = files(~[files.isdir]);
files = natsortfiles(files);  % Natural sort

% === Step 2: Load depth and coordinate info ===
depth = ncread(fullfile(files(1).folder, files(1).name), 'depth');  % Full column
P = 1040 * 9.8 * depth / 1e4;  % Convert to hPa

lat = ncread(fullfile(files(1).folder, files(1).name), 'latitude');
lon = ncread(fullfile(files(1).folder, files(1).name), 'longitude');
m1 = length(lon);  % Number of longitude points
n1 = length(lat);  % Number of latitude points
origin = datetime(1950,1,1);

% === Constants ===
rho = 1040; g = 9.8;
bp_all = [];
bpga_all = [];
time_all = [];
bpssh_all = [];

% === Step 3: Loop through each file ===
for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    fprintf('Processing %s (%d of %d)\n', files(i).name, i, length(files));

    sal = ncread(filePath, 'so');        % [lon x lat x depth x time]
    tem = ncread(filePath, 'thetao');
    ssh = ncread(filePath, 'zos');
    time = ncread(filePath, 'time');     % [time]
    time_vec = hours(time) + datetime(1950,1,1);

    for t = 1:length(time)
        bp_t = nan(m1, n1);
        bpga_t = nan(m1, n1);
        for x = 1:m1
            for y = 1:n1
                profile_s = squeeze(sal(x,y,:,t));
                profile_t = squeeze(tem(x,y,:,t));
                if all(isnan(profile_s))
                    continue;
                end
                nan_idx = find(isnan(profile_s), 1);
                if isempty(nan_idx)
                    % All valid levels
                    ga = sw_gpan(profile_s, profile_t, P);
                elseif nan_idx > 1
                    % Partially valid
                    ga = sw_gpan(profile_s(1:nan_idx-1), profile_t(1:nan_idx-1), P(1:nan_idx-1));
                else
                    continue;
                end
                bpga_t(x,y) = rho * g * ga(end) / g;
                bp_t(x,y)   = rho * g * (ssh(x,y,t) - ga(end) / g);
                bpssh = rho * g * ssh(x,y,t);
            end
        end
        bp_all(:,:,end+1) = bp_t;
        bpga_all(:,:,end+1) = bpga_t;
        bpssh_all(:,:,end+1) = bpssh;
         
        time_all(end+1) = days(time_vec(t) - origin);  % stores as days since 1950-01-01

    end
end

fprintf("✅ Finished full-column processing. Total time steps: %d\n", length(time_all));
%%
%%
bp2_N = nan(size(bp_all)); 
for x1=1:m1
    for y1=1:n1
        bp2_N(x1,y1,:)=bp_all(x1,y1,1:end)-nanmean(bp_all(x1,y1,1:end));
    end
end


% bpssh2 = nan(size(bpssh_all)); 
% for x1=1:m1
%     for y1=1:n1
%         bpssh2(x1,y1,:)=bpssh_all(x1,y1,1:end)-nanmean(bpssh_all(x1,y1,1:end));
%     end
% end

bpga2 = nan(size(bpga_all)); 
for x1=1:m1
    for y1=1:n1
        bpga2(x1,y1,:)=bpga_all(x1,y1,1:end)-nanmean(bpga_all(x1,y1,1:end));
    end
end
%%
lat = lat';
desiredSize = [4320, 2041]; %lon,lat

% Replicate the columns to create Lon with the desired size
Lon1 = repmat(lon, 1, desiredSize(2));
Lat1 = repmat(lat, desiredSize(1), 1);
%%
yr = 3;
eff_num = sum(sum(~isnan(bp2_N(:,:,yr))));
eff_sta = zeros(eff_num,2);
eff_bpga = zeros(eff_num,size(bpga2,3));
% eff_bpssh = zeros(eff_num,size(bpssh2,3));
eff_bpt = zeros(eff_num,size(bp2_N,3));
num=1;

for is1 = 1:size(Lon1,1)
    for is2 = 1:size(Lon1,2)
        if ~isnan(bp2_N(is1,is2,yr))
            eff_sta(num,:) = [Lon1(is1,is2), Lat1(is1,is2)];
            eff_bpga(num,:) = squeeze(bpga2(is1,is2,:)) / 10^2;
            % eff_bpssh(num,:) = squeeze(bpssh2(is1,is2,:)) / 10^2;
            eff_bpt(num,:) = squeeze(bp2_N(is1,is2,:)) / 10^2;
            num = num + 1;
        end
    end
end
%%
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
% eff_bpssh_nearest = eff_bpssh(nearest_idx, :);
eff_bpt_nearest = eff_bpt(nearest_idx, :);
%%
referencetime = datetime('1950-01-01 00:00:00'); 
Timedate = referencetime + hours(time_all);
%%
figure;
plot( eff_bpt_nearest / 100, 'r', 'LineWidth', 1);
legend('HYCOM','prs\_bsp4')
xlabel('Time');
ylabel('Effective Bottom Pressure (dbar)');
title(['Bottom Pressure at Lat: ', num2str(nearest_lat), ', Lon: ', num2str(nearest_lon-360)]);
grid on;