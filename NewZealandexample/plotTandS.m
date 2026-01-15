T_feb = ncread("H:\PhD work\bp_cross\NewZealand\data\temperature\8.nc",'water_temp');
S_feb = ncread("H:\PhD work\bp_cross\NewZealand\data\salinity\8.nc",'salinity');
time_t =  ncread("H:\PhD work\bp_cross\NewZealand\data\temperature\8.nc",'time');
time_T = hours(time_t) + referenceDate;
T_profile = squeeze(T_feb(ilon,ilat,:,:)); 
T_masked = T_profile(valid_depth_idx, :);
S_profile = squeeze(S_feb(ilon,ilat,:,:)); 
S_masked = S_profile(valid_depth_idx, :);
%%
figure
subplot(312)
imagesc(time_T, depth_masked,T_masked);
colorbar
title(sprintf('Temperature profile at at (Lat: %.3f, Lon: %.3f) in 2015 Feb', target_lat, target_lon));
xlabel('Time'); ylabel('Depth (m)');
%xlim([datetime(2015,2,1), datetime(2015,2,28)]);

subplot(313)
imagesc(time_T, depth_masked,S_masked);
colorbar
title(sprintf('Salinity profile at (Lat: %.3f, Lon: %.3f) in 2015 Feb', target_lat, target_lon));
xlabel('Time'); ylabel('Depth (m)');
%xlim([datetime(2015,2,1), datetime(2015,2,28)]);

subplot(311)
imagesc(Timedate, depth_masked, profile_masked);
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth (m)');
title(sprintf('Cumulative Steric Pressure (hPa) profile at (Lat: %.3f, Lon: %.3f) in 2015 Feb', target_lat, target_lon));
xlim([datetime(2015,2,1), datetime(2015,2,28)]);
%title('Steric Pressure in Feb 2015')
colorbar
caxis([-10,10])
%% read all the T and S files 
% Initialize empty arrays for time concatenation
T_all = [];
S_all = [];

% Loop through files 1.nc to 12.nc
for i = 1:12
    % Construct file paths
    T_file = fullfile("H:\PhD work\bp_cross\NewZealand\data\temperature", sprintf('%d.nc', i));
    S_file = fullfile("H:\PhD work\bp_cross\NewZealand\data\salinity", sprintf('%d.nc', i));

    % Read temperature and salinity
    T = ncread(T_file, 'water_temp');  % [lon, lat, depth, time]
    S = ncread(S_file, 'salinity');    % [lon, lat, depth, time]

    % Concatenate along 4th (time) dimension
    T_all = cat(4, T_all, T);
    S_all = cat(4, S_all, S);
end

%% Extract profile at (ilon, ilat)
T_profile = squeeze(T_all(ilon, ilat, :, :));  % [depth, time]
S_profile = squeeze(S_all(ilon, ilat, :, :));  % [depth, time]

% Find valid (non-zero or non-NaN) depths — pick your masking criteria
valid_depth_idx = all(~isnan(T_profile), 2) & all(T_profile ~= 0, 2);

% Apply depth mask
T_masked = T_profile(valid_depth_idx, :);
S_masked = S_profile(valid_depth_idx, :);

%% Plot to verify
figure;
subplot(3,1,1)
imagesc(Timedate, depth_masked,T_masked);
title(sprintf('Temperature profile at at (Lat: %.3f, Lon: %.3f)', target_lat, target_lon));
colorbar;
colorbar;
ax=gca;
ax.FontSize = 14; 

subplot(3,1,2)
imagesc(Timedate, depth_masked,S_masked);
title(sprintf('Salinity profile at (Lat: %.3f, Lon: %.3f)', target_lat, target_lon));
colorbar;
ax=gca;
ax.FontSize = 14; 

subplot(3,1,3)
imagesc(Timedate, depth_masked, profile_masked);
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth (m)');
title(sprintf('Cumulative Steric Pressure (hPa) profile at (Lat: %.3f, Lon: %.3f)',target_lat,target_lon));
colorbar;
ax=gca;
ax.FontSize = 14; 