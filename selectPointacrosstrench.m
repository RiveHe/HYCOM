clear all;

% Read bathymetry data from netCDF file
filename_nc = '/Volumes/My Data/PhD work/bp_cross/Chile/dataset/ncmap/chile2.nc';
gebconc = netcdf.open(filename_nc, 'NOWRITE');

lat = double(netcdf.getVar(gebconc,0));
lon = double(netcdf.getVar(gebconc,1));
dep = double(netcdf.getVar(gebconc,2));
dep = dep';
netcdf.close(gebconc);

% Define a series of latitudinal lines across which you want to find points.
lat_lines = min(lat):1:max(lat);  % You can adjust the spacing based on your needs

% Interpolate depth values along these lines
interpolated_depths = zeros(length(lat_lines), length(lon));
for idx = 1:length(lat_lines)
    interpolated_depths(idx, :) = interp2(lon, lat, dep, lon, repmat(lat_lines(idx), length(lon), 1));
end

% For each latitudinal line, find where the interpolated depth is within the threshold of the desired isobath
depth_range = [0, -5000];
selected_lons = [];
selected_lats = [];

for idx = 1:length(lat_lines)
    % Find indices where interpolated depth is within the desired range
    valid_depth_indices = find(interpolated_depths(idx, :) <= depth_range(1) & interpolated_depths(idx, :) >= depth_range(2));
    
    selected_lats = [selected_lats; repmat(lat_lines(idx), length(valid_depth_indices), 1)];
    selected_lons = [selected_lons; lon(valid_depth_indices)];
end

downsampling_rate = 10; % For instance, take every 5th point. You can adjust this rate based on your needs.

% New containers for downsampled data
downsampled_lons = [];
downsampled_lats = [];

for idx = 1:length(lat_lines)
    % Find indices where interpolated depth is within the desired range
    valid_depth_indices = find(interpolated_depths(idx, :) <= depth_range(1) & interpolated_depths(idx, :) >= depth_range(2));

    % Downsampling
    valid_depth_indices = valid_depth_indices(1:downsampling_rate:end); % Take every nth point

    downsampled_lats = [downsampled_lats; repmat(lat_lines(idx), length(valid_depth_indices), 1)];
    downsampled_lons = [downsampled_lons; lon(valid_depth_indices)];
end

mask = downsampled_lons > -73;
selected_lats = downsampled_lats(mask);
selected_lons = downsampled_lons(mask);

% For latitudes < -31.5, modify the longitude mask
mask_lon2 = -33.5 < selected_lats & selected_lats < -31.5 & selected_lons <= -72.7;

% Remove the points based on the second mask
selected_lats(mask_lon2) = [];
selected_lons(mask_lon2) = [];

mask_lon3 = selected_lats > -31.5 & selected_lons <= -72.5;

selected_lats(mask_lon3) = [];
selected_lons(mask_lon3) = [];

% Get indices of selected_lats values that are less than -32.5 and greater than -33.5
indices = (selected_lats < -33.5)&(selected_lats > -34.5);

% Filter the data using these indices
filtered_lons = selected_lons(indices);
filtered_lats = selected_lats(indices);

% Check if any values match the condition
if ~isempty(filtered_lats)
    filepath = 'downsampled_points_34.xlsx';
    
    % Delete existing file if it exists
    if isfile(filepath)
        delete(filepath);
    end

    % Create a new table with filtered data and proper column names
    T_filtered = table(filtered_lons, filtered_lats, 'VariableNames', {'Longitude', 'Latitude'});

    % Write the filtered data to the Excel file
    writetable(T_filtered, filepath);
end


% Plotting
figure;
contour_levels = [-500, -1000:-1000:-6000];
[C, h] = contour(lon, lat, dep, contour_levels);
clabel(C, h, contour_levels, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 1000);  
colormap(flipud(gray));
axis equal;
axis([-73.5 -71 -34 -28]);

hold on;
scatter(selected_lons, selected_lats, 'MarkerFaceColor', rand(1,3));

% Adding title and labels
title('Points across Isobaths');
xlabel('Longitude');
ylabel('Latitude');
set(gca, 'YDir', 'normal');
grid on;

hold off;
