%clear all

% Load the data
load('/Volumes/My Data/PhD work/bp_cross/Chile/result/2014.mat', 'bp2');
load('/Volumes/My Data/PhD work/bp_cross/Chile/result/timedate.mat', 'Timedate');
Timedate = Timedate;


% Load lat and lon data
lat = ncread('/Volumes/My Data/PhD work/bp_cross/Chile/data/salinity/1.nc', 'lat');
lon = ncread('/Volumes/My Data/PhD work/bp_cross/Chile/data/salinity/1.nc', 'lon');
dep = ncread('/Volumes/My Data/PhD work/Miguel/data2/salinity/2/1.nc', 'depth');
bptdata = lp_demean_bpt3;
%%
% Initialize a cell array to hold the data

% List of file names
%filenames = {'downsampled_points_500m.xlsx','downsampled_points_1000m.xlsx','downsampled_points_2000m.xlsx',...
    %'downsampled_points_3000m.xlsx','downsampled_points_4000m.xlsx','downsampled_points_5000m.xlsx','downsampled_points_6000m.xlsx'};
%filenames = {'alongtrench-500m.xlsx', 'alongtrench-1000m.xlsx','alongtrench-2000m.xlsx','alongtrench-3000m.xlsx','alongtrench-4000m.xlsx','alongtrench-5000m.xlsx','alongtrench-6000m.xlsx'};
%names = {'Alongtrench-500m', 'Alongtrench-1000m', 'Alongtrench-2000m', 'Alongtrench-3000m', 'Alongtrench-4000m','Alongtrench-5000m','Alongtrench-6000m'};
filenames = {'acrosstrench-28.xlsx','acrosstrench-30.xlsx','acrosstrench-32.xlsx'}
names ={'Acrosstrench-28','Acrosstrench-30','Acrosstrench-32'}
%filenames = {'alongtrench-500m.xlsx', 'alongtrench-1000m.xlsx'};
%names = {'Alongtrench-500m', 'Alongtrench-1000m'};
%filenames = {'alongtrench-2000m.xlsx','alongtrench-2000m-1.xlsx'};
%names = {'Alongtrench-2000m-1', 'Alongtrench-2000m-2',};
%filenames = {'alongtrench-3000m.xlsx','alongtrench-3000m-1.xlsx','alongtrench-3000m-2.xlsx','alongtrench-3000m-3.xlsx'};
%names = {'Alongtrench-3000m-1', 'Alongtrench-3000m-2','Alongtrench-3000m-3','Alongtrench-3000m-4'};
%filenames = {'alongtrench-4000m.xlsx','alongtrench-4000m-1.xlsx','alongtrench-4000m-2.xlsx'};
%names = {'Alongtrench-4000m-1', 'Alongtrench-4000m-2','Alongtrench-4000m-3'};
%filenames = {'alongtrench-5000m.xlsx','alongtrench-5000m-1.xlsx','alongtrench-5000m-2.xlsx'};
%names = {'Alongtrench-5000m-1', 'Alongtrench-5000m-2','Alongtrench-5000m-3'};
%filenames = {'alongtrench-6000m.xlsx','alongtrench-6000m-1.xlsx','alongtrench-6000m-2.xlsx','alongtrench-6000m-3.xlsx'};
%names = {'Alongtrench-6000m-1', 'Alongtrench-6000m-2','Alongtrench-6000m-3','Alongtrench-6000m-4'};
num = length(filenames);
allData = cell(1, num);
correlations = cell(num, 1);
distances = cell(num, 1);

% Loop through each file and load the data


% Initialize cell arrays to store latitudes and longitudes
allLatitudes = cell(1, length(filenames));
allLongitudes = cell(1, length(filenames));
% Loop through the cell array to extract latitudes and longitudes

%%
for k = 1:num
    [data, ~] = xlsread(filenames{k});
    allData{k} = data;
    currentData = allData{k};  % Retrieve data from the i-th file
    allLatitudes{k} = currentData(:, 2);  % latitude is the second column
    allLongitudes{k} = currentData(:, 1);
    % Create subplot for the k-th plot

    % Extract the latitudes and longitudes for the k-th file
    latitudes = allLatitudes{k};
    longitudes = allLongitudes{k};

    time_series_data_all = zeros(length(latitudes), length(Timedate));

    % Initialize a 2D array to store the time series data for all locations
    time_series_data_all = zeros(length(latitudes), size(bptdata, 2));
    closest_points = cell(length(latitudes), 3);  % [Lat, Lon, Index]

    for loc_idx = 1:length(latitudes)
        target_lat = latitudes(loc_idx);
        target_lon = longitudes(loc_idx);

        % Find the index in eff_sta that matches the target latitude and longitude
        % Use a threshold for matching, as exact matches are unlikely
        threshold = 0.1; % Adjust this threshold as needed
        distance = sqrt((eff_sta(:,1) - target_lon).^2 + (eff_sta(:,2) - target_lat).^2);
        [~, closest_index] = min(distance);

        if distance(closest_index) <= threshold
            % Extract the time series data for this location from eff_bpt
            time_series_data_all(loc_idx, :) = bptdata(closest_index, :);
        else
            % Handle case where no close match is found
            disp(['No close match found for location ' num2str(loc_idx)]);
        end

        data(loc_idx, 3) = eff_sta(closest_index, 2);%lat
        data(loc_idx, 4) = eff_sta(closest_index, 1);%lon
    end

    %% Demeaned data
    % Initialize arrays to store the correlations and distances
    num_pairs = nchoosek(size(data, 1), 2);  % Number of unique pairs
    temp_correlations = zeros(num_pairs, 1);
    temp_distances = zeros(num_pairs, 1);

    % Radius of the Earth in kilometers
    R = 6371;

    % Counter for storing correlations and distances
    count = 1;

    % Loop through all pairs of points along the trench
    for i = 1:size(data, 1)
        for j = i+1:size(data, 1)
            fprintf('i: %d, j: %d\n', i, j);

            % Convert coordinates from degrees to radians
            lat1 = deg2rad(data(i,3));
            lon1 = deg2rad(data(i, 4));
            lat2 = deg2rad(data(j,3));
            lon2 = deg2rad(data(j, 4));

            % Calculate differences
            dlat = lat2 - lat1;
            dlon = lon2 - lon1;

            % Haversine formula
            a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            distance = R * c;

            % Perform cross-correlation between points i and j
            data1 = time_series_data_all(i, :);
            data2 = time_series_data_all(j, :);

            % Check if one or both of the time series are constant
            if std(data1) == 0 || std(data2) == 0
                % Skip this pair of locations
                continue;
            end

            % Calculate correlation coefficient
            coeff_matrix = corrcoef(data1, data2);
            correlation = coeff_matrix(1, 2);  % Extracting the correlation coefficient

            % Store the correlation and distance
            temp_correlations(count) = correlation;
            temp_distances(count) = distance;

            % Increment the counter
            count = count + 1;
        end
    end

    % Trim the pre-allocated arrays if needed
    temp_correlations = temp_correlations(1:count-1);
    temp_distances = temp_distances(1:count-1);

    % Store the temp_correlations and temp_distances arrays into the cell arrays
    correlations{k} = temp_correlations;
    distances{k} = temp_distances;


    %%
    % Create a new figure
    figure;

    % Subplot for the original time series data
    subplot(3, 1, 1);  % 3 rows, 1 column, 1st subplot
    % Loop through each latitude index to plot its corresponding time series
    for loc_idx = 1:length(latitudes)
        plot(Timedate, time_series_data_all(loc_idx, :));
        hold on;  % Keep the plot, so new lines are added to the same figure
    end
    plot_title_ = ['Time Series Data -' num2str(names{k})];
    ax = gca;  % Get handle to current axes
    ax.FontSize = 14;  % Set tick label font size to 14
    xlabel('Time',FontSize=15);
    ylabel('Pressure (hPa)',FontSize=15);
    title(plot_title_,FontSize=20);
    hold off;

    % Subplot for the demeaned time series data
    subplot(3, 1, 2);  % 3 rows, 1 column, 2nd subplot
    % Loop through each latitude index to plot its corresponding time series
    for loc_idx = 1:length(latitudes)
        plot(Timedate, time_series_data_all_demeaned(loc_idx, :));
        hold on;  % Keep the plot, so new lines are added to the same figure
    end
    ax = gca;  % Get handle to current axes
    ax.FontSize = 14;  % Set tick label font size to 14
    xlabel('Time',FontSize=15);
    ylabel('Demeaned Pressure (hPa)',FontSize=15);
    plot_title_ = ['Demeaned Time Series Data -' num2str(names{k})];
    title(plot_title_,FontSize=20);
    hold off;

    % Subplot for the spatial mean time series
    subplot(3, 1, 3);  % 3 rows, 1 column, 3rd subplot
    plot(Timedate, spatial_mean_time_series);
    ax = gca;  % Get handle to current axes
    ax.FontSize = 14;  % Set tick label font size to 14
    xlabel('Time',FontSize=15);
    ylabel('Spatial Mean Pressure (hPa)',FontSize=15);
    plot_title_ = ['Spatial Mean Time Series Data-' num2str(names{k})];
    title(plot_title_,FontSize=20);
    close(gcf);
end
%%
% Plot the data

figure
for k = 1:num
    subplot(4,1,k)
    plot(distances{k}, correlations{k}, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);
    % Add a red dashed line at the correlation level of 0.7
    hold on;
    yline(0, '--r',LineWidth=1.5);
    hold off;
    % Set y-axis limits to [0, 1]
    ylim([-1, 1]);
    xlim([0, 100]);
    ax = gca;  % Get handle to current axes
    ax.FontSize = 14;  % Set tick label font size to 14
    % Add labels and title
    plot_title_ = ['Demeaned-Cross-Correlation-' num2str(names{k})];
    title(plot_title_,Fontsize=20);
    xlabel('Distance (km)',Fontsize=15);
    ylabel('Correlation',Fontsize=15);

end
close(gcf);
%% Plot polyfit lines
figure
for k = 1:num
    subplot(1,4,k)
    plot(distances{k}, correlations{k}, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);

    % Polyfit for a linear regression
    p = polyfit(distances{k}, correlations{k}, 1); % 1 is for linear fit
    fittedY = polyval(p, distances{k});
    hold on;
    plot(distances{k}, fittedY, '-b', 'LineWidth', 1.5); % Plot the fitted line

    % Add a red dashed line at the correlation level of 0.5
    xline(20, '--r', 'LineWidth', 1.5);

    ylim([-1, 1]);
    xlim([0, 100]);
    ax = gca;  % Get handle to current axes
    ax.FontSize = 14;  % Set tick label font size to 14

    % Add labels and title
    plot_title_ = [num2str(names{k})];
    title(plot_title_, 'FontSize', 20);
    if k == 2
        xlabel('Distance (km)',Fontsize=15);
    elseif k == 1
        ylabel('Correlation',Fontsize=15);
    end

    % Annotate the gradient (slope) of the linear regression
    gradient_str = sprintf('Gradient: %.4f', p(1));
    text(0.99*ax.XLim(2), 0.97*ax.YLim(2), gradient_str, 'FontSize', 12, 'HorizontalAlignment', 'right');
end
close(gcf);
%% Plot polyfit lines alongtrench
figure
for k = 1:num
    subplot(4,1,k)
    plot(distances{k}, correlations{k}, 'o', 'MarkerSize', 10, 'LineWidth', 1.5);

    % Polyfit for a linear regression
    p = polyfit(distances{k}, correlations{k}, 1); % 1 is for linear fit
    fittedY = polyval(p, distances{k});
    hold on;
    plot(distances{k}, fittedY, '-b', 'LineWidth', 1.5); % Plot the fitted line

    % Add a red dashed line at the correlation level of 0.5
    yline(0, '--r', 'LineWidth', 1.5);

    ylim([-1, 1]);
    xlim([0, 120]);
    ax = gca;  % Get handle to current axes
    ax.FontSize = 14;  % Set tick label font size to 14

    % Add labels and title
    plot_title_ = ['Demeaned-'  num2str(names{k})];
    title(plot_title_, 'FontSize', 20);

    xlabel('Distance (km)',Fontsize=15);

    ylabel('Correlation',Fontsize=15);


    % Annotate the gradient (slope) of the linear regression
    gradient_str = sprintf('Gradient: %.4f', p(1));
    text(0.99*ax.XLim(2), 0.95*ax.YLim(2), gradient_str, 'FontSize', 12, 'HorizontalAlignment', 'right');
end

%%
%%
figure;
hold on;

selected_points = [1];
orangeYellow = [1, 0.7, 0]; % Orange-yellow color
pur = [0.6, 1, 1];
legendHandles = gobjects(0);
legendLabels = {'Acrosstrench-28','Acrosstrench-30','Acrosstrench-32'};

% Define a set of colors for different files
fileColors = { orangeYellow, 'g', 'b', 'm', 'c','r', pur}; % Add or change colors as needed

% Assuming distances and correlations are cell arrays with data from each file
for k = 1:num  % Loop through each file
    [data, ~] = xlsread(filenames{k});
    allData{k} = data;
    currentData = allData{k};
    for idx = 1:length(selected_points)
        i = selected_points(idx);

        % Initialize arrays to hold the x and y values for each selected point and its pairs
        x_selected = [];
        y_selected = [];

        for j = i + 1:size(currentData, 1)
            point_pair_str = sprintf('%d & %d', i, j);
            index = (i-1)*(size(allData{k}, 1)) - (i-1)*i/2 + (j-i);

            % Populating the x_selected and y_selected arrays
            x_selected = [x_selected, distances{k}(index)];
            y_selected = [y_selected, correlations{k}(index)];

            % Assign color based on the file index k
            color = fileColors{k};

            % Plot each point
            h = scatter(distances{k}(index), correlations{k}(index), 100, ... % Size of marker
                'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'LineWidth', 1.5);
            hold on;



            % Add legend entry for the first data point of each file
            if idx == 1 && j == i + 1
                legendHandles(end+1) = h;
                legendLabels{end+1} = sprintf('File %d', k);
            end

            % Annotate each point
            %text(distances{k}(index), correlations{k}(index), point_pair_str, ...
            %'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            %'Color', color, 'FontSize', 14);
        end

        % Plot lines connecting each selected point with its pairs
        %plot(x_selected, y_selected, 'Color', color, 'LineWidth', 1.5);
    end
end


ax = gca;
ax.FontSize = 14;
xlabel('Distance (km)', 'FontSize', 20);
ylabel('Correlation', 'FontSize', 20);
xline(40, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
grid on

legend(legendHandles, legendLabels, 'Location', 'northeast');
xlim([0, 100]);
ylim([-1, 1]);
title('Cross-Correlation Analysis: First Point vs. Remaining Dataset (100-365 days)', 'FontSize', 20);

