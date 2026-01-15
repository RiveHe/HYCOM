% Define file path and variable name
load('E:\HYCOMdata\Chile\2014-2015SSEChile.mat');
load('F:\PhD work\bp_cross\Chile\result\timedate.mat');

%% bandpass filter for data
% order = 4;
% filtT = [3, 100]; % Filter periods in days
% lp_bpt = NaN(size(eff_bpt));
% ftype = 'bandpassiir';
%
% % The sample rate, assuming one sample per day
% sampRate = 8; % eight sample per day
%
% for il = 1:size(eff_bpt,1)
%     % The normalized cutoff frequencies must be between 0 and 0.5, where 0.5 corresponds to the Nyquist frequency
%     highPassCutoffNormalized = (1 / (filtT(1) * 2)); % Normalize by Nyquist frequency
%     lowPassCutoffNormalized = (1 / (filtT(2) * 2));
%
%     % Check if the cutoff frequencies are within the valid range
%     if highPassCutoffNormalized <= 0.5 && lowPassCutoffNormalized >= 0 && highPassCutoffNormalized > lowPassCutoffNormalized
%         % Design the bandpass filter
%         d = designfilt(ftype, 'FilterOrder', order, ...
%             'HalfPowerFrequency1', lowPassCutoffNormalized, 'HalfPowerFrequency2', highPassCutoffNormalized, ...
%             'SampleRate', sampRate);
%
%         % Apply the filter to the data
%         lp_bpt(il,:) = filtfilt(d, eff_bpt(il,:));
%     else
%         warning('Filter frequencies are outside the valid range for il=%d', il);
%         % If the frequencies are invalid, skip this iteration
%         continue;
%     end
% end
%%

% Define target station coordinates
target_coords = [
    -71.767082, -27.197904;
    -71.677419, -27.268363;
    -71.537950, -27.298967;
    -71.348672, -27.329567;
    -71.239088, -27.369796;
    -71.079696, -27.389978;
    -70.930266, -27.419983;
    ];

% Find the index where both latitude and longitude match within the specified tolerance
tolerance = 1e-1; 
indices = zeros(size(target_coords, 1), 1);
for k = 1:size(target_coords, 1)
    target_lon = target_coords(k, 1);
    target_lat = target_coords(k, 2);
    index = find(abs(eff_sta(:, 1) - target_lon) < tolerance & abs(eff_sta(:, 2) - target_lat) < tolerance, 1);
    if ~isempty(index)
        indices(k) = index(1);
    else
        indices(k) = NaN;
    end
end

%% Plot the time series of Bottom Pressure at the specified locations

% Define the start time for the SSE signal
sse_duration_weeks = 5; % duration of the SSE in weeks
% Convert duration to time steps (each time step is 3 hours)
sse_duration_steps = sse_duration_weeks * 7 * 8; % 7 days per week, 8 time steps per day (3-hour intervals)

sse_start_time = datetime(2015, 1, 4, 18, 0, 0); 
sse_start_index = find(Timedate >= sse_start_time, 1);
sse_end_index = sse_start_index + sse_duration_steps - 1;

figure;
hold on;

% Initialize an array to store the standard deviations
std_values = zeros(length(indices), 1);
y_offset = 5; % Define the Y-step offset value

% Plot the time series of Bottom Pressure at the specified locations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Plot the bottom pressure time series with an offset
        plot(Timedate, lp_bpt(indices(k), :) + y_offset * (k - 1), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

        % Calculate the standard deviation for the current time series
        std_values(k) = std(lp_bpt(indices(k), :));
    end
end

% Add dashed lines to indicate the start and end of the SSE ramp
xline(Timedate(sse_start_index), '--k', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
xline(Timedate(sse_end_index), '--k', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

% Add scale bar
y_position = max(max(lp_bpt(indices, :))+20); % Set the y position for the scale bar
plot([max(Timedate)-25 max(Timedate)-25], [y_position y_position+5], 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Scale bar
text(max(Timedate)-23, y_position+1, '5 hPa \DeltaP', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k'); % Label

% Set the labels and title of the plot
ylim([-10 50])
xlabel('Time');
ylabel('Bottom Pressure (hPa)');
title('Bottom Pressure Time Series of HYCOM output');
l1 = legend('show');
l1.Location = 'northwest';
ax = gca;
ax.FontSize = 16;

offset = 1.2;
% Annotate the plot with the standard deviations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Position the text annotation at the end of the time series plot
        text(max(Timedate), y_offset * (k - 1), sprintf('STD Location %d: %.2f hPa', k, std_values(k)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
    end
end

hold off;



% Save the figure
%saveas(gcf, 'bottom_pressure_timeseries.png');
%% add sse to original signals

% Define the SSE signal parameters
sse_duration_weeks = 5; % duration of the SSE in weeks
sse_amplitude_cm = 1; % amplitude of the SSE in cm
sse_amplitude_hPa = sse_amplitude_cm; % convert cm to hPa

% Convert duration to time steps (each time step is 3 hours)
sse_duration_steps = sse_duration_weeks * 7 * 8; % 7 days per week, 8 time steps per day (3-hour intervals)

% Generate the SSE ramp signal
sse_signal = linspace(0, sse_amplitude_hPa, sse_duration_steps);

% Define the start time for the SSE signal
sse_start_time = datetime(2015, 1, 4, 18, 0, 0); %  start time yyyy-m-d-hh-mm-ss
sse_start_index = find(Timedate >= sse_start_time, 1);
sse_end_index = sse_start_index + sse_duration_steps - 1;

% Add 1 cm signal after the end of the SSE
post_sse_signal = ones(1, length(Timedate) - sse_end_index) * sse_amplitude_cm;

% Ensure the end index does not exceed the time series length
if sse_end_index > length(Timedate)
    sse_end_index = length(Timedate);
    sse_signal = linspace(0, sse_amplitude_hPa, sse_end_index - sse_start_index + 1);
end

% Plot the time series of Bottom Pressure at the specified locations
y_offset = 5;
figure;
hold on;

% Initialize an array to store the standard deviations
std_values = zeros(length(indices), 1);

% Plot the time series of Bottom Pressure at the specified locations with SSE signal
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Corrected bottom pressure with SSE signal
        corrected_bp_with_sse = lp_bpt(indices(k), :);
        corrected_bp_with_sse(sse_start_index:sse_end_index) = ...
            corrected_bp_with_sse(sse_start_index:sse_end_index) + sse_signal;
        % Add the post-SSE signal
        corrected_bp_with_sse(sse_end_index+1:end) = ...
            corrected_bp_with_sse(sse_end_index+1:end) + post_sse_signal;

        % Plot the bottom pressure time series
        plot(Timedate, corrected_bp_with_sse+ y_offset * (k - 1), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

        % Calculate the standard deviation for the current time series
        std_values(k) = std(corrected_bp_with_sse);
    end
end

% Add dashed lines to indicate the start and end of the SSE ramp
xline(Timedate(sse_start_index), '--k', 'Start of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
xline(Timedate(sse_end_index), '--k', 'End of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

% Add scale bar
y_position = max(max(lp_bpt(indices, :))+20); % Set the y position for the scale bar
plot([max(Timedate)-25 max(Timedate)-25], [y_position y_position+5], 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Scale bar
text(max(Timedate)-23, y_position+1, '5 hPa \DeltaP', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k'); % Label

xlabel('Time');
ylabel('Bottom Pressure (hPa)');
title('Bottom Pressure Time Series of HYCOM output with SSE');
l1 = legend('show');
l1.Location = 'northwest';
ax = gca; 
ax.FontSize = 16;

offset = 1.2;
% Annotate the plot with the standard deviations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Position the text annotation at the end of the time series plot
        text(max(Timedate), y_offset * (k - 1), sprintf('STD Location %d: %.2f hPa', k, std_values(k)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
    end
end

hold off;

% Save the figure
saveas(gcf, 'bottom_pressure_corrected_timeseries_with_sse.png');

%% corrected by reference stations
% If change the reference station, also need to change the figure title
% ref_coords = [
%     -71.960000, -27.657000;
%     -71.889000, -27.730000;
%     -71.722000, -27.798000;
%     -71.498000, -27.831000;
%     -71.390000, -27.834000;
%     -71.198000, -27.857000;
%     -71.110000, -27.878000
% ];%50km

ref_coords = [
    -71.495000, -26.269000;
    -71.439000, -26.308000;
    -71.251000, -26.335000;
    -71.172000, -26.353000;
    -71.084000, -26.422000;
    -70.887000, -26.482000;
    -70.708000, -26.524000
    ];%100km

% Find the index where both latitude and longitude match within the specified tolerance
tolerance = 1e-1; 
target_indices = zeros(size(target_coords, 1), 1);
ref_indices = zeros(size(ref_coords, 1), 1);

for k = 1:size(target_coords, 1)
    target_lon = target_coords(k, 1);
    target_lat = target_coords(k, 2);
    ref_lon = ref_coords(k, 1);
    ref_lat = ref_coords(k, 2);

    target_index = find(abs(eff_sta(:, 1) - target_lon) < tolerance & abs(eff_sta(:, 2) - target_lat) < tolerance, 1);
    ref_index = find(abs(eff_sta(:, 1) - ref_lon) < tolerance & abs(eff_sta(:, 2) - ref_lat) < tolerance, 1);

    if ~isempty(target_index)
        target_indices(k) = target_index(1);
    else
        target_indices(k) = NaN; 
    end

    if ~isempty(ref_index)
        ref_indices(k) = ref_index(1);
    else
        ref_indices(k) = NaN;
    end
end
%%
% Save each eff_sta(ref_indices) and eff_sta(target_indices) to an Excel file
% filename = 'eff_sta_data_50km.xlsx';
%
% for k = 1:length(target_indices)
%     if ~isnan(target_indices(k)) && ~isnan(ref_indices(k))
%         target_data = eff_sta(target_indices(k), :);
%         ref_data = eff_sta(ref_indices(k), :);
%
%         % Write target data to Excel
%         sheet_name = 'Eff_Sta_Data';
%         writematrix(target_data, filename, 'Sheet', sheet_name, 'Range', sprintf('A%d', k * 2 - 1));
%         writematrix(ref_data, filename, 'Sheet', sheet_name, 'Range', sprintf('A%d', k * 2));
%     end
% end

% Plot the time series of Bottom Pressure corrected by the reference station
y_offset = 5; % Define the Y-step offset value
figure;
hold on;
std_values = [];
for k = 1:length(target_indices)
    if ~isnan(target_indices(k)) && ~isnan(ref_indices(k))
        corrected_bp = lp_bpt(target_indices(k), :) - lp_bpt(ref_indices(k), :);
        plot(Timedate, corrected_bp + y_offset * (k - 1), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);
        % Calculate the standard deviation
        std_val = std(corrected_bp);
        std_values = [std_values, std_val];
    end
end

% Add dashed lines to indicate the start and end of the SSE ramp
xline(Timedate(sse_start_index), '--k', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
xline(Timedate(sse_end_index), '--k', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

y_position = max(max(lp_bpt(indices, :))+20); % Set the y position for the scale bar
plot([max(Timedate)-25 max(Timedate)-25], [y_position y_position+5], 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Scale bar
text(max(Timedate)-23, y_position+1, '5 hPa \DeltaP', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k'); % Label

xlabel('Time');
ylabel('Bottom Pressure (hPa)');
title('Bottom Pressure Time Series Corrected by Reference Station in 100km');
l1 = legend('show');
l1.Location = 'northwest';
ax = gca; 
%ylim([-6 6])
ax.FontSize = 16;

offset = 1.2;
% Annotate the plot with the standard deviations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Position the text annotation at the end of the time series plot
        text(max(Timedate), y_offset * (k - 1), sprintf('STD Location %d: %.2f hPa', k, std_values(k)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
    end
end


% Save the figure
saveas(gcf, 'bottom_pressure_corrected_timeseries.png');
%% add SSE singal to the time series

% Define the SSE signal parameters
sse_duration_weeks = 5; % duration of the SSE in weeks
sse_amplitude_hPa = 1; % convert cm to hPa

% Convert duration to time steps (each time step is 3 hours)
sse_duration_steps = sse_duration_weeks * 7 * 8; % 7 days per week, 8 time steps per day (3-hour intervals)

% Generate the SSE ramp signal
sse_signal = linspace(0, sse_amplitude_hPa, sse_duration_steps);


% Define the start time for the SSE signal
sse_start_time = datetime(2015, 1, 4, 18, 0, 0); 
sse_start_index = find(Timedate >= sse_start_time, 1);
sse_end_index = sse_start_index + sse_duration_steps - 1;

% Add 1 cm signal after the end of the SSE
post_sse_signal = ones(1, length(Timedate) - sse_end_index) * sse_amplitude_cm;

y_offset = 5; % Define the Y-step offset value
figure;
hold on;
std_values = [];
for k = 1:length(target_indices)
    if ~isnan(target_indices(k)) && ~isnan(ref_indices(k))
        corrected_bp = lp_bpt(target_indices(k), :) - lp_bpt(ref_indices(k), :);

        % Add the SSE signal to the corrected bottom pressure time series within the correct duration
        corrected_bp_with_sse = corrected_bp;
        corrected_bp_with_sse(sse_start_index:sse_end_index) = ...
            corrected_bp_with_sse(sse_start_index:sse_end_index) + sse_signal;
        % Add the post-SSE signal
        corrected_bp_with_sse(sse_end_index+1:end) = ...
            corrected_bp_with_sse(sse_end_index+1:end) + post_sse_signal;

        plot(Timedate, corrected_bp_with_sse + y_offset * (k - 1), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

        % Calculate the standard deviation
        std_val = std(corrected_bp_with_sse);
        std_values = [std_values, std_val];
    end
end



% Add dashed lines to indicate the start and end of the SSE ramp
xline(Timedate(sse_start_index), '--k', 'Start of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
xline(Timedate(sse_end_index), '--k', 'End of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

% Add scale bar
y_position = max(max(lp_bpt(indices, :))+20); % Set the y position for the scale bar
plot([max(Timedate)-25 max(Timedate)-25], [y_position y_position+5], 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Scale bar
text(max(Timedate)-23, y_position+1, '5 hPa \DeltaP', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k'); % Label

xlabel('Time');
ylabel('Bottom Pressure (hPa)');
title('Bottom Pressure Time Series Corrected by Reference Station in 100km with SSE');
l1 = legend('show');
l1.Location = 'northwest';
ax = gca; % Get current axis
ax.FontSize = 16;
offset = 1.2;
% Annotate the plot with the standard deviations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Position the text annotation at the end of the time series plot
        text(max(Timedate), y_offset * (k - 1), sprintf('STD Location %d: %.2f hPa', k, std_values(k)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
    end
end

hold off;
%saveas(gcf, 'bottom_pressure_corrected_timeseries_with_sse.png');

%% decomposite the corrected time series into original and reference locations

figure;
std_values = [];
subplot(211)
for k = 4 % change this to plot each location
    hold on
    if ~isnan(target_indices(k)) && ~isnan(ref_indices(k))
        corrected_bp = (lp_bpt(target_indices(k), :)) - (lp_bpt(ref_indices(k), :));
        plot(Timedate, corrected_bp, 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

        % Calculate the standard deviation
        std_val = std(corrected_bp);
        std_values = [std_values, std_val];
    end
    offset = 1.2;
    % Annotate the plot with the standard deviations


    % Position the text annotation at the end of the time series plot
    text(max(Timedate)-20,  0.5, sprintf('STD Location %d: %.2f hPa', k, std_values), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
end

l1 = legend(sprintf('corrected location %d',k));
l1.Location = 'northwest';
title(sprintf('100km for station %d', k))
subplot(212)
hold on
if ~isnan(target_indices(k)) && ~isnan(ref_indices(k))
    plot(Timedate, lp_bpt(target_indices(k), :), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);
    plot(Timedate, lp_bpt(ref_indices(k), :), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

    % Calculate the standard deviation
    std_val = std(corrected_bp);
    std_values = [std_values, std_val];
end
l1 = legend(sprintf('Location %d', k), sprintf('Reference station for Location %d', k));
l1.Location = 'northwest';

xlabel('Time');
ylabel('Bottom Pressure (hPa)');
%% PSD of SSE +HYCOM
% Define target station coordinates
target_coords = [
    -71.767082, -27.197904;
    -71.677419, -27.268363;
    -71.537950, -27.298967;
    -71.348672, -27.329567;
    -71.239088, -27.369796;
    -71.079696, -27.389978;
    -70.930266, -27.419983;
];

% Find the index where both latitude and longitude match within the specified tolerance
tolerance = 1e-1;
indices = zeros(size(target_coords, 1), 1);
for k = 1:size(target_coords, 1)
    target_lon = target_coords(k, 1);
    target_lat = target_coords(k, 2);
    index = find(abs(eff_sta(:, 1) - target_lon) < tolerance & abs(eff_sta(:, 2) - target_lat) < tolerance, 1);
    if ~isempty(index)
        indices(k) = index(1);
    else
        indices(k) = NaN;
    end
end

% Define the SSE signal parameters
sse_duration_weeks = 5; % duration of the SSE in weeks
sse_amplitude_cm = 1; % amplitude of the SSE in cm
sse_amplitude_hPa = sse_amplitude_cm; % convert cm to hPa

% Convert duration to time steps (each time step is 3 hours)
sse_duration_steps = sse_duration_weeks * 7 * 8; % 7 days per week, 8 time steps per day (3-hour intervals)

% Generate the SSE ramp signal
sse_signal = linspace(0, sse_amplitude_hPa, sse_duration_steps);

% Define the start time for the SSE signal
sse_start_time = datetime(2015, 1, 4, 18, 0, 0); % start time yyyy-m-d-hh-mm-ss
sse_start_index = find(Timedate >= sse_start_time, 1);
sse_end_index = sse_start_index + sse_duration_steps - 1;

% Add 1 cm signal after the end of the SSE
post_sse_signal = ones(1, length(Timedate) - sse_end_index) * sse_amplitude_cm;

% Ensure the end index does not exceed the time series length
if sse_end_index > length(Timedate)
    sse_end_index = length(Timedate);
    sse_signal = linspace(0, sse_amplitude_hPa, sse_end_index - sse_start_index + 1);
end

% Plot the time series of Bottom Pressure at the specified locations with SSE signal
y_offset = 5;
figure;
hold on;

% Initialize an array to store the standard deviations
std_values = zeros(length(indices), 1);

% Plot the time series of Bottom Pressure at the specified locations with SSE signal
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Corrected bottom pressure with SSE signal
        corrected_bp_with_sse = lp_bpt(indices(k), :);
        corrected_bp_with_sse(sse_start_index:sse_end_index) = ...
            corrected_bp_with_sse(sse_start_index:sse_end_index) + sse_signal;
        % Add the post-SSE signal
        corrected_bp_with_sse(sse_end_index+1:end) = ...
            corrected_bp_with_sse(sse_end_index+1:end) + post_sse_signal;

        % Plot the bottom pressure time series
        plot(Timedate, corrected_bp_with_sse + y_offset * (k - 1), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

        % Calculate the standard deviation for the current time series
        std_values(k) = std(corrected_bp_with_sse);
    end
end

% Add dashed lines to indicate the start and end of the SSE ramp
xline(Timedate(sse_start_index), '--k', 'Start of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
xline(Timedate(sse_end_index), '--k', 'End of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

% Add scale bar
y_position = max(max(lp_bpt(indices, :))+20); % Set the y position for the scale bar
plot([max(Timedate)-25 max(Timedate)-25], [y_position y_position+5], 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Scale bar
text(max(Timedate)-23, y_position+1, '5 hPa \DeltaP', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k'); % Label

xlabel('Time');
ylabel('Bottom Pressure (hPa)');
title('Bottom Pressure Time Series of HYCOM output with SSE');
l1 = legend('show');
l1.Location = 'northwest';
ax = gca;
ax.FontSize = 16;

offset = 1.2;
% Annotate the plot with the standard deviations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Position the text annotation at the end of the time series plot
        text(max(Timedate), y_offset * (k - 1), sprintf('STD Location %d: %.2f hPa', k, std_values(k)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
    end
end

hold off;

% Calculate and plot PSD for the corrected bottom pressure with SSE
figure;
hold on;
for k = 1:length(indices)
    if ~isnan(indices(k))
        corrected_bp_with_sse = lp_bpt(indices(k), :);
        corrected_bp_with_sse(sse_start_index:sse_end_index) = ...
            corrected_bp_with_sse(sse_start_index:sse_end_index) + sse_signal;
        % Add the post-SSE signal
        corrected_bp_with_sse(sse_end_index+1:end) = ...
            corrected_bp_with_sse(sse_end_index+1:end) + post_sse_signal;

        % Calculate PSD
        [pxx, f] = pwelch(corrected_bp_with_sse, [], [], [], 1/10800); % 3 hours = 10800 seconds
        plot(f, 10*log10(pxx), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);
    end
end

% Set the labels and title of the plot
xlabel('Frequency (Hz)');
ylabel('Noise Spectrum (dB relative to 1 (m/s)^2/Hz^{1/2})');
title('Noise Spectra of HYCOM output with SSE');
legend('show');
ax = gca;
ax.FontSize = 16;
hold off;
grid on;
%% PSD analysis 

% Define target station coordinates
target_coords = [
    -71.767082, -27.197904;
    -71.677419, -27.268363;
    -71.537950, -27.298967;
    -71.348672, -27.329567;
    -71.239088, -27.369796;
    -71.079696, -27.389978;
    -70.930266, -27.419983;
];

% Find the index where both latitude and longitude match within the specified tolerance
tolerance = 1e-1; 
indices = zeros(size(target_coords, 1), 1);
for k = 1:size(target_coords, 1)
    target_lon = target_coords(k, 1);
    target_lat = target_coords(k, 2);
    index = find(abs(eff_sta(:, 1) - target_lon) < tolerance & abs(eff_sta(:, 2) - target_lat) < tolerance, 1);
    if ~isempty(index)
        indices(k) = index(1);
    else
        indices(k) = NaN;
    end
end

% Define the start time for the SSE signal
sse_duration_weeks = 5; % duration of the SSE in weeks
% Convert duration to time steps (each time step is 3 hours)
sse_duration_steps = sse_duration_weeks * 7 * 8; % 7 days per week, 8 time steps per day (3-hour intervals)

sse_start_time = datetime(2015, 1, 4, 18, 0, 0); 
sse_start_index = find(Timedate >= sse_start_time, 1);
sse_end_index = sse_start_index + sse_duration_steps - 1;

% Plot the time series of Bottom Pressure at the specified locations
figure;
hold on;

% Initialize an array to store the standard deviations
std_values = zeros(length(indices), 1);
y_offset = 5; % Define the Y-step offset value

% Plot the time series of Bottom Pressure at the specified locations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Plot the bottom pressure time series with an offset
        plot(Timedate, lp_bpt(indices(k), :) + y_offset * (k - 1), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);

        % Calculate the standard deviation for the current time series
        std_values(k) = std(lp_bpt(indices(k), :));
    end
end

% Add dashed lines to indicate the start and end of the SSE ramp
xline(Timedate(sse_start_index), '--k', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
xline(Timedate(sse_end_index), '--k', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

% Add scale bar
y_position = max(max(lp_bpt(indices, :))+20); % Set the y position for the scale bar
plot([max(Timedate)-25 max(Timedate)-25], [y_position y_position+5], 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Scale bar
text(max(Timedate)-23, y_position+1, '5 hPa \DeltaP', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k'); % Label

% Set the labels and title of the plot
ylim([-10 50])
xlabel('Time');
ylabel('Bottom Pressure (hPa)');
title('Bottom Pressure Time Series of HYCOM output');
l1 = legend('show');
l1.Location = 'northwest';
ax = gca;
ax.FontSize = 16;

offset = 1.2;
% Annotate the plot with the standard deviations
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Position the text annotation at the end of the time series plot
        text(max(Timedate), y_offset * (k - 1), sprintf('STD Location %d: %.2f hPa', k, std_values(k)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');
    end
end

hold off;

% Compute and plot the noise spectra
figure;
hold on;
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Compute the power spectral density
        [pxx_without_sse, f] = pwelch(lp_bpt(indices(k), :), [], [], [], 1/(3*3600)); % 3-hour intervals

        % Plot the power spectral density
        loglog(f, 10*log10(pxx_without_sse), 'DisplayName', sprintf('Location %d', k), 'LineWidth', 1.5);
    end
end

xlabel('Frequency (Hz)');
ylabel('Noise spectrum (dB relative to 1 (m/s)/Hz^{1/2})');
title('Noise Spectra of HYCOM output');
legend('show');
grid on;
ax = gca;
ax.FontSize = 16;
hold off;
%%
%% diff between SSE+HYCOM and HYCOM 
% Assuming psd_without_sse and psd_with_sse are the PSD values
% and frequency is the corresponding frequency array

% Perform paired t-test
[h, p] = ttest(pxx_without_sse, pxx);

% Display the p-value
disp(['P-value: ', num2str(p)]);

% Compute the difference in PSD values
psd_diff = pxx - pxx_without_sse;

% Plot the difference in PSD
figure;
plot(f, psd_diff, 'k', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Difference in Noise Spectrum (dB)');
title('Difference in Noise Spectrum between HYCOM with SSE and HYCOM');
grid on
ax = gca;
ax.FontSize = 16;
xlim([0 1e-5]);

% Add the p-value as text on the plot
text(0.9e-5, 0.9 * max(psd_diff), ['P-value: ', num2str(p)], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'b');




%% time-frequency analysis 


% Define target station coordinates
target_coords = [
    -71.767082, -27.197904;
    -71.677419, -27.268363;
    -71.537950, -27.298967;
    -71.348672, -27.329567;
    -71.239088, -27.369796;
    -71.079696, -27.389978;
    -70.930266, -27.419983;
    ];

% Find the index where both latitude and longitude match within the specified tolerance
tolerance = 1e-1;
indices = zeros(size(target_coords, 1), 1);
for k = 1:size(target_coords, 1)
    target_lon = target_coords(k, 1);
    target_lat = target_coords(k, 2);
    index = find(abs(eff_sta(:, 1) - target_lon) < tolerance & abs(eff_sta(:, 2) - target_lat) < tolerance, 1);
    if ~isempty(index)
        indices(k) = index(1);
    else
        indices(k) = NaN;
    end
end

% Define the SSE signal parameters
sse_duration_weeks = 5; % duration of the SSE in weeks
sse_amplitude_cm = 1; % amplitude of the SSE in cm
sse_amplitude_hPa = sse_amplitude_cm; % convert cm to hPa

% Convert duration to time steps (each time step is 3 hours)
sse_duration_steps = sse_duration_weeks * 7 * 8; % 7 days per week, 8 time steps per day (3-hour intervals)

% Generate the SSE ramp signal
sse_signal = linspace(0, sse_amplitude_hPa, sse_duration_steps);

% Define the start time for the SSE signal
sse_start_time = datetime(2015, 1, 4, 18, 0, 0); %  start time yyyy-m-d-hh-mm-ss
sse_start_index = find(Timedate >= sse_start_time, 1);
sse_end_index = sse_start_index + sse_duration_steps - 1;

% Add 1 cm signal after the end of the SSE
post_sse_signal = ones(1, length(Timedate) - sse_end_index) * sse_amplitude_cm;

% Ensure the end index does not exceed the time series length
if sse_end_index > length(Timedate)
    sse_end_index = length(Timedate);
    sse_signal = linspace(0, sse_amplitude_hPa, sse_end_index - sse_start_index + 1);
end

% Convert datetime to day of year for x-axis
day_of_year = day(Timedate, 'dayofyear');

% Perform time-frequency analysis for each location
for k = 1:length(indices)
    if ~isnan(indices(k))
        % Original signal
        signal_orig = lp_bpt(indices(k), :);

        % Corrected signal with SSE
        corrected_bp_with_sse = lp_bpt(indices(k), :);
        corrected_bp_with_sse(sse_start_index:sse_end_index) = ...
            corrected_bp_with_sse(sse_start_index:sse_end_index) + sse_signal;
        corrected_bp_with_sse(sse_end_index+1:end) = ...
            corrected_bp_with_sse(sse_end_index+1:end) + post_sse_signal;

        % Plot spectrogram for original signal

        figure;
        subplot(2,1,1);
        spectrogram(signal_orig, 256, 250, 256, 1/(3*3600), 'yaxis');
        title(sprintf('Spectrogram of Original Signal at Location %d', k));
        xlabel('Time');
        ylabel('Frequency (Hz)');
        xline((sse_start_index-1)*3/24, '--k', 'Start of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
        xline((sse_end_index-1)*3/24, '--k', 'End of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
        ylim([0 1e-5])
        caxis([40,80])

        % Plot spectrogram of the signal with SSE
        subplot(2,1,2);
        spectrogram(corrected_bp_with_sse, 256, 250, 256, 1/(3*3600), 'yaxis');
        title(sprintf('Spectrogram of Signal with SSE at Location %d', k));
        colorbar;
        xlabel('Time');
        ylabel('Frequency (Hz)');
        xline((sse_start_index-1)*3/24, '--k', 'Start of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
        xline((sse_end_index-1)*3/24, '--k', 'End of SSE', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
        ylim([0 1e-5])
        caxis([40,80])

    end
end





