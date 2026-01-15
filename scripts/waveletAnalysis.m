%%
for il = 1:size(eff_bpt_combined,1)
    if any(isnan(eff_bpt_combined(il,:)))
        % Linearly interpolate over NaNs for each row
        data_no_nan = eff_bpt_combined(il,:);
        nan_indices = isnan(data_no_nan);
        data_no_nan(nan_indices) = interp1(find(~nan_indices), data_no_nan(~nan_indices), find(nan_indices), 'linear', 'extrap');
        
        % Replace the original row with the interpolated data
        eff_bpt_combined(il,:) = data_no_nan;
    end
end
%%
% Assuming 'eff_bpt' is your loaded data matrix where rows represent different
% observation points/stations and columns represent time series data points.

% Choose the mother wavelet and the number of decomposition levels
motherWavelet = 'db40'; % Daubechies wavelet of order 40
N = 8; % Number of decomposition levels
originalSamplingFrequency = 1 / 3; % samples per hour (3 hours per sample)

% Calculate the Nyquist frequency for the original sampling rate
nyquistFrequency = originalSamplingFrequency / 2;

% Preallocate arrays for the reconstructed components
A = zeros(size(eff_bpt_combined)); % Approximation coefficients at each level
D = zeros([size(eff_bpt_combined), N]); % Detail coefficients at each level

% Calculate frequency range for each level
frequencyRanges = zeros(N, 2); % Initialize matrix to hold frequency ranges
cyclesPerHourToDays = 1 / 24;
for level = 1:N
     % Calculate the lower and upper frequency limits for each level
    lowerFrequency = nyquistFrequency / (2^level);
    upperFrequency = nyquistFrequency / (2^(level-1));
    
    % Convert frequency limits to periods in hours
    lowerPeriodHours = 1 / lowerFrequency;
    upperPeriodHours = 1 / upperFrequency;
    
    % Convert periods in hours to days
    lowerPeriodDays = lowerPeriodHours * cyclesPerHourToDays;
    upperPeriodDays = upperPeriodHours * cyclesPerHourToDays;
    
    % Store the period range in days in the frequencyRanges matrix
    frequencyRanges(level, :) = [lowerPeriodDays, upperPeriodDays];
end

% Perform DWT for each observation point/station
for i = sta
% for i = 1:size(eff_sta,1)
    [C, L] = wavedec(eff_bpt_combined(i, :), N, motherWavelet);
    
    % Reconstruct the approximation and detail components separately
    A(i, :) = wrcoef('a', C, L, motherWavelet, N);
    for j = 1:N
        D(i, :, j) = wrcoef('d', C, L, motherWavelet, j);
    end
    disp(['Finish station: ', num2str(i), ' out of ', num2str(size(eff_sta, 1))]);
end
%%
% Plotting the original signal and the components for the first observation point as an example
 % Replace with actual time vector if available
figure;
subplot(N+2, 1, 1);
plot(combineTime, eff_bpt_combined(sta, :));
title('Original Signal');

subplot(N+2, 1, 2);
plot(combineTime, A(sta, :));
title('Approximation Component');

for j = 1:N
    subplot(N+2, 1, j+2);
    plot(combineTime, D(sta, :, j));
    title(['Detail Component at Level ', num2str(j), ' (', num2str(frequencyRanges(j, 1)), '-', num2str(frequencyRanges(j, 2)), ' days)']);
end

% Configure the layout
xlabel('Time');
ylabel('Pressure');
%%
%filepath = '/Volumes/My Data/PhD work/bp_cross/Cascadia/data/ssh/1.nc';
%timeData = ncread(filepath, 'time');
% Perform Continuous Wavelet Transform (CWT) on the data
load('/Volumes/My Data/PhD work/bp_cross/Cascadia/result/timedate.mat', 'Timedate');
sta = 7941; 
sta2 = 7402;
sta3 = 6770;


[cfs, frequencies] = cwt(eff_bpt_combined(sta, :), originalSamplingFrequency);
[cfs2, frequencies2] = cwt(eff_bpt_combined(sta2, :), originalSamplingFrequency);
[cfs3, frequencies3] = cwt(eff_bpt_combined(sta3, :), originalSamplingFrequency);

% Calculate the power spectrum, transpose it to match the time and frequency dimensions
powerSpectrum = abs(cfs).^2;
powerSpectrum2 = abs(cfs2).^2;
powerSpectrum3 = abs(cfs3).^2;

% Ensure the dimensions match the expectations for contourf
% The powerSpectrum should be transposed since MATLAB expects the rows to correspond to the y-axis (frequencies)


timeAxis = datenum(combineTime);


%%
figure;
subplot(311)
contourf(timeAxis, frequencies, powerSpectrum, 'LineColor', 'none');
% Adjust the x-axis to display datetime values
% Choose the date and time format for the x-axis labels
dateFormat = 'yy/mm/dd'; % This is an example format, change according to your needs
ax = gca; % Get current axis
ax.FontSize = 14;
ax.XTick = linspace(min(timeAxis), max(timeAxis), 12); % Set 10 ticks on x-axis
datetick(ax, 'x', dateFormat, 'keepticks', 'keeplimits'); % Apply the date format to the x-axis ticks
caxis([0,5]);
ylim([3.9970e-04,0.04]); 
title('Wavelet Spectrum for Station 1');
xlabel('Time');
ylabel('Frequency (Hz)');
% Adding a horizontal color bar and positioning it
cb = colorbar;
cb.Location = 'eastoutside'; % This places the color bar at the top of the plot
cb.Orientation = 'vertical';
% Manually adjust the color bar position to move it towards the northeast corner
%cb.Position = [0.7 0.89 0.02 1]; 

subplot(312)
contourf(timeAxis, frequencies2, powerSpectrum2, 'LineColor', 'none');
% Adjust the x-axis to display datetime values
% Choose the date and time format for the x-axis labels
dateFormat = 'mm/dd'; % This is an example format, change according to your needs
ax = gca; % Get current axis
ax.FontSize = 14;
ax.XTick = linspace(min(timeAxis), max(timeAxis), 12); % Set 10 ticks on x-axis
datetick(ax, 'x', dateFormat, 'keepticks', 'keeplimits'); % Apply the date format to the x-axis ticks
caxis([0,5]);
ylim([3.9970e-04,0.04]); 
title('Wavelet Spectrum for Station 2');
xlabel('Time');
ylabel('Frequency (Hz)');
cb = colorbar;
cb.Location = 'eastoutside'; % This places the color bar at the top of the plot
cb.Orientation = 'vertical';
%colorbar; % Adds a color bar to indicate the power spectrum values

subplot(313)
contourf(timeAxis, frequencies3, powerSpectrum3, 'LineColor', 'none');
% Adjust the x-axis to display datetime values
% Choose the date and time format for the x-axis labels
dateFormat = 'mm/dd'; % This is an example format, change according to your needs
ax = gca; % Get current axis
ax.FontSize = 14;
ax.XTick = linspace(min(timeAxis), max(timeAxis), 12); % Set 10 ticks on x-axis
datetick(ax, 'x', dateFormat, 'keepticks', 'keeplimits'); % Apply the date format to the x-axis ticks
caxis([0,5]);
ylim([3.9970e-04,0.04]); 
title('Wavelet Spectrum for Station 3');
xlabel('Time');
ylabel('Frequency (Hz)');
cb = colorbar;
cb.Location = 'eastoutside'; % This places the color bar at the top of the plot
cb.Orientation = 'vertical';
%colorbar; % Adds a color bar to indicate the power spectrum values

%%
% Assuming 'eff_v2' is your loaded data matrix where rows represent different
% observation points/stations and columns represent time series data points.

% Choose the mother wavelet and the number of decomposition levels
motherWavelet = 'db40'; % Daubechies wavelet of order 40
N = 8; % Number of decomposition levels
originalSamplingFrequency = 1 / 3; % samples per hour (3 hours per sample)

% Calculate the Nyquist frequency for the original sampling rate
nyquistFrequency = originalSamplingFrequency / 2;

% Preallocate arrays for the reconstructed components
A1 = zeros(size(eff_v2)); % Approximation coefficients at each level
D1 = zeros([size(eff_v2), N]); % Detail coefficients at each level

% Calculate frequency range for each level
frequencyRanges = zeros(N, 2); % Initialize matrix to hold frequency ranges
cyclesPerHourToDays = 1 / 24;
for level = 1:N
     % Calculate the lower and upper frequency limits for each level
    lowerFrequency = nyquistFrequency / (2^level);
    upperFrequency = nyquistFrequency / (2^(level-1));
    
    % Convert frequency limits to periods in hours
    lowerPeriodHours = 1 / lowerFrequency;
    upperPeriodHours = 1 / upperFrequency;
    
    % Convert periods in hours to days
    lowerPeriodDays = lowerPeriodHours * cyclesPerHourToDays;
    upperPeriodDays = upperPeriodHours * cyclesPerHourToDays;
    
    % Store the period range in days in the frequencyRanges matrix
    frequencyRanges(level, :) = [lowerPeriodDays, upperPeriodDays];
end

% Perform DWT for each observation point/station
for i = 1:10
    [C1, L1] = wavedec(eff_v2(i, :), N, motherWavelet);
    
    % Reconstruct the approximation and detail components separately
    A1(i, :) = wrcoef('a', C1, L1, motherWavelet, N);
    for j = 1:N
        D1(i, :, j) = wrcoef('d', C1, L1, motherWavelet, j);
    end
end
%%
% Assuming 'eff_u2' is your loaded data matrix where rows represent different
% observation points/stations and columns represent time series data points.

% Choose the mother wavelet and the number of decomposition levels
motherWavelet = 'db40'; % Daubechies wavelet of order 40
N = 8; % Number of decomposition levels
originalSamplingFrequency = 1 / 3; % samples per hour (3 hours per sample)

% Calculate the Nyquist frequency for the original sampling rate
nyquistFrequency = originalSamplingFrequency / 2;

% Preallocate arrays for the reconstructed components
A2 = zeros(size(eff_u2)); % Approximation coefficients at each level
D2 = zeros([size(eff_u2), N]); % Detail coefficients at each level

% Calculate frequency range for each level
frequencyRanges = zeros(N, 2); % Initialize matrix to hold frequency ranges
cyclesPerHourToDays = 1 / 24;
for level = 1:N
     % Calculate the lower and upper frequency limits for each level
    lowerFrequency = nyquistFrequency / (2^level);
    upperFrequency = nyquistFrequency / (2^(level-1));
    
    % Convert frequency limits to periods in hours
    lowerPeriodHours = 1 / lowerFrequency;
    upperPeriodHours = 1 / upperFrequency;
    
    % Convert periods in hours to days
    lowerPeriodDays = lowerPeriodHours * cyclesPerHourToDays;
    upperPeriodDays = upperPeriodHours * cyclesPerHourToDays;
    
    % Store the period range in days in the frequencyRanges matrix
    frequencyRanges(level, :) = [lowerPeriodDays, upperPeriodDays];
end

% Perform DWT for each observation point/station
for i = 1:10
    [C2, L2] = wavedec(eff_u2(i, :), N, motherWavelet);
    
    % Reconstruct the approximation and detail components separately
    A2(i, :) = wrcoef('a', C2, L2, motherWavelet, N);
    for j = 1:N
        D2(i, :, j) = wrcoef('d', C2, L2, motherWavelet, j);
    end
end
%%
% 0 pressure 1 is V and 2 is U
figure;

subplot(N+2, 1, 1);
plot(Timedate, eff_bpt(1, :),'b');
hold on; 
plot(Timedate, eff_v2(1, :),'r');
plot(Timedate, eff_u2(1, :),'g');
hold off
title('Original Signal(blue:pressure red:V green:U)',FontSize=15);
xlabel('Time');
ylabel('hPa/(cm/s)');

subplot(N+2, 1, 2);
plot(Timedate, A(1, :),'b');
hold on;
plot(Timedate, A1(1, :),'r');
plot(Timedate, A2(1, :),'g');
ylabel('hPa/(cm/s)');
title('Approximation Component',FontSize=15);

for j = 1:N
    subplot(N+2, 1, j+2);
    plot(Timedate, D(1, :, j),'b');
    hold on;
    plot(Timedate, D1(1, :, j),'r');
    plot(Timedate, D2(1, :, j),'g');

    xlabel('Time');
    ylabel('hPa/(cm/s)');
    title(['Detail Component at Level ', num2str(j), ' (', num2str(frequencyRanges(j, 1)), '-', num2str(frequencyRanges(j, 2)), ' days)'],FontSize=15);
end

% Configure the layout


