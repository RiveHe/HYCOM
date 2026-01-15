% A = '/Volumes/My Data/PhD work/Chile/data/test/sal/2017-4/';
% B = '/Volumes/My Data/PhD work/Chile/data/test/tem/2017-4/';
A = "H:\PhD work\bp_cross\Chile\data\salinity\2014-2015";
B = "H:\PhD work\bp_cross\Chile\data\temperature\2014-2015";



addpath("H:\PhD work\Miguel\seawater_ver3_3.1")
Lat = ncread('H:\PhD work\bp_cross\Chile\data\salinity\2014-2015\1.nc','lat');
Lon = ncread('H:\PhD work\bp_cross\Chile\data\salinity\2014-2015\1.nc','lon');
Depth = ncread('H:\PhD work\bp_cross\Chile\data\salinity\2014-2015\1.nc','depth');
%ssh = ncread('/Volumes/My Data/PhD work/Chile/data/test/ssh/2017-4/1.nc','surf_el');
ssh = ncread("H:\PhD work\bp_cross\Chile\data\ssh\2014-2015\1.nc",'surf_el');
mon = 12;

length(ssh(1,1,:))
m = size(Lon,1); 
n = size(Lat,1);
bp1 = nan(size(ssh)); 
Temp = nan(size(ssh));
bpssh = nan(size(ssh));
bpga = nan(size(ssh)); % max(l1)*12 gives a rough estimate for the total number of time steps


rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
P = rho*Depth.*g/10^4; %hpa 
num = 1;
l1 = zeros(mon,1);
l2 = zeros(mon,1);
%%
all_ga = NaN(m, n, 40, length(ssh(1,1,:)));

for mo = 1:mon
    %for mo = 1:mon;
    %     ncdisp([A(2).folder '/' A(mo).name])
    %     ncdisp([A(2).folder '/' B(mo).name])
    fileName = [num2str(mo) '.nc'];
    S = ncread(fullfile(A, fileName),'salinity');
    T = ncread(fullfile(B, fileName), 'water_temp');
    l1(mo)=length(S(1,1,1,:));
    l2(mo)=length(T(1,1,1,:));
    for yr = 1:length(S(1,1,1,:)) % every month have different time steps
        for x = 1:m
            for y = 1:n %every point in the region(x,y)
                if ~isnan(S(x,y,end,yr))
                    [ga] = sw_gpan(squeeze(S(x,y,1:end,yr)),...
                        squeeze(T(x,y,1:end,yr)),P(1:end));
                    bp1(x,y,num) = rho*g*(ssh(x,y,num)-ga(end)/g);
                    Temp(x,y,num) = squeeze(T(x,y,end,yr));
                    bpga(x,y,num) = rho*g*ga(end)/g; 
                    all_ga(x, y, 1:end, num) = rho*g*ga(1:end)/g;
                 

                else
                    for depth = 2:length(S(1,1,:,1))
                        if isnan(S(x,y,depth,yr)) && ~isnan(S(x,y,depth-1,yr))
                            % this layer is nan, last layer is valuable.
                            [ga] = sw_gpan(squeeze(S(x,y,1:(depth-1),yr)),...
                                squeeze(T(x,y,1:(depth-1),yr)),P(1:(depth-1)));
                            bp1(x,y,num) = rho*g*(ssh(x,y,num)-ga(end)/g);
                            Temp(x,y,num) = squeeze(T(x,y,(depth-1),yr));
                            bpga(x,y,num) = rho*g*ga(end)/g;
                            all_ga(x, y, 1:depth-1, num) = rho*g*ga(1:depth-1)/g;
                        

                        end
                    end
                end
            end
        end
        num = num+1;
    end
end

sum(l1)
sum(l2)
return
% remove mean value of each year
%%
% Initialize containers
S_all = [];
T_all = [];

% Loop through each month
for mo = 1:mon
    fileName = [num2str(mo) '.nc'];
    
    % Read salinity and temperature from NetCDF
    S = ncread(fullfile(A, fileName), 'salinity');      % [m x n x 40 x t1]
    T = ncread(fullfile(B, fileName), 'water_temp');    % [m x n x 40 x t2]
    
    % Sanity check: make sure S and T are the same shape
    if ~isequal(size(S), size(T))
        warning('Shape mismatch in month %d. Skipping...', mo);
        continue;
    end

    % Concatenate along the 4th (time) dimension
    S_all = cat(4, S_all, S);  % Size: [m x n x 40 x total_time]
    T_all = cat(4, T_all, T);
    
    fprintf('Month %d loaded. Total time steps so far: %d\n', mo, size(S_all, 4));
end
%%
% Given dimensions
[m, n, old_z] = size(bpga);  % m = 126, n = 76, old_z = 228
new_z = 248;  % Target third dimension

% Initialize the output array
bpga_int = nan(m, n, new_z);

% Define the old and new grid for the third dimension
old_indices = linspace(1, old_z, old_z);
new_indices = linspace(1, old_z, new_z);

% Loop through each element in the 2D grid
for x = 1:m
    for y = 1:n
        % Extract the data along the third dimension for this (x,y)
        data_slice = squeeze(bpga(x, y, :));

        % Check for non-NaN data to avoid errors in interpolation
        if any(~isnan(data_slice))
            % Interpolate data to new grid
            bpga_int(x, y, :) = interp1(old_indices, data_slice, new_indices, 'linear', 'extrap');
        end
    end
end

%% Now bpga_int is the interpolated version of bpga with size 126x76x248
bpssh = rho*g*ssh; 
% Given dimensions
[m, n, old_z] = size(bpssh);  % m = 126, n = 76, old_z = 228
new_z = 248;  % Target third dimension

% Initialize the output array
bpssh_int = nan(m, n, new_z);

% Define the old and new grid for the third dimension
old_indices = linspace(1, old_z, old_z);
new_indices = linspace(1, old_z, new_z);

% Loop through each element in the 2D grid
for x = 1:m
    for y = 1:n
        % Extract the data along the third dimension for this (x,y)
        data_slice = squeeze(bpssh(x, y, :));

        % Check for non-NaN data to avoid errors in interpolation
        if any(~isnan(data_slice))
            % Interpolate data to new grid
            bpssh_int(x, y, :) = interp1(old_indices, data_slice, new_indices, 'linear', 'extrap');
        end
    end
end
%%
% Loop through each element in the 2D grid
bp1_int = nan(m, n, new_z);
for yr = 1:length(S(1,1,1,:))
    for x = 1:m
        for y = 1:n
            bp1_int(x,y,num) = bpssh_int(x,y,num) - bpga_int(x,y,num);
        end
    end
    num = num +1;
end

%%
%filepath = '/Volumes/My Data/PhD work/bp_cross/Chile/data/salinity/2016/1.nc';

% filepath = '/Volumes/My Data/PhD work/Chile/data/test/ssh/2017-4/1.nc';
filepath = "H:\PhD work\bp_cross\Chile\data\ssh\2014-2015\1.nc";
timeData = ncread(filepath, 'time');
referenceDate = datetime('2000-01-01 00:00:00');
Timedate = referenceDate +hours(timeData);
%%
sta = 7941;
sta2 = 7402;
sta3 = 6770;

bp2 = nan(size(bp1)); 
for x=1:m
    for y=1:n
        bp2(x,y,:)=bp1(x,y,1:end)-nanmean(bp1(x,y,1:end));
    end
end
%%
bpssh = rho*g*ssh; 
bpssh2 = nan(size(bpssh)); 
for x=1:m
    for y=1:n
        bpssh2(x,y,:)=bpssh(x,y,1:end)-nanmean(bpssh(x,y,1:end));
    end
end

bpga2 = nan(size(bpga)); 
for x=1:m
    for y=1:n
        bpga2(x,y,:)=bpga(x,y,1:end)-nanmean(bpga(x,y,1:end));
    end
end

all_ga2 = nan(size(all_ga)); 
for x=1:m
    for y=1:n
        for num = 1:40
        all_ga2(x,y,num,:)=all_ga(x,y,num,1:end)-nanmean(all_ga(x,y,num,1:end));
        end
    end
end
%%
profile(isnan(profile)) = 0;


%%
% Choose your specific location (e.g., lat=50, lon=30)
ilat = 51;% -30
ilon = 101; %-72

% Extract the [depth x time] profile
profile = squeeze(all_ga2(ilon, ilat, :, :))/100;  % Now size is [40 x 2854]
ssh_pressure = squeeze(bpssh2(ilon,ilat,:)/100);
% Ensure ssh_pressure is a row vector [1 x 2854]
ssh_pressure_row = ssh_pressure(:)';  % Transpose to row if needed

% Replicate along the depth dimension to match profile shape
ssh_matrix = repmat(ssh_pressure_row, size(profile, 1), 1);  % [40 x 2854]

% Subtract: (ssh at surface - profile at each depth)
anomaly = ssh_matrix - profile;  % [40 x 2854]

% Optional: actual depth vector if available
depths = 1:40;  % Replace with actual depth values if you have them

%% Plotting
figure;
subplot(311)
imagesc(Timedate, Depth, profile);  % time on x, depth on y
set(gca, 'YDir', 'reverse');         % depth increases downward
xlabel('Time');
ylabel('Depth');
title('Cumulative Steric Pressure from Surface to Bottom');
% cb = colorbar;
%cb.Label.String = 'Pressure Anomaly (hPa)';  % Optional: add a label
% Place colorbar below the plot


% if Timedate is datetime, otherwise format manually
subplot(312)
imagesc(Timedate, Depth, anomaly);  % Replace 1:40 with actual depth values if available
set(gca, 'YDir', 'reverse');  % So deeper depths appear lower
xlabel('Time');
ylabel('Depth');
title('Bottom Pressure Profile');
% Place colorbar below the plot
% cb = colorbar;
%cb.Label.String = 'Pressure Anomaly (hPa)';  % Optional: add a label

datetick('x', 'keeplimits');

subplot(313)
plot(Timedate,squeeze(bpga2(ilon,ilat,:))/100)
hold on
plot(Timedate,squeeze(bpssh2(ilon,ilat,:))/100)
plot(Timedate,squeeze(bp2(ilon,ilat,:))/100)
legend('Steric pressure','SSH pressure','All pressure','Location','best')