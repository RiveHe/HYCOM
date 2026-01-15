sta = 7941;
sta2 = 7402;
sta3 = 6770;
load('/Volumes/My Data/PhD work/bp_cross/Chile/result/2016-newresult.mat');
Lat = ncread('/Volumes/My Data/PhD work/bp_cross/Chile/data/salinity/2016/1.nc','lat');
Lon = ncread('/Volumes/My Data/PhD work/bp_cross/Chile/data/salinity/2016/1.nc','lon');
Lat = Lat';
desiredSize = [126, 76]; %lon,lat

% Replicate the columns to create Lon with the desired size
Lon = repmat(Lon, 1, desiredSize(2));
Lat = repmat(Lat, desiredSize(1), 1);
yr =1;
addpath('/Users/river/Downloads/cmocean-main')
cmap = cmocean('delta');

%%

A = '/Volumes/My Data/PhD work/bp_cross/Chile/data/vel/2016';
u = [];
v = [];

for mo = 1:10
    fileName = [num2str(mo) '.nc'];
    %ncdisp([A(2).folder '/' A(mo).name])
    udata = ncread(fullfile(A, fileName),'water_u_bottom');
    vdata = ncread(fullfile(A, fileName),'water_v_bottom');
    u = cat(3, u, udata);
    v = cat(3, v, vdata);
end

eff_num = sum(sum(~isnan(u(:,:,yr))));
eff_u2 = zeros(eff_num,size(u,3));
eff_v2 = zeros(eff_num,size(v,3));
num=1;
for is1 = 1:size(Lon,1)
    for is2 = 1:size(Lon,2)
        if ~isnan(u(is1,is2,yr))
            plot(Lon(is1,is2),Lat(is1,is2),'ro','MarkerFaceColor','r')
            text(Lon(is1,is2),Lat(is1,is2),num2str(num));

            eff_u2(num,:) = squeeze(u(is1,is2,:))*10^2;
            eff_v2(num,:) = squeeze(v(is1,is2,:))*10^2;

            num=num+1;
        end
    end
end
%%
% Desired size: same number of rows, more columns
new_rows = 8096;
new_cols = 2406;

% Resizing the matrix using bicubic interpolation
eff_u2 = imresize(eff_u2, [new_rows new_cols], 'bicubic');
eff_v2 = imresize(eff_v2, [new_rows new_cols], 'bicubic');
%%
% Assuming u and v are your original 3D matrices of size 126x76x2399

% Define the size of the new data
new_size = [126, 76, 2406];

% Assuming u and v are originally of size [numX, numY, numZ]
[numX, numY, numZ] = size(u);

% Create original grid
[oX, oY, oZ] = meshgrid(1:numY, 1:numX, 1:numZ);

% Define new grid size
new_size = [126, 76, 2406];

% Create new grid using linspace to properly distribute the coordinates
[X, Y, Z] = meshgrid(1:new_size(2), 1:new_size(1), linspace(1, numZ, new_size(3)));
% Interpolate u and v
u_new = interp3(oX, oY, oZ, u, X, Y, Z, 'linear');
v_new = interp3(oX, oY, oZ, v, X, Y, Z, 'linear');
%%
filepath = '/Volumes/My Data/PhD work/bp_cross/Chile/data/ssh/2016/1.nc';
timeData = ncread(filepath, 'time');
referenceDate = datetime('2000-01-01 00:00:00');
Timedate = referenceDate +hours(timeData);
time = 1:1:size(eff_u2,2);
%%
% Read bathymetry data from netCDF file
filename_nc = '/Volumes/My Data/PhD work/bp_cross/Chile/dataset/ncmap/chile2.nc';
gebconc = netcdf.open('/Volumes/My Data/PhD work/bp_cross/Chile/dataset/ncmap/chile2.nc', 'NOWRITE');

XGRID = netcdf.getVar(gebconc,0);
YGRID = netcdf.getVar(gebconc,1);
BATHY = netcdf.getVar(gebconc,2);
netcdf.close(gebconc);
%% make video with current and pressure on the right panels 
sta = 7941;
sta2 = 7402;
sta3 = 6770;
% Assuming Timedate is an array of datetime objects
% Define the specific date and time you're looking for
targetDateTime = datetime('03-10-2016 00:00:00', 'InputFormat', 'MM-d-yyyy HH:mm:ss');

% Find the index
index = find(Timedate == targetDateTime);
index = 1942;

% Display the index
disp('Index of the specific date and time:');
disp(index);

videoFilePath = '/Volumes/My Data/PhD work/bp_cross/Chile/dataset/2Dbpuv1_2016.mp4';

% Initialize VideoWriter with correct output file format and frame rate
vid0bj1 = VideoWriter('./2Dbpuv1_2016.avi');
vid0bj1.FrameRate = 10;
open(vid0bj1);

% Set up the figure with fixed dimensions
f1 = figure;
set(f1, 'Position', [100, 100, 1400, 670]);  % [left, bottom, width, height]
addpath('/Users/river/Downloads/cmocean-main')
cmap = cmocean('delta');

tic;
for yr = index
    clf;  % Clear the figure to ensure clean slate for each frame
    subplot(131);
    selectedDate = Timedate(yr);  % Assuming Timedate is a function you've defined earlier

    % Your plotting code
    pcolor(Lon, Lat, (bp2(:,:,yr)/10^2));
    set(gca, 'YDir', 'normal');
    shading interp;
    caxis([-10 10]);
    hold on;
    colormap(cmap);
    colorbar('southoutside');
    quiver(Lon,Lat,u(:,:,yr),v(:,:,yr),'Color','r','AutoScaleFactor',4,...
        'LineWidth',1.0,'MarkerFaceColor','r');
    quiver(-71.5,-33.4,0.05,0,'Color','r','AutoScaleFactor',4,...
        'LineWidth',1.0,'MarkerFaceColor','r');%plot scaler
    text(-71.5,-33.5,'5cm/s','FontSize',13)

    % Plot bathymetry
    contour_levels = [-500, -1000:-1000:-6000];
    [C, h] = contour(YGRID, XGRID, BATHY', contour_levels);
    clabel(C, h, contour_levels, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 1000);
    colormap(cmap);
    plot(eff_sta(sta2,1),eff_sta(sta2,2),'o','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor','k')%plot staion position
    text(eff_sta(sta2,1)+0.05,eff_sta(sta2,2)+0.03,num2str(2),'FontSize',14)
    plot(eff_sta(sta,1),eff_sta(sta,2),'o','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor','k')%plot staion position
    text(eff_sta(sta,1)+0.05,eff_sta(sta,2)+0.03,num2str(1),'FontSize',14)
    plot(eff_sta(sta3,1),eff_sta(sta3,2),'o','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor','k')%plot staion position
    text(eff_sta(sta3,1)+0.05,eff_sta(sta3,2)+0.03,num2str(3),'FontSize',14)
    axis([-73.5 -71 -34 -28]);
    title(['Bottom pressure anomaly on ', datestr(selectedDate, 'yyyy-mm-dd HH:MM:SS')]);

    subplot(3,9,5:9)
    plot(time,eff_bpt(sta,:),'b','LineWidth',1)
    hold on

    for i=1:1:size(eff_v2,2)
        p = plot([i,i+eff_u2(sta,i)],[0,eff_v2(sta,i)],'r','LineWidth',1);%plot current at station 73
        set(p,'LineStyle','-','Marker','none');
        hold on
    end
    plot(yr, eff_bpt(sta, yr), 'ko', 'MarkerSize', 10, 'LineWidth', 2);

    changeMonthIdx = find(diff([0; month(Timedate)]) ~= 0);
    firstOfMonthIdx = find(day(Timedate(changeMonthIdx)) == 1);
    firstOfMonth = changeMonthIdx(firstOfMonthIdx);
    firstOfMonthLabels = datestr(Timedate(firstOfMonth), 'dd-mmm');
    hold off
    ax = gca; % Get current axis
    ax.FontSize = 16;
    ax.XTick = firstOfMonth; % Set ticks at the first of each month
    ax.XTickLabel = firstOfMonthLabels; % Set labels for these ticks
    ax.XTickLabelRotation = 0;
    xlim([1, length(Timedate)]);
    l1 = legend('Pressure','Current');
    l1.FontSize=14;l1.Orientation='horizontal';
    l1.Location = 'southeast';
    title('Pressure and current from 2015-2016 at station 3');


    subplot(3,9,14:18)
    plot(time,eff_bpt(sta2,:),'b','LineWidth',1)
    hold on

    for i=1:1:size(eff_v2,2)
        p = plot([i,i+eff_u2(sta2,i)],[0,eff_v2(sta2,i)],'r','LineWidth',1);%plot current at station 73
        set(p,'LineStyle','-','Marker','none');
        hold on
    end
    plot(yr, eff_bpt(sta2, yr), 'ko', 'MarkerSize', 10, 'LineWidth', 2);

    changeMonthIdx = find(diff([0; month(Timedate)]) ~= 0);
    firstOfMonthIdx = find(day(Timedate(changeMonthIdx)) == 1);
    firstOfMonth = changeMonthIdx(firstOfMonthIdx);
    firstOfMonthLabels = datestr(Timedate(firstOfMonth), 'dd-mmm');
    hold off
    ax = gca; % Get current axis
    ax.FontSize = 16;
    ax.XTick = firstOfMonth; % Set ticks at the first of each month
    ax.XTickLabel = firstOfMonthLabels; % Set labels for these ticks
    ax.XTickLabelRotation = 0;
    xlim([1, length(Timedate)]);
    l1 = legend('Pressure','Current');
    l1.FontSize=14;l1.Orientation='horizontal';
    l1.Location = 'southeast';
    title('Pressure and current from 2015-2016 at station 3');


    subplot(3,9,23:27)

    % Plotting eff_bpt over time for station 3
    plot(time,eff_bpt(sta3,:),'b','LineWidth',1)
    hold on

    for i=1:1:size(eff_v2,2)
        p = plot([i,i+eff_u2(sta3,i)],[0,eff_v2(sta3,i)],'r','LineWidth',1);%plot current at station 73
        set(p,'LineStyle','-','Marker','none');
        hold on
    end
    plot(yr, eff_bpt(sta3, yr), 'ko', 'MarkerSize', 10, 'LineWidth', 2);

    changeMonthIdx = find(diff([0; month(Timedate)]) ~= 0);
    firstOfMonthIdx = find(day(Timedate(changeMonthIdx)) == 1);
    firstOfMonth = changeMonthIdx(firstOfMonthIdx);
    firstOfMonthLabels = datestr(Timedate(firstOfMonth), 'dd-mmm');
    hold off
    ax = gca; % Get current axis
    ax.FontSize = 16;
    ax.XTick = firstOfMonth; % Set ticks at the first of each month
    ax.XTickLabel = firstOfMonthLabels; % Set labels for these ticks
    ax.XTickLabelRotation = 0;
    xlim([1, length(Timedate)]);
    l1 = legend('Pressure','Current');
    l1.FontSize=14;l1.Orientation='horizontal';
    l1.Location = 'southeast';
    title('Pressure and current from 2015-2016 at station 3');


    % Capture the frame
    currFrame = getframe(f1);
    writeVideo(vid0bj1, currFrame);
    toc;
end

% Close the video file to finalize
close(vid0bj1);

% Display completion message
disp(['Video saved to ' videoFilePath]);

%% Only bottom pressure on the right panels
targetDateTime = datetime('01-1-2016 00:00:00', 'InputFormat', 'MM-d-yyyy HH:mm:ss');

% Find the index
index = find(Timedate == targetDateTime);

% Display the index
disp('Index of the specific date and time:');
disp(index);

videoFilePath = '/Users/river/Downloads/chilevideo/';

% Initialize VideoWriter with correct output file format and frame rate
videoFileName = '2Dbpuv1_2016-2.avi';
fullVideoPath = fullfile(videoFilePath, videoFileName);

vid0bj1 = VideoWriter(fullVideoPath);
vid0bj1.FrameRate = 10;
open(vid0bj1);

% Set up the figure with fixed dimensions
f1 = figure;
set(f1, 'Position', [100, 100, 1440, 840]);  % [left, bottom, width, height]

tic;
for yr = index:index2
    clf;  % Clear the figure to ensure clean slate for each frame
    selectedDate = Timedate(yr);  % Assuming Timedate is a function you've defined earlier
    subplot(131);
    pcolor(Lon, Lat, (bp2(:,:,yr)/10^2));
    set(gca, 'YDir', 'normal');
    shading interp;
    caxis([-10 10]);
    hold on;
    colormap(cmap);
    % Create colorbar at the specified location
    cb = colorbar('southoutside');

    % Add label to the colorbar
    xlabel(cb, 'Pressure Anomaly (hPa)');
    % Your plotting code
    quiver(Lon,Lat,u(:,:,yr),v(:,:,yr),'Color','r','AutoScaleFactor',4,...
        'LineWidth',1.0,'MarkerFaceColor','r');
    quiver(-71.5,-33.4,0.05,0,'Color','r','AutoScaleFactor',4,...
        'LineWidth',1.0,'MarkerFaceColor','r');%plot scaler
    text(-71.5,-33.5,'5cm/s','FontSize',13)

    % Plot bathymetry
    contour_levels = [-500, -1000:-1000:-6000];
    [C, h] = contour(YGRID, XGRID, BATHY', contour_levels);
    clabel(C, h, contour_levels, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 1000);
    colormap(cmap);
    plot(eff_sta(sta2,1),eff_sta(sta2,2),'o','MarkerFaceColor','r','MarkerSize',8,'MarkerEdgeColor','k')%plot staion position
    text(eff_sta(sta2,1)+0.05,eff_sta(sta2,2)+0.03,num2str(2),'FontSize',18,'FontWeight', 'bold')
    plot(eff_sta(sta,1),eff_sta(sta,2),'o','MarkerFaceColor','r','MarkerSize',8,'MarkerEdgeColor','k')%plot staion position
    text(eff_sta(sta,1)+0.05,eff_sta(sta,2)+0.03,num2str(1),'FontSize',18,'FontWeight', 'bold')
    plot(eff_sta(sta3,1),eff_sta(sta3,2),'o','MarkerFaceColor','r','MarkerSize',8,'MarkerEdgeColor','k')%plot staion position
    text(eff_sta(sta3,1)+0.05,eff_sta(sta3,2)+0.03,num2str(3),'FontSize',18,'FontWeight', 'bold')
    axis([-73.5 -71 -34 -28]);
    ax = gca; % Get current axis
    ax.FontSize = 16;
    title(['BP Anomaly and UV on ', datestr(selectedDate, 'yyyy-mm-dd HH:MM:SS')]);

    subplot(3,9,5:9)
    rmsValue = sqrt(nanmean(eff_bpt(sta,:).^2));
    plot(time,eff_bpt(sta,:),'b','LineWidth',1)
    hold on
   
    plot(yr, eff_bpt(sta, yr), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    
    % Plot horizontal line at y=0 for the range of x values specified by yr
    yline(0, 'k--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    yline(-rmsValue, 'r--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    yline(rmsValue, 'r--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    changeMonthIdx = find(diff([0; month(Timedate)]) ~= 0);
    firstOfMonthIdx = find(day(Timedate(changeMonthIdx)) == 1);
    firstOfMonth = changeMonthIdx(firstOfMonthIdx);
    firstOfMonthLabels = datestr(Timedate(firstOfMonth), 'dd-mmm');
   text(max(xlim)-150, max(ylim)-0.5, ['RMS: ', num2str(rmsValue, '%.2f'), ' hPa'], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'color','b','FontSize', 14, 'FontWeight', 'bold');
  
    hold off
    ax = gca; % Get current axis
    ax.FontSize = 16;
    ax.XTick = firstOfMonth; % Set ticks at the first of each month
    ax.XTickLabel = firstOfMonthLabels; % Set labels for these ticks
    ax.XTickLabelRotation = 0;
    xlim([1, length(Timedate)]);
    xlabel('Time');
    ylabel('Pressure Anomaly (hPa)')
    title('Pressure from 2015-2016 at station 1');

    subplot(3,9,14:18)
    rmsValue = sqrt(nanmean(eff_bpt(sta2,:).^2));
    plot(time,eff_bpt(sta2,:),'b','LineWidth',1)
    hold on
   
    plot(yr, eff_bpt(sta2, yr), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    
    % Plot horizontal line at y=0 for the range of x values specified by yr
    yline(0, 'k--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    yline(-rmsValue, 'r--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    yline(rmsValue, 'r--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    changeMonthIdx = find(diff([0; month(Timedate)]) ~= 0);
    firstOfMonthIdx = find(day(Timedate(changeMonthIdx)) == 1);
    firstOfMonth = changeMonthIdx(firstOfMonthIdx);
    firstOfMonthLabels = datestr(Timedate(firstOfMonth), 'dd-mmm');
    text(max(xlim)-150, max(ylim)-0.5, ['RMS: ', num2str(rmsValue, '%.2f'), ' hPa'], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'color','b','FontSize', 14, 'FontWeight', 'bold');
  
    hold off
    ax = gca; % Get current axis
    ax.FontSize = 16;
    ax.XTick = firstOfMonth; % Set ticks at the first of each month
    ax.XTickLabel = firstOfMonthLabels; % Set labels for these ticks
    ax.XTickLabelRotation = 0;
    xlim([1, length(Timedate)]);
    xlabel('Time');
    ylabel('Pressure Anomaly (hPa)')
    title('Pressure from 2015-2016 at station 2');

    subplot(3,9,23:27)
    rmsValue = sqrt(nanmean(eff_bpt(sta3,:).^2));

    plot(time,eff_bpt(sta3,:),'b','LineWidth',1)
    hold on
   
    plot(yr, eff_bpt(sta3, yr), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    
    % Plot horizontal line at y=0 for the range of x values specified by yr
    yline(0, 'k--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    yline(-rmsValue, 'r--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    yline(rmsValue, 'r--', 'LineWidth', 1); % Plotting a dashed line ('k--')
    changeMonthIdx = find(diff([0; month(Timedate)]) ~= 0);
    firstOfMonthIdx = find(day(Timedate(changeMonthIdx)) == 1);
    firstOfMonth = changeMonthIdx(firstOfMonthIdx);
    firstOfMonthLabels = datestr(Timedate(firstOfMonth), 'dd-mmm');
    text(max(xlim)-150, max(ylim)-0.5, ['RMS: ', num2str(rmsValue, '%.2f'), ' hPa'], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'color','b','FontSize', 14, 'FontWeight', 'bold');
  
    hold off
    ax = gca; % Get current axis
    ax.FontSize = 16;
    ax.XTick = firstOfMonth; % Set ticks at the first of each month
    ax.XTickLabel = firstOfMonthLabels; % Set labels for these ticks
    ax.XTickLabelRotation = 0;
    xlim([1, length(Timedate)]);
    xlabel('Time');
    ylabel('Pressure Anomaly (hPa)')
    title('Pressure from 2015-2016 at station 3');
    
    % Capture the frame
    currFrame = getframe(f1);
    writeVideo(vid0bj1, currFrame);
    toc;
end

% Close the video file to finalize
close(vid0bj1);

% Display completion message
disp(['Video saved to ' videoFilePath]);