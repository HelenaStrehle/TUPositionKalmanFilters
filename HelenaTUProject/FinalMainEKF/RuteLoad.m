function [imuTable, gnssTable, bleTable] = RuteLoad(filename)
%% Read file and clean data
fid = fopen(filename, 'r');
logData = textscan(fid, '%s', 'Delimiter', '\n');
logData = logData{1};
fclose(fid);

% correct '0025' til '2025'
logData = regexprep(logData, '0025-', '2025-');

% dropped message
drop_seq = 'uart:~\$ uart:~\$ uart:~\$';

gnss_pattern = 'GNSS Date & Time'; % Find gnss data
ble_pattern = 'MAC: ([0-9A-F:]+) RSSI: (-?\d+) dBm'; %Find ble data

outputLog = {};
latestTimestamp = '';
currentSecond = '';
ms_counter = -10;

for i = 1:length(logData)
    currentLine = logData{i};
    timestamp = regexp(currentLine, '(2025-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})', 'match', 'once');
      isBLE = ~isempty(regexp(currentLine, ble_pattern, 'once')); 

    if ~isempty(timestamp)
        latestTimestamp = timestamp;
        if ~strcmp(timestamp, currentSecond)
            currentSecond = timestamp;
            ms_counter = 0;
        elseif ~isBLE % samling ble
            ms_counter = ms_counter + 10;
        end
    else
        if ~isBLE % sampling ble
            ms_counter = ms_counter + 10;
        end
    end

    if ~isempty(regexp(currentLine, drop_seq, 'once'))
        outputLog{end+1,1} = [latestTimestamp, sprintf('.%03d', ms_counter), ' message dropped'];
        ms_counter = ms_counter + 10;
    end

    lineWithoutTimestamp = regexprep(currentLine, '^.*?(\[.*?UTC\])?\s*', '');

    if ~isempty(regexp(currentLine, gnss_pattern, 'once'))
        ms_counter = ms_counter + 10;
    end

    outputLog{end+1,1} = sprintf('%s.%03d %s', latestTimestamp, ms_counter, lineWithoutTimestamp);
end

fid = fopen('cleaned_log_output.log', 'w');
fprintf(fid, '%s\n', outputLog{:});
fclose(fid);

%Split data
imuPattern = '(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}).*Accel: X=(-?\d+) Y=(-?\d+) Z=(-?\d+) Gyro: X=(-?\d+) Y=(-?\d+) Z=(-?\d+)';
blePattern = '(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}).*MAC: ([0-9A-F:]+) RSSI: (-?\d+) dBm';
gnssPattern = '(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}).*Lat: ([0-9.E+-]+), Lon: ([0-9.E+-]+)';

imuData = [];
gnssData = [];
bleData = [];

for i = 1:length(outputLog)
    currentLine = outputLog{i};
    imuTokens = regexp(currentLine, imuPattern, 'tokens');
    bleTokens = regexp(currentLine, blePattern, 'tokens');
    gnssTokens = regexp(currentLine, gnssPattern, 'tokens');

    if ~isempty(imuTokens)
        imuTokens = imuTokens{1};
       
        imuData = [imuData; {datetime(imuTokens{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'TimeZone', 'UTC'), ...
                    str2double(imuTokens{2}), str2double(imuTokens{3}), str2double(imuTokens{4}), ...
                    str2double(imuTokens{5}), str2double(imuTokens{6}), str2double(imuTokens{7})}];
    elseif ~isempty(bleTokens)
        bleTokens = bleTokens{1};
        bleData = [bleData; {
            datetime(bleTokens{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'TimeZone', 'UTC'), ...
            bleTokens{2}, str2double(bleTokens{3})}];

    elseif ~isempty(gnssTokens)
        gnssTokens = gnssTokens{1};
        timestamp = datetime(gnssTokens{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'TimeZone', 'UTC');
        lat = str2double(gnssTokens{2}) / 1e9;
        lon = str2double(gnssTokens{3}) / 1e9;
        gnssData = [gnssData; {timestamp, lat, lon}];
    end
end

% if size(gnssData,1) > 2
%     gnssData(1,:) = [];  % remove first unreliable gnss
% 
% end


imuTable = cell2table(imuData, 'VariableNames', {'Timestamp', 'AccelX', 'AccelY', 'AccelZ', 'GyroX', 'GyroY', 'GyroZ'});

if ~isempty(gnssData)
    gnssTable = cell2table(gnssData, 'VariableNames', {'Timestamp', 'Latitude', 'Longitude'});
else
    gnssTable = table();
end

if ~isempty(bleData)
    bleTable = cell2table(bleData, 'VariableNames', {'Timestamp', 'MAC', 'RSSI'});
else
    bleTable = table();
end

% Time in ms
imuTable.Timestamp.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
if ~isempty(gnssTable), gnssTable.Timestamp.Format = 'yyyy-MM-dd HH:mm:ss.SSS'; end
if ~isempty(bleTable), bleTable.Timestamp.Format = 'yyyy-MM-dd HH:mm:ss.SSS'; end

% Time correction
t = imuTable.Timestamp;
dt = seconds(diff(t));
valid_rows = [true; dt > seconds(0)];
imuTable = imuTable(valid_rows, :);
end
