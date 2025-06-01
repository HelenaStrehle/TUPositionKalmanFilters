function [acc_offset, gyro_offset, acc_noise, gyro_noise] = calibrate_imu(filename)

    fid = fopen(filename, 'r');
    accDataKal = textscan(fid, '%s', 'Delimiter', '\n');
    accDataKal = accDataKal{1};
    fclose(fid);

    % Retsæt '0025' til '2025'
    accDataKal = regexprep(accDataKal, '0025-', '2025-');

    imuPatternKal = '(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) .*?Accel: X=(-?\d+) Y=(-?\d+) Z=(-?\d+) Gyro: X=(-?\d+) Y=(-?\d+) Z=(-?\d+)';

    imuDataKal = {};

    for i = 1:length(accDataKal)
        tokens = regexp(accDataKal{i}, imuPatternKal, 'tokens');
        if ~isempty(tokens)
            t = tokens{1};
            if numel(t) == 7
                imuDataKal(end+1, :) = { ...
                    datetime(t{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC'), ...
                    str2double(t{2}), ...
                    str2double(t{3}), ...
                    str2double(t{4}), ...
                    str2double(t{5}), ...
                    str2double(t{6}), ...
                    str2double(t{7})};
            end
        end
    end

    % Hvis der ikke er nogen data, vis fejl
    if isempty(imuDataKal)
        error('Ingen gyldige IMU-data fundet i logfilen.');
    end


    imuTableKal = cell2table(imuDataKal, 'VariableNames', ...
        {'Timestamp', 'AccelX', 'AccelY', 'AccelZ', 'GyroX', 'GyroY', 'GyroZ'});

    imuTableKal.AccelX = imuTableKal.AccelX * 9.81 / 1000;
    imuTableKal.AccelY = imuTableKal.AccelY * 9.81 / 1000;
    imuTableKal.AccelZ = imuTableKal.AccelZ * 9.81 / 1000;
    imuTableKal.GyroX = deg2rad(imuTableKal.GyroX / 1000);
    imuTableKal.GyroY = deg2rad(imuTableKal.GyroY / 1000);
    imuTableKal.GyroZ = deg2rad(imuTableKal.GyroZ / 1000);

    % Forventet acceleration
    expected = [-9.81, 0, 0];

    % Beregn offset
    measured = [mean(imuTableKal.AccelX), mean(imuTableKal.AccelY), mean(imuTableKal.AccelZ)];
    acc_offset = measured - expected;

    gyro_offset = [ ...
        mean(imuTableKal.GyroX), ...
        mean(imuTableKal.GyroY), ...
        mean(imuTableKal.GyroZ)];

    % Beregn støj
    acc_noise = [ ...
        std(imuTableKal.AccelX - measured(1)), ...
        std(imuTableKal.AccelY - measured(2)), ...
        std(imuTableKal.AccelZ - measured(3))];

    gyro_noise = [ ...
        std(imuTableKal.GyroX - gyro_offset(1)), ...
        std(imuTableKal.GyroY - gyro_offset(2)), ...
        std(imuTableKal.GyroZ - gyro_offset(3))];

    % Udskriv resultater
    fprintf('Accelerometer offset (m/s²):\nX: %.4f  Y: %.4f  Z: %.4f\n', acc_offset);
    fprintf('Accelerometer noise σ (m/s²):\nX: %.4f  Y: %.4f  Z: %.4f\n', acc_noise);
    fprintf('Gyroscope offset (rad/s):\nX: %.6f  Y: %.6f  Z: %.6f\n', gyro_offset);
    fprintf('Gyroscope noise σ (rad/s):\nX: %.6f  Y: %.6f  Z: %.6f\n', gyro_noise);
end