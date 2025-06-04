clear all;
%% Calibrating IMU
[acc_offset, gyro_offset, acc_noise, gyro_noise] = calibrate_imu('StartNorthNotMove.log');
%% Load data from TU700

[imuTable, gnssTable, bleTable] = RuteLoad('TrolleyTest.log');
gnssTable(1,:) = [];  % remove first unreliable gnss

%% Load Geo tracker ref and edit time

% Ref fix
RuteRef = GeotrackerRef('TrolleyTestRef.kml');
[yGNSSRef, xGNSSRef, utmzone] = deg2utm(RuteRef.Latitude, RuteRef.Longitude);
RuteRef.TimeUTC = datetime(RuteRef.TimeUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
RuteRef.TimeUTC.TimeZone = '';

% Time zone settings
imuTable.Timestamp.TimeZone = '';
gnssTable.Timestamp.TimeZone = '';
bleTable.Timestamp.TimeZone = '';

t_start = min(gnssTable.Timestamp);
t_stop = max(RuteRef.TimeUTC);

%%
% Trim data
imuTable(imuTable.Timestamp < t_start | imuTable.Timestamp > t_stop, :) = [];
gnssTable(gnssTable.Timestamp < t_start | gnssTable.Timestamp > t_stop, :) = [];
bleTable(bleTable.Timestamp < t_start | bleTable.Timestamp > t_stop, :) = [];

%% Data correcion

% Remove offset
[accX_corr, accY_corr, accZ_corr, GyroX_corr, GyroY_corr, GyroZ_corr] = preProcessSensorData(imuTable, acc_offset, gyro_offset);

%gnss data to utm
[yGNSS, xGNSS, utmzone] = deg2utm(gnssTable.Latitude, gnssTable.Longitude);
%% Global

N = height(imuTable);

t = imuTable.Timestamp;
dt_all = seconds(diff(t));

Rpi2 = [0 -1; 1 0];

theta0 = deg2rad(0);
h0 = [cos(theta0); sin(theta0)];
h0 = h0 / norm(h0);  % Initiate heading vector

% State transition
f = @(x, u, dt) [
    x(1) + x(3)*x(4)*dt;  % x + v*h_x*dt
    x(2) + x(3)*x(5)*dt;  % y + v*h_y*dt
    x(3) + u(1)*dt;       % v + acc*dt
    x(4) + u(2)*(-x(5))*dt; % h_x + ω*(-h_y)*dt
    x(5) + u(2)*x(4)*dt;     % h_y + ω*h_x*dt

    ];

% Jacobian Matrix
computeF = @(x, u, dt) [
    1, 0, x(4)*dt, x(3)*dt, 0;
    0, 1, x(5)*dt, 0, x(3)*dt;
    0, 0, 1, 0, 0;
    0, 0, 0, 1, -u(2)*dt;
    0, 0, 0, u(2)*dt, 1
    ];

%State Transition bias
f_bias = @(x, u, dt) [
    x(1) + x(3)*x(4)*dt;
    x(2) + x(3)*x(5)*dt;
    x(3) + (u(1) - x(6)) * dt;
    x(4) + u(2)*(-x(5))*dt;
    x(5) + u(2)*x(4)*dt;
    x(6);
    %
    ];

% Jacobian Bias
computeF_bias = @(x, u, dt) [
    1, 0, x(4)*dt, x(3)*dt, 0, 0;
    0, 1, x(5)*dt, 0, x(3)*dt, 0;
    0, 0, 1, 0, 0, -dt;
    0, 0, 0, 1, -u(2)*dt, 0;
    0, 0, 0, u(2)*dt, 1, 0;
    0, 0, 0, 0, 0, 1,

    ];

%% Dead reckoning IMU

%Initialise
x_ekf = zeros(5, N);
x_ekf(:,1) = [xGNSS(1); yGNSS(1); 0; h0];
P_pred = eye(5);
x_pred = zeros(5,N);
x_pred(:,1) = x_ekf(:,1);


% Q scaling constants
alpha_v = 1;
alpha_h = 1;

% prediction loop
for k = 2:N
    dt=dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q(x_ekf(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt);

    [x_pred, P_pred, Q] = ekf_prediction_step(x_ekf(:,k-1), P_pred, u, Q, dt, f, computeF);

    x_ekf(:,k) = x_pred;
end

[mse_IMU, rmse_IMU] = compare_estimate_to_gnss( ...
    x_ekf(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');

figure;
plot(imuTable.Timestamp, x_ekf(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('Prediction step with IMU: Velocity over time');
set(gca, 'FontSize', 18); grid on;


[x_offset, y_offset, ax] = plot_ekf_result_with_arrows( ...
    x_ekf(1,:), ...
    x_ekf(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    5, ...            % arrow_length
    10, ...              % arrow_width
    rmse_IMU ...
    );

title(ax, 'Prediction step with IMU: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);

%% GNSS first iteration

%Initialise
x_ekf_gnss = zeros(5, N);
x_ekf_gnss(:,1) = [xGNSS(1); yGNSS(1); 0; h0];
P_pred = eye(5);
x_pred = zeros(5,N);
x_pred(:,1) = x_ekf(:,1);


% Q scaling constants
alpha_v = 100000;
alpha_h = 10000;

% Update step initialising
H_gnss = [1 0 0 0 0;
    0 1 0 0 0];
R_gnss = diag([1.5 1.5])^2;

% gnss check
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

% EKF loop
for k = 2:N
    dt=dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q(x_ekf_gnss(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt);

    [x_pred, P_pred, Q] = ekf_prediction_step(x_ekf_gnss(:,k-1), P_pred, u, Q, dt, f, computeF);

    if k == 2
        gnss_used = false(height(gnssTable), 1);
    end

    for idx_gnss = 1:height(gnssTable)
        if gnss_used(idx_gnss)
            continue; 
        end

        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];
            [x_pred, P_pred, y_tilde] = ekf_gnss_correction(x_pred, P_pred, z, H_gnss, R_gnss);

            gnss_used(idx_gnss) = true;
            gnss_hits = gnss_hits + 1;
            fprintf('[k=%d] GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));

            break;
        end
    end

    x_ekf_gnss(:,k) = x_pred;

end

[mse_IMU_gnss, rmse_IMU_gnss] = compare_estimate_to_gnss( ...
    x_ekf_gnss(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');

% Velocity plot
figure;
plot(imuTable.Timestamp, x_ekf_gnss(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('EKF with GNSS Correction: Velocity over time');
set(gca, 'FontSize', 18); grid on;

%%Position plot
[x_offset, y_offset, ax] = plot_ekf_result_with_arrows( ...
    x_ekf_gnss(1,:), ...
    x_ekf_gnss(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1.5, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU_gnss ...
    );

title(ax, 'EKF with GNSS Correction: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);

%% EKF gnss with bias state


%Initialise
x_ekf_gnss_bias = zeros(6, N);
x_ekf_gnss_bias(:,1) = [xGNSS(1); yGNSS(1); 0.0; h0; 0.3];
P_pred = eye(6);
x_pred = zeros(6,N);
x_pred(:,1) = x_ekf_gnss_bias(:,1);


% Q scaling constants
alpha_v = 1000;
alpha_h = 1000;
alpha_bias = 1e-6;

% Update step initialising
H_gnss = [1 0 0 0 0 0;
    0 1 0 0 0 0];
R_gnss = diag([1.5 1.5])^2;

% gnss check
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

% EKF loop
for k = 2:N
    dt=dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q_bias = compute_Q_bias(x_ekf_gnss_bias(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    [x_pred, P_pred, Q] = ekf_prediction_step(x_ekf_gnss_bias(:,k-1), P_pred, u, Q_bias, dt, f_bias, computeF_bias);
    x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

    if k == 2
        gnss_used = false(height(gnssTable), 1);
    end

    for idx_gnss = 1:height(gnssTable)
        if gnss_used(idx_gnss)
            continue;
        end

        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];
            [x_pred, P_pred, y_tilde] = ekf_gnss_correction(x_pred, P_pred, z, H_gnss, R_gnss);
            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction
            gnss_used(idx_gnss) = true;
            gnss_hits = gnss_hits + 1;
            fprintf('[k=%d] GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));

            break;
        end
    end

    x_ekf_gnss_bias(:,k) = x_pred;

end

[mse_IMU_gnss_bias, rmse_IMU_gnss_bias] = compare_estimate_to_gnss( ...
    x_ekf_gnss_bias(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');
% Velocity plot
figure;
plot(imuTable.Timestamp, x_ekf_gnss_bias(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('EKF with GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% Velocity plot
figure;
plot(imuTable.Timestamp, x_ekf_gnss_bias(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('EKF with GNSS Correction and Bias State: Bias over time');
set(gca, 'FontSize', 18); grid on;

%%Position plot
[x_offset, y_offset, ax] = plot_ekf_result_with_arrows( ...
    x_ekf_gnss_bias(1,:), ...
    x_ekf_gnss_bias(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1.5, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU_gnss_bias ...
    );

title(ax, 'EKF with GNSS Correction and Bias State: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);
%% EKF ble with bias state

%Initialise
x_ekf_ble_bias = zeros(6, N);
x_ekf_ble_bias(:,1) = [xGNSS(1); yGNSS(1); 0.0; h0; 0.1];
P_pred = eye(6);
x_pred = zeros(6,N);
x_pred(:,1) = x_ekf_ble_bias(:,1);

% Q scaling constants
alpha_v = 1000;
alpha_h = 10;
alpha_bias = 1e-6;

% Bluetooth set up
mac_to_id = containers.Map({ ...
    '29:B9:A5:C2:3B:58', ...  % SN 10014802
    '66:97:31:C2:3B:58', ...  % SN 10014565
    'DE:52:34:C2:3B:58' ...  % SN 10014520
    }, [1, 2, 3]);

idx_valid = ismember(bleTable.MAC, keys(mac_to_id));
ble_used = bleTable(idx_valid, :);

ble_used.BeaconID = zeros(height(ble_used),1);
for i = 1:height(ble_used)
    mac = ble_used.MAC{i};
    ble_used.BeaconID(i) = mac_to_id(mac);
end

ble_used_flags = false(height(ble_used), 1);
ble_correction_indices = [];


[yGNSS_10014565, xGNSS_10014565, utmzone] = deg2utm(57.0480867, 9.9486683);
[yGNSS_10014802, xGNSS_10014802, utmzone] = deg2utm(57.0480133, 9.948805);
[yGNSS_10014520, xGNSS_10014520, utmzone] = deg2utm(57.0480467, 9.94856);

% Beacons placering
beaconPosList = [
    xGNSS_10014565, yGNSS_10014565;
    xGNSS_10014802, yGNSS_10014802;
    xGNSS_10014520, yGNSS_10014520
    ];

sigma_shadowing = 10;
R_ble = sigma_shadowing^2;

last_rssi = containers.Map('KeyType', 'double', 'ValueType', 'double');
last_time = containers.Map('KeyType', 'double', 'ValueType', 'any');
n = 1.28;
A = -77.11;

% EKF loop
for k = 2:N
    dt=dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q_bias = compute_Q_bias(x_ekf_ble_bias(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    [x_pred, P_pred, Q] = ekf_prediction_step(x_ekf_ble_bias(:,k-1), P_pred, u, Q_bias, dt, f_bias, computeF_bias);
    x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction


    for idx_ble = 1:height(ble_used)
        if ble_used_flags(idx_ble)
            continue;
        end

        deltaT_ble = seconds(ble_used.Timestamp(idx_ble) - imuTable.Timestamp(k));

        if deltaT_ble >= 0 && deltaT_ble < 1

            beacon_id = ble_used.BeaconID(idx_ble);
            rssi_now = ble_used.RSSI(idx_ble);
            time_now = ble_used.Timestamp(idx_ble);

            rssi = rssi_now;
            if isKey(last_rssi, beacon_id)
                time_diff = time_now - last_time(beacon_id);
                if time_diff <= seconds(2)
                    rssi = mean([last_rssi(beacon_id), rssi_now]);
                end
            end
            last_rssi(beacon_id) = rssi_now;
            last_time(beacon_id) = time_now;


            ref = beaconPosList(beacon_id, :);
            [x_pred, P_pred, y_ble, rssi_expected, d_measured] = ekf_ble_correction_rssi(x_pred, P_pred, ref, rssi, R_ble, n, A);

            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

            fprintf('[k=%d] BLE-korrektion #%d fra beacon #%d, ændring = %.2f dB\n, expected = %.2f m,målt = %.2f m', ...
                k, idx_ble, beacon_id, y_ble,rssi_expected, rssi);

            ble_used_flags(idx_ble) = true;
            ble_correction_indices(end+1) = k;
            break;
        end
    end

    x_ekf_ble_bias(:,k) = x_pred;

end

[mse_IMU_ble, rmse_IMU_ble] = compare_estimate_to_gnss( ...
    x_ekf_ble_bias(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');
% Vel plot
figure;
plot(imuTable.Timestamp, x_ekf_ble_bias(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('EKF with BLE-, GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% Bias plot
figure;
plot(imuTable.Timestamp, x_ekf_ble_bias(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('EKF with BLE-, GNSS Correction and Bias State: Bias over time');
set(gca, 'FontSize', 18); grid on;

%ekf plot with beacon circles
[x_offset, y_offset, ax] = plot_ekf_result_with_arrows( ...
    x_ekf_ble_bias(1,:), ...
    x_ekf_ble_bias(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1, ...            % arrow_length
    1.5, ...              % arrow_width
    rmse_IMU_ble ...
    );

title(ax, 'EKF with BLE Correction and Bias State: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);

xGNSS_10014565_shifted = xGNSS_10014565 - x_offset;
xGNSS_10014802_shifted = xGNSS_10014802 - x_offset;
xGNSS_10014520_shifted = xGNSS_10014520 - x_offset;
yGNSS_10014565_shifted = yGNSS_10014565 - y_offset;
yGNSS_10014802_shifted = yGNSS_10014802 - y_offset;
yGNSS_10014520_shifted = yGNSS_10014520 - y_offset;

hold on;
plot(xGNSS_10014565_shifted, yGNSS_10014565_shifted, 'kx', ...
    'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Bluetooth transmitter');
plot(xGNSS_10014802_shifted, yGNSS_10014802_shifted, 'kx', ...
    'LineWidth', 2, 'MarkerSize', 5, 'HandleVisibility', 'off');
plot(xGNSS_10014520_shifted, yGNSS_10014520_shifted, 'kx', ...
    'LineWidth', 2, 'MarkerSize', 5, 'HandleVisibility', 'off');

radius = 15;
theta = linspace(0, 2*pi, 100);
fill(xGNSS_10014565_shifted + radius*cos(theta), ...
    yGNSS_10014565_shifted + radius*sin(theta), ...
    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', ...
    'DisplayName', 'BLE coverage (15 m radius)');
fill(xGNSS_10014802_shifted + radius*cos(theta), ...
    yGNSS_10014802_shifted + radius*sin(theta), ...
    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(xGNSS_10014520_shifted + radius*cos(theta), ...
    yGNSS_10014520_shifted + radius*sin(theta), ...
    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

legend('Location', 'best');

%% EKF ble and gnss with bias state


%Initialise
x_ekf_all = zeros(6, N);
x_ekf_all(:,1) = [xGNSS(1); yGNSS(1); 0.01; h0; 0.3];
P_pred = eye(6);
x_pred = zeros(6,N);
x_pred(:,1) = x_ekf_all(:,1);

% Q scaling constants
alpha_v = 1000;
alpha_h = 1000;
alpha_bias = 1e-6;

% gnss update initialize
H_gnss = [1 0 0 0 0 0;
    0 1 0 0 0 0];
R_gnss = diag([4 4])^2;

% gnss check
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

% Bluetooth set up
mac_to_id = containers.Map({ ...
    '29:B9:A5:C2:3B:58', ...  % SN 10014802
    '66:97:31:C2:3B:58', ...  % SN 10014565
    'DE:52:34:C2:3B:58' ...  % SN 10014520
    }, [1, 2, 3]);

idx_valid = ismember(bleTable.MAC, keys(mac_to_id));
ble_used = bleTable(idx_valid, :);

ble_used.BeaconID = zeros(height(ble_used),1);
for i = 1:height(ble_used)
    mac = ble_used.MAC{i};
    ble_used.BeaconID(i) = mac_to_id(mac);
end

ble_used_flags = false(height(ble_used), 1);
ble_correction_indices = [];


[yGNSS_10014565, xGNSS_10014565, utmzone] = deg2utm(57.0480867, 9.9486683);
[yGNSS_10014802, xGNSS_10014802, utmzone] = deg2utm(57.0480133, 9.948805);
[yGNSS_10014520, xGNSS_10014520, utmzone] = deg2utm(57.0480467, 9.94856);

% Beacons placering
beaconPosList = [
    xGNSS_10014565, yGNSS_10014565;
    xGNSS_10014802, yGNSS_10014802;
    xGNSS_10014520, yGNSS_10014520
    ];


sigma_shadowing = 8;
R_ble = sigma_shadowing^2;

n = 1.28;
A = -77.11;

last_rssi = containers.Map('KeyType', 'double', 'ValueType', 'double');
last_time = containers.Map('KeyType', 'double', 'ValueType', 'any');

% EKF loop
for k = 2:N
    dt=dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q_bias = compute_Q_bias(x_ekf_all(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    [x_pred, P_pred, Q] = ekf_prediction_step(x_ekf_all(:,k-1), P_pred, u, Q_bias, dt, f_bias, computeF_bias);
    x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

    if k == 2
        gnss_used = false(height(gnssTable), 1);
    end

    for idx_gnss = 1:height(gnssTable)
        if gnss_used(idx_gnss)
            continue;
        end

        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];
            [x_pred, P_pred, y_tilde] = ekf_gnss_correction(x_pred, P_pred, z, H_gnss, R_gnss);
            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction
            gnss_used(idx_gnss) = true;
            gnss_hits = gnss_hits + 1;
            fprintf('[k=%d] GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));

            break;
        end
    end


    allowed_indices = [1,2,3,4,5,6,7,8,9,10,11,12, 13, 14];
    for idx_ble = 1:height(ble_used)
        if ble_used_flags(idx_ble)
            continue;
        end

        deltaT_ble = seconds(ble_used.Timestamp(idx_ble) - imuTable.Timestamp(k));

        if deltaT_ble >= 0 && deltaT_ble < 1

            beacon_id = ble_used.BeaconID(idx_ble);
            rssi_now = ble_used.RSSI(idx_ble);
            time_now = ble_used.Timestamp(idx_ble);

            rssi = rssi_now;
            if isKey(last_rssi, beacon_id)
                time_diff = time_now - last_time(beacon_id);
                if time_diff <= seconds(1)
                    rssi = mean([last_rssi(beacon_id), rssi_now]);
                end
            end
            last_rssi(beacon_id) = rssi_now;
            last_time(beacon_id) = time_now;


            ref = beaconPosList(beacon_id, :);
            [x_pred, P_pred, y_ble, rssi_expected, d_measured] = ekf_ble_correction_rssi(x_pred, P_pred, ref, rssi, R_ble, n, A);
            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

            fprintf('[k=%d] BLE-korrektion #%d fra beacon #%d, ændring = %.2f dB\n, expected = %.2f m,målt = %.2f m', ...
                k, idx_ble, beacon_id, abs(y_ble),rssi_expected, rssi);

            ble_used_flags(idx_ble) = true;
            ble_correction_indices(end+1) = k;
            break;
        end
    end

    x_ekf_all(:,k) = x_pred;

end

[mse_IMU_all, rmse_IMU_all] = compare_estimate_to_gnss( ...
    x_ekf_all(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');
% Vel plot
figure;
plot(imuTable.Timestamp, x_ekf_all(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('ekf with BLE-, GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% Bias plot
figure;
plot(imuTable.Timestamp, x_ekf_all(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('ekf with BLE-, GNSS Correction and Bias State: Bias over time');
set(gca, 'FontSize', 18); grid on;

%ekf plot with beacon circles
[x_offset, y_offset, ax] = plot_ekf_result_with_arrows( ...
    x_ekf_all(1,:), ...
    x_ekf_all(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU_all ...
    );

title(ax, 'EKF with BLE Correction and Bias State: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);

xGNSS_10014565_shifted = xGNSS_10014565 - x_offset;
xGNSS_10014802_shifted = xGNSS_10014802 - x_offset;
xGNSS_10014520_shifted = xGNSS_10014520 - x_offset;
yGNSS_10014565_shifted = yGNSS_10014565 - y_offset;
yGNSS_10014802_shifted = yGNSS_10014802 - y_offset;
yGNSS_10014520_shifted = yGNSS_10014520 - y_offset;

hold on;
plot(xGNSS_10014565_shifted, yGNSS_10014565_shifted, 'kx', ...
    'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Bluetooth transmitter');
plot(xGNSS_10014802_shifted, yGNSS_10014802_shifted, 'kx', ...
    'LineWidth', 2, 'MarkerSize', 5, 'HandleVisibility', 'off');
plot(xGNSS_10014520_shifted, yGNSS_10014520_shifted, 'kx', ...
    'LineWidth', 2, 'MarkerSize', 5, 'HandleVisibility', 'off');

radius = 15;
theta = linspace(0, 2*pi, 100);
fill(xGNSS_10014565_shifted + radius*cos(theta), ...
    yGNSS_10014565_shifted + radius*sin(theta), ...
    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', ...
    'DisplayName', 'BLE coverage (15 m radius)');
fill(xGNSS_10014802_shifted + radius*cos(theta), ...
    yGNSS_10014802_shifted + radius*sin(theta), ...
    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(xGNSS_10014520_shifted + radius*cos(theta), ...
    yGNSS_10014520_shifted + radius*sin(theta), ...
    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

legend('Location', 'best');


%% EKF functions


function [x_pred, P_pred, Q] = ekf_prediction_step(x_prev, P_prev, u, Q, dt, f, computeF)

% prediction
x_pred = f(x_prev, u, dt);
x_pred(4:5) = x_pred(4:5) / norm(x_pred(4:5) + 1e-6);
x_pred_k = x_pred;

% 3.covariance pred
F = computeF(x_prev, u, dt);
P_pred = F * P_prev * F' + Q;
end

function Q = compute_Q(x_prev, alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt)

% Q construction
h = x_prev(4:5);
h = h / norm(h + 1e-6);  % Undgå div-by-zero
H_rot = Rpi2 * (h * h') * Rpi2';

Q = zeros(5,5);
Q(3,3) = alpha_v * (acc_noise(3)^2 * dt^2);
Q_h = gyro_noise(1)^2 * dt^2 * H_rot;
Q(4:5,4:5) = alpha_h * Q_h;

end

function Q = compute_Q_bias(x_prev, alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt,alpha_bias)

% Q construction
h = x_prev(4:5);
h = h / norm(h + 1e-6);
H_rot = Rpi2 * (h * h') * Rpi2';

Q = zeros(6,6);
Q(3,3) = alpha_v * (acc_noise(3)^2 * dt^2);
Q_h = gyro_noise(1)^2 * dt^2 * H_rot;
Q(4:5,4:5) = alpha_h * Q_h;
Q(6,6) = alpha_bias;

end

function [x_upd, P_upd, y_tilde] = ekf_gnss_correction(x_pred, P_pred, z, H, R)
y_tilde = z - H * x_pred;

S = H * P_pred * H' + R;
K = P_pred * H' / S;

x_upd = x_pred + K * y_tilde;
x_upd(4:5) = x_upd(4:5) / norm(x_upd(4:5) + 1e-6);

I = eye(size(P_pred));
P_upd = (I - K * H) * P_pred * (I - K * H)' + K * R * K';
end

function [x_upd, P_upd, y_ble, rssi_expected, d_measured] = ekf_ble_correction_rssi( ...
    x_pred, P_pred, ref_pos, rssi_now, R_ble, n, A)

dx = x_pred(1) - ref_pos(1);
dy = x_pred(2) - ref_pos(2);
d_expected = norm([dx; dy]) + 1e-6;
rssi_expected = A - 10 * n * log10(d_expected);
d_measured = 10^((A - rssi_now) / (10 * n));

y_ble = -(rssi_now - rssi_expected);

H_ble = zeros(1,6);
H_ble(1) = (10 * n / log10(10)) * dx / d_expected^2;
H_ble(2) = (10 * n / log10(10)) * dy / d_expected^2;

S_ble = H_ble * P_pred * H_ble' + R_ble;
K_ble = (P_pred * H_ble') / S_ble;

x_upd = x_pred + K_ble * y_ble;
I = eye(size(P_pred));
P_upd = (I - K_ble * H_ble) * P_pred * (I - K_ble * H_ble)' + K_ble * R_ble * K_ble';
end


%% Functions for comparison

function [mse_total, rmse_total] = compare_estimate_to_gnss(x_est, timestamp_est, x_ref, y_ref, timestamp_ref, name)

if nargin < 6
    name = 'Estimat';
end

x_ref_interp = interp1(timestamp_ref, x_ref, timestamp_est, 'linear', 'extrap');
y_ref_interp = interp1(timestamp_ref, y_ref, timestamp_est, 'linear', 'extrap');

squared_error = (x_est(1,:)' - x_ref_interp).^2 + (x_est(2,:)' - y_ref_interp).^2;

mse_over_time = movmean(squared_error, 5);

mse_total = mean(squared_error, 'omitnan');
rmse_total = sqrt(mse_total);

fprintf('--- %s vs. GNSS ---\n', name);
fprintf('Samlet MSE: %.3f m²\n', mse_total);
fprintf('Samlet RMSE: %.3f m\n', rmse_total);

figure;

subplot(3,1,1)
plot(timestamp_est, x_est(1,:), 'b', timestamp_est, x_ref_interp, 'r--');
xlabel('Tid'); ylabel('x [m]');
legend([name ' estimat'], 'GNSS-reference');
title('x-position');

subplot(3,1,2)
plot(timestamp_est, x_est(2,:), 'b', timestamp_est, y_ref_interp, 'r--');
xlabel('Tid'); ylabel('y [m]');
legend([name ' estimat'], 'GNSS-reference');
title('y-position');

subplot(3,1,3)
plot(timestamp_est, mse_over_time, 'k');
xlabel('Tid'); ylabel('MSE [m²]');
title('Mean Square Error over tid');

sgtitle([name ' vs. GNSS-reference']);
end

%%

function [x_offset, y_offset, ax] = plot_ekf_result_with_arrows( ...
    x_ekf, y_ekf, xGNSS, yGNSS, xRef, yRef, ...
    arrow_step, arrow_length, arrow_width, rmse_val)

x_ekf = movmean(x_ekf, 5);
y_ekf = movmean(y_ekf, 5);

x_offset = min([x_ekf(:); xGNSS(:); xRef(:)]);
y_offset = min([y_ekf(:); yGNSS(:); yRef(:)]);
x_ekf = x_ekf - x_offset;
y_ekf = y_ekf - y_offset;
xGNSS = xGNSS - x_offset;
yGNSS = yGNSS - y_offset;
xRef = xRef - x_offset;
yRef = yRef - y_offset;

fig = figure;
ax = gca;

plot(ax, x_ekf, y_ekf, 'b-', 'LineWidth', 1.5, 'DisplayName', 'EKF position estimate'); hold(ax, 'on');
plot(ax, xGNSS, yGNSS, 'ro', 'DisplayName', 'GNSS position');
plot(ax, xRef, yRef, '-', 'Color', [0.5 0 0.5], 'LineWidth', 1.5, 'DisplayName', 'Geo Tracker Ref');

idx = 1:arrow_step:(length(x_ekf) - 1);
u = x_ekf(idx + 1) - x_ekf(idx);
v = y_ekf(idx + 1) - y_ekf(idx);
mag = sqrt(u.^2 + v.^2);
mag(mag == 0) = eps;
u = u ./ mag;
v = v ./ mag;

for i = 1:length(idx)
    x0 = x_ekf(idx(i));
    y0 = y_ekf(idx(i));
    dx = u(i); dy = v(i);
    perp = [-dy, dx];
    tip = [x0 + dx * arrow_length, y0 + dy * arrow_length];
    base1 = [x0 + perp(1) * arrow_width / 2, y0 + perp(2) * arrow_width / 2];
    base2 = [x0 - perp(1) * arrow_width / 2, y0 - perp(2) * arrow_width / 2];

    fill(ax, [base1(1), base2(1), tip(1)], ...
        [base1(2), base2(2), tip(2)], ...
        'b', 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

idx_ref = 1:5:(length(xRef) - 1);
u_ref = xRef(idx_ref + 1) - xRef(idx_ref);
v_ref = yRef(idx_ref + 1) - yRef(idx_ref);
mag_ref = sqrt(u_ref.^2 + v_ref.^2);
u_ref = u_ref ./ mag_ref;
v_ref = v_ref ./ mag_ref;

for i = 1:length(idx_ref)
    x0 = xRef(idx_ref(i));
    y0 = yRef(idx_ref(i));
    dx = u_ref(i); dy = v_ref(i);
    perp = [-dy, dx];
    tip = [x0 + dx * arrow_length, y0 + dy * arrow_length];
    base1 = [x0 + perp(1) * arrow_width / 2, y0 + perp(2) * arrow_width / 2];
    base2 = [x0 - perp(1) * arrow_width / 2, y0 - perp(2) * arrow_width / 2];

    fill(ax, [base1(1), base2(1), tip(1)], ...
        [base1(2), base2(2), tip(2)], ...
        [0.5 0 0.5], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xlabel(ax, 'Distance East [m]');
ylabel(ax, 'Distance North [m]');
title(ax, 'ekf with GNSS Correction and Bias State');
legend(ax, 'Location', 'best');
axis(ax, 'equal');
grid(ax, 'on');
set(ax, 'FontSize', 18);

drawnow;
xlims = xlim(ax);
ylims = ylim(ax);
x_txt = xlims(2) - 0.01 * range(xlims);
y_txt = ylims(2) - 0.02 * range(ylims);
text(ax, x_txt, y_txt, sprintf('RMSE = %.2f m', rmse_val), ...
    'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

end
