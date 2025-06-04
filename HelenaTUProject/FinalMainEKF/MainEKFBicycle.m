clear all;
%% Calibrating IMU
[acc_offset, gyro_offset, acc_noise, gyro_noise] = calibrate_imu('StartNorthNotMove.log');
%% Load data from TU700

[imuTable, gnssTable, bleTable] = RuteLoad('BicycleTest.log');
gnssTable(1,:) = [];  % remove first unreliable gnss

%% Load Geo tracker ref and edit time

% Ref fix
RuteRef = GeotrackerRef('BicycleTestRef.kml');
[yGNSSRef, xGNSSRef, utmzone] = deg2utm(RuteRef.Latitude, RuteRef.Longitude);
RuteRef.TimeUTC = datetime(RuteRef.TimeUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
RuteRef.TimeUTC.TimeZone = '';

% Apply time zone settings to other tables
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
x_pred(:,1) = x_ekf_gnss(:,1);


% Q scaling constants
alpha_v = 10000;
alpha_h = 1000;

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
x_ekf_gnss_bias(:,1) = [xGNSS(1); yGNSS(1); 0.0; h0; 0.5];
P_pred = eye(6);

% log
x_pred = zeros(6,N);
x_pred(:,1) = x_ekf_gnss_bias(:,1);

x_pred_all = zeros(6, N);
x_pred_all(:,1) = x_ekf_gnss_bias(:,1); 

% Q scaling constants
alpha_v = 100000;
alpha_h = 1000000000;

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
    x_pred(3) = max(x_pred(3), 0); % Have to cap velocity to keep direction
    
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

%
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
x_ekf_ble_bias(:,1) = [xGNSS(1); yGNSS(1); 0.0; h0; 0.5];
P_pred = eye(6);
x_pred = zeros(6,N);
x_pred(:,1) = x_ekf_ble_bias(:,1);

% Q scaling constants
alpha_v = 100000;
alpha_h = 100000;

alpha_bias = 1e-6;

% Bluetooth set up
mac_to_id = containers.Map({ ...
    '29:B9:A5:C2:3B:58', ...  % SN 10014802
    'DE:52:34:C2:3B:58' ...  % SN 10014520
    '66:97:31:C2:3B:58', ...  % SN 10014565
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


[yGNSS_10014802, xGNSS_10014802, utmzone] = deg2utm(57.0481, 9.9487);
[yGNSS_10014520, xGNSS_10014520, utmzone] = deg2utm(57.04803, 9.94811);
[yGNSS_10014565, xGNSS_10014565, utmzone] = deg2utm(57.0477667, 9.948115);


% Beacons placering
beaconPosList = [
    xGNSS_10014802, yGNSS_10014802;
    xGNSS_10014520, yGNSS_10014520
    xGNSS_10014565, yGNSS_10014565;

    ];


sigma_shadowing = 6;
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
title('EKF with BLE Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% Bias plot
figure;
plot(imuTable.Timestamp, x_ekf_ble_bias(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('EKF with BLE and Bias State: Bias over time');
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
x_ekf_all(:,1) = [xGNSS(1); yGNSS(1); 0.01; h0; 0.5];
P_pred = eye(6);
x_pred = zeros(6,N);
x_pred(:,1) = x_ekf_all(:,1);

% Q scaling constants
alpha_v = 100000;
alpha_h = 1000000;
alpha_bias = 1e-6;

% gnss update initialize
H_gnss = [1 0 0 0 0 0;
    0 1 0 0 0 0];
R_gnss = diag([3 3])^2;

% gnss check
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

% Bluetooth set up
mac_to_id = containers.Map({ ...
    '29:B9:A5:C2:3B:58', ...  % SN 10014802
    'DE:52:34:C2:3B:58' ...  % SN 10014520
    '66:97:31:C2:3B:58', ...  % SN 10014565
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


[yGNSS_10014802, xGNSS_10014802, utmzone] = deg2utm(57.0481, 9.9487);
[yGNSS_10014520, xGNSS_10014520, utmzone] = deg2utm(57.04803, 9.94811);
[yGNSS_10014565, xGNSS_10014565, utmzone] = deg2utm(57.0477667, 9.948115);


% Beacons placering
beaconPosList = [
    xGNSS_10014802, yGNSS_10014802;
    xGNSS_10014520, yGNSS_10014520
    xGNSS_10014565, yGNSS_10014565;

    ];


sigma_shadowing = 16;
R_ble = sigma_shadowing^2;

last_rssi = containers.Map('KeyType', 'double', 'ValueType', 'double');
last_time = containers.Map('KeyType', 'double', 'ValueType', 'any');
n = 1.28;
A = -77.11;

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
