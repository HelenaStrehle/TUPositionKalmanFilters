clear all;
%% Calibrating IMU
[acc_offset, gyro_offset, acc_noise, gyro_noise] = calibrate_imu('StartNorthNotMove.log');
%% Load data from TU700
[imuTable, gnssTable, bleTable] = RuteLoad('BicycleTest.log');
gnssTable(1,:) = [];  % remove first unreliable gnss

%% Load Geo tracker ref and edit time start + stop time
% Ref fix
RuteRef = GeotrackerRef('BicycleTestRef.kml');
[yGNSSRef, xGNSSRef, utmzone] = deg2utm(RuteRef.Latitude, RuteRef.Longitude);
RuteRef.TimeUTC = datetime(RuteRef.TimeUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
RuteRef.TimeUTC.TimeZone = '';

imuTable.Timestamp.TimeZone = '';
gnssTable.Timestamp.TimeZone = '';
bleTable.Timestamp.TimeZone = '';

t_start = min(gnssTable.Timestamp);
t_stop = max(RuteRef.TimeUTC);

imuTable(imuTable.Timestamp < t_start | imuTable.Timestamp > t_stop, :) = [];
gnssTable(gnssTable.Timestamp < t_start | gnssTable.Timestamp > t_stop, :) = [];
bleTable(bleTable.Timestamp < t_start | bleTable.Timestamp > t_stop, :) = [];

%% Ready data

%gnssTsble to utm koordinates
[yGNSS, xGNSS, utmzone] = deg2utm(gnssTable.Latitude, gnssTable.Longitude);

% Remove offset
[accX_corr, accY_corr, accZ_corr, GyroX_corr, GyroY_corr, GyroZ_corr] = preProcessSensorData(imuTable, acc_offset, gyro_offset);
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

% State transition bias

f_bias = @(x, u, dt) [
    x(1) + x(3)*x(4)*dt;
    x(2) + x(3)*x(5)*dt;
    x(3) + (u(1) - x(6)) * dt;
    x(4) + u(2)*(-x(5))*dt;
    x(5) + u(2)*x(4)*dt;
    x(6);
    ];


%% UKF IMU prediction

% UKF parametre
alpha = 1e-3;
kappa = 0;
beta = 2;
x_ukf = zeros(6, N);
x_ukf(:,1) = [xGNSS(1); yGNSS(1); 0; h0; 0.45];
P_pred = eye(6);

% Q scaling constants
alpha_v = 10000000;
alpha_h = 10000000;
% alpha_v = 1000;
% alpha_h = 100000;
alpha_bias = 1e-6;

for k = 2:N
    dt = dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q_bias(x_ukf(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    % Generate sigma points
    [X_sigma, Wm, Wc] = generate_sigma_points(x_ukf(:,k-1), P_pred, alpha, beta, kappa);

    [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f_bias);
    x_pred(3) = max(x_pred(3),0);
    x_ukf(:,k) = x_pred;
end

[mse_IMU, rmse_IMU] = compare_estimate_to_gnss( ...
    x_ukf(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');

figure;
plot(imuTable.Timestamp, x_ukf(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('Prediction step with IMU: Velocity over time');
set(gca, 'FontSize', 18); grid on;


[x_offset, y_offset, ax] = plot_ukf_result_with_arrows( ...
    x_ukf(1,:), ...
    x_ukf(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1.5, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU ...
    );

title(ax, 'Prediction step with IMU: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);
%% UKF GNSS first itt

% UKF-parametre
alpha = 1e-3;
kappa = 0;
beta  = 2;
x_ukf_gnss = zeros(5, N);
x_ukf_gnss(:,1) = [xGNSS(1); yGNSS(1); 0; h0];
P_pred = eye(5);

% Q scaling constants
alpha_v = 1;
alpha_h = 1;

% gnss check
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;
R_gnss = diag([1.5 1.5])^2
H_gnss = [1 0 0 0 0;
    0 1 0 0 0];

%UKF loop
for k = 2:N
    dt = dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q(x_ukf_gnss(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt);

    [X_sigma, Wm, Wc] = generate_sigma_points(x_ukf_gnss(:,k-1), P_pred, alpha, beta, kappa);

    [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f);

    % GNSS-measurement update
    for idx_gnss = 1:height(gnssTable)
        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];

            % GNSS update step for UKF
            [x_pred, P_pred, y_tilde] = ukf_gnss_update(x_pred, P_pred, X_pred, Wm, Wc, z, H_gnss, R_gnss);

            gnss_hits = gnss_hits + 1;
            fprintf('[k=%d] UKF-GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));
        end
    end

    x_ukf_gnss(:,k) = x_pred;
end


[mse_IMU, rmse_IMU_gnss] = compare_estimate_to_gnss( ...
    x_ukf_gnss(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');

% Velocity plot
figure;
plot(imuTable.Timestamp, x_ukf_gnss(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('UKF with GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

%%Position plot
[x_offset, y_offset, ax] = plot_ukf_result_with_arrows( ...
    x_ukf_gnss(1,:), ...
    x_ukf_gnss(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1.5, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU_gnss ...
    );

title(ax, 'UKF with GNSS Correction: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);
%% UKF GNSS Bias state

% UKF-parametre
alpha = 1e-3;
kappa = 0;
beta  = 2;
x_ukf_gnss_bias = zeros(6, N);
x_ukf_gnss_bias(:,1) = [xGNSS(1); yGNSS(1); 0.0; h0; 0.5];
P_pred = eye(6);

% Q scaling constants
alpha_v = 1000;
alpha_h = 100000;
alpha_bias = 1e-6;

% gnss check
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

R_gnss = diag([4 4])^2
H_gnss = [1 0 0 0 0 0;
    0 1 0 0 0 0];

for k = 2:N
    dt = dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q_bias(x_ukf_gnss_bias(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);
    x_pred(3) = max(x_pred(3),0);
    [X_sigma, Wm, Wc] = generate_sigma_points(x_ukf_gnss_bias(:,k-1), P_pred, alpha, beta, kappa);

    [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f_bias);
    x_pred(3) = max(x_pred(3),0);

    if k == 2
        gnss_used = false(height(gnssTable), 1);
    end

    % GNSS-korrektion
    for idx_gnss = 1:height(gnssTable)
        if gnss_used(idx_gnss)
            continue;
        end

        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];
            [x_pred, P_pred, y_tilde] = ukf_gnss_update(x_pred, P_pred, X_pred, Wm, Wc, z, H_gnss, R_gnss);
            x_pred(3) = max(x_pred(3),0);
            gnss_used(idx_gnss) = true;
            gnss_hits = gnss_hits + 1;

            fprintf('[k=%d] GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));

            break;
        end
    end

    x_ukf_gnss_bias(:,k) = x_pred;
end
Q

[mse_IMU_gnss_bias, rmse_IMU_gnss_bias] = compare_estimate_to_gnss( ...
    x_ukf_gnss_bias(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'UKF');


% Vel plot
figure;
plot(imuTable.Timestamp, x_ukf_gnss_bias(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('UKF with GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% Bias plot
figure;
plot(imuTable.Timestamp, x_ukf_gnss_bias(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('UKFwith GNSS Correction and Bias State: Bias over time');
set(gca, 'FontSize', 18); grid on;

%UKF plot
[x_offset, y_offset, ax] = plot_ukf_result_with_arrows( ...
    x_ukf_gnss_bias(1,:), ...
    x_ukf_gnss_bias(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU_gnss_bias ...
    );

title(ax, 'UKF with GNSS Correction and Bias State: Estimated position', 'FontSize', 22);
set(ax, 'FontSize', 16);


%% UKF ble with bias state
%UKF paramters
alpha = 1e-2;
kappa = 0;
beta  = 2;

%Initialise
x_ukf_ble_bias = zeros(6, N);
x_ukf_ble_bias(:,1) = [xGNSS(1); yGNSS(1); 0.01; h0; 0.5];
P_pred = eye(6);


% Q scaling constants
alpha_v = 1000;
alpha_h = 10;
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
last_rssi = containers.Map('KeyType', 'double', 'ValueType', 'double');
last_time = containers.Map('KeyType', 'double', 'ValueType', 'any');

%BLE parameters
sigma_shadowing = 6;
R_ble = sigma_shadowing^2;
n = 1.28;
A = -77.11;

for k = 2:N
    dt = dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q_bias(x_ukf_ble_bias(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    [X_sigma, Wm, Wc] = generate_sigma_points(x_ukf_ble_bias(:,k-1), P_pred, alpha, beta, kappa);

    [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f_bias);
    x_pred(3) = max(x_pred(3),0);

    %BLE update
    allowed_indices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];  %
    for i = 1:length(allowed_indices)
        idx_ble = allowed_indices(i);
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

            [x_pred, P_pred, y_ble, rssi_expected, d_measured] = ukf_ble_correction_rssi(x_pred, P_pred, X_pred, Wm, Wc, ref, rssi_now, R_ble, n, A);
            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

            fprintf('[k=%d] BLE-korrektion #%d fra beacon #%d, ændring = %.2f dB\n, expected = %.2f m,målt = %.2f m', ...
                k, idx_ble, beacon_id, y_ble,rssi_expected, rssi);

            ble_used_flags(idx_ble) = true;
            ble_correction_indices(end+1) = k;
            break;
        end
    end

    x_ukf_ble_bias(:,k) = x_pred;

end

[mse_IMU_ble, rmse_IMU_ble] = compare_estimate_to_gnss( ...
    x_ukf_ble_bias(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'ukf');

% Vel plot
figure;
plot(imuTable.Timestamp, x_ukf_ble_bias(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('UKF with BLE-, GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% Bias plot
figure;
plot(imuTable.Timestamp, x_ukf_ble_bias(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('UKF with BLE-, GNSS Correction and Bias State: Bias over time');
set(gca, 'FontSize', 18); grid on;

%UKF plot with beacon circles
[x_offset, y_offset, ax] = plot_ukf_result_with_arrows( ...
    x_ukf_ble_bias(1,:), ...
    x_ukf_ble_bias(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1.5, ...            % arrow_length
    2, ...              % arrow_width
    rmse_IMU_ble ...
    );

title(ax, 'UKF with BLE Correction and Bias State: Estimated position', 'FontSize', 22);
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


%% UKF ble and gnss

%UKF parameters
alpha = 1e-2;
kappa = 0;
beta  = 2;

%Initialise
x_ukf_all = zeros(6, N);
x_ukf_all(:,1) = [xGNSS(1); yGNSS(1); 0.0; h0; 0.5];
P_pred = eye(6);

% Q scaling constants
alpha_v = 1000;
alpha_h = 1000;
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

beaconPosList = [
    xGNSS_10014802, yGNSS_10014802;
    xGNSS_10014520, yGNSS_10014520
    xGNSS_10014565, yGNSS_10014565;

    ];

%BLE PArameters
sigma_shadowing = 6;
R_ble = sigma_shadowing^2;

last_rssi = containers.Map('KeyType', 'double', 'ValueType', 'double');
last_time = containers.Map('KeyType', 'double', 'ValueType', 'any');
n = 1.28;
A = -77.11;

% gnss Parameters
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

R_gnss = diag([1.5 1.5])^2
H_gnss = [1 0 0 0 0 0;
    0 1 0 0 0 0];

for k = 2:N
    dt = dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q_bias(x_ukf_all(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    [X_sigma, Wm, Wc] = generate_sigma_points(x_ukf_all(:,k-1), P_pred, alpha, beta, kappa);

    [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f_bias);
    x_pred(3) = max(x_pred(3),0);
    if k == 2
        gnss_used = false(height(gnssTable), 1);
    end

    % GNSS-korrektion
    for idx_gnss = 1:height(gnssTable)
        if gnss_used(idx_gnss)
            continue;
        end

        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];
            [x_pred, P_pred, y_tilde] = ukf_gnss_update(x_pred, P_pred, X_pred, Wm, Wc, z, H_gnss, R_gnss);
            x_pred(3) = max(x_pred(3),0);
            gnss_used(idx_gnss) = true;
            gnss_hits = gnss_hits + 1;

            fprintf('[k=%d] GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));

            break;
        end
    end

    %BLE update
    allowed_indices = [1,2,3,4,5,6,7,8,9,10,11,12, 13, 14];
    for idx_ble = 1:height(ble_used)
        if ble_used_flags(idx_ble)
            continue;
        end

        deltaT_ble = seconds(ble_used.Timestamp(idx_ble) - imuTable.Timestamp(k));

        if deltaT_ble >= 0 && deltaT_ble < 0.1

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

            [x_pred, P_pred, y_ble, rssi_expected, d_measured] = ukf_ble_correction_rssi(x_pred, P_pred, X_pred, Wm, Wc, ref, rssi_now, R_ble, n, A);
            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

            fprintf('[k=%d] BLE-korrektion #%d fra beacon #%d, ændring = %.2f dB\n, expected = %.2f m,målt = %.2f m', ...
                k, idx_ble, beacon_id, abs(y_ble),rssi_expected, rssi);

            ble_used_flags(idx_ble) = true;
            ble_correction_indices(end+1) = k;
            break;
        end
    end

    x_ukf_all(:,k) = x_pred;
end


[mse_IMU_all, rmse_IMU_all] = compare_estimate_to_gnss( ...
    x_ukf_all(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'ukf');

% vel plot
figure;
plot(imuTable.Timestamp, x_ukf_all(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Velocity [m/s]');
title('UKF with GNSS Correction and Bias State: Velocity over time');
set(gca, 'FontSize', 18); grid on;

% bias plot
figure;
plot(imuTable.Timestamp, x_ukf_all(6,:), 'b-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Bias');
title('UKF with GNSS Correction and Bias State: Bias over time');
set(gca, 'FontSize', 18); grid on;

%UKF plot with BLE circles
[x_offset, y_offset, ax] = plot_ukf_result_with_arrows( ...
    x_ukf_all(1,:), ...
    x_ukf_all(2,:), ...
    xGNSS, yGNSS, ...
    xGNSSRef, yGNSSRef, ...
    1800, ...           % arrow_step
    1.5, ...            % arrow_length
    1, ...              % arrow_width
    rmse_IMU_all ...
    );

title(ax, 'UKF with BLE Correction and Bias State: Estimated position', 'FontSize', 22);
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


%% Functions for UKF
function [x_upd, P_upd, y_tilde, rssi_expected, d_measured] = ukf_ble_correction_rssi(x_pred, P_pred, X_pred, Wm, Wc, ref_pos, rssi_now, R_ble, n, A)

n_sigma = size(X_pred, 2);
Z_sigma = zeros(1, n_sigma);
d_measured = 10^((A - rssi_now) / (10 * n));

for i = 1:n_sigma
    dx = X_pred(1,i) - ref_pos(1);
    dy = X_pred(2,i) - ref_pos(2);
    d_exp = norm([dx; dy]) + 1e-6;
    Z_sigma(i) = A - 10 * n * log10(d_exp);
end


z_pred = Z_sigma * Wm';
rssi_expected = z_pred;

y_tilde = -(rssi_now - z_pred);

P_zz = R_ble;
for i = 1:n_sigma
    dz = Z_sigma(i) - z_pred;
    P_zz = P_zz + Wc(i) * (dz * dz');
end

P_xz = zeros(size(x_pred,1), 1);
for i = 1:n_sigma
    dx = X_pred(:,i) - x_pred;
    dz = Z_sigma(i) - z_pred;
    P_xz = P_xz + Wc(i) * (dx * dz);
end

K = P_xz / P_zz;

x_upd = x_pred + K * y_tilde;
P_upd = P_pred - K * P_zz * K';
end

%% Funtions for UKF

function Q = compute_Q(x_prev, alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt)

h = x_prev(4:5);
h = h / norm(h + 1e-6);
H_rot = Rpi2 * (h * h') * Rpi2';

Q = zeros(5,5);
Q(3,3) = alpha_v * (acc_noise(3)^2 * dt^2);
Q_h = gyro_noise(1)^2 * dt^2 * H_rot;
Q(4:5,4:5) = alpha_h * Q_h;

end

function Q = compute_Q_bias(x_prev, alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt,alpha_bias)

h = x_prev(4:5);
h = h / norm(h + 1e-6);
H_rot = Rpi2 * (h * h') * Rpi2';

Q = zeros(5,5);
Q(3,3) = alpha_v * (acc_noise(3)^2 * dt^2);
Q_h = gyro_noise(1)^2 * dt^2 * H_rot;
Q(4:5,4:5) = alpha_h * Q_h;
Q(6,6) = alpha_bias;

end

function [X, Wm, Wc] = generate_sigma_points(x, P, alpha, beta, kappa)
n = numel(x);
lambda = alpha^2 * (n + kappa) - n;
gamma = sqrt(n + lambda);

X = zeros(n, 2*n + 1);
Wm = zeros(1, 2*n + 1);
Wc = zeros(1, 2*n + 1);

X(:,1) = x;
Wm(1) = lambda / (n + lambda);
Wc(1) = Wm(1) + (1 - alpha^2 + beta);


P = (P + P') / 2;
eigvals = eig(P);
if any(eigvals <= 0)
    P = P + eye(size(P)) * (abs(min(eigvals)) + 1e-6);  % Make positive definit
end

sqrtP = chol(P, 'lower');
for i = 1:n
    X(:, i+1)     = x + gamma * sqrtP(:,i);
    X(:, i+1+n)   = x - gamma * sqrtP(:,i);
    Wm(i+1)       = 1 / (2*(n + lambda));
    Wm(i+1+n)     = Wm(i+1);
    Wc(i+1)       = Wm(i+1);
    Wc(i+1+n)     = Wm(i+1);
end
end

function [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f)


n_sigma = size(X_sigma, 2);
X_pred = zeros(size(X_sigma));

for i = 1:n_sigma
    X_pred(:,i) = f(X_sigma(:,i), u, dt);
    X_pred(4:5,i) = X_pred(4:5,i) / norm(X_pred(4:5,i) + 1e-6);
end

x_pred = X_pred * Wm';
x_pred(4:5) = x_pred(4:5) / norm(x_pred(4:5) + 1e-6);

P_pred = zeros(size(Q));
for i = 1:n_sigma
    dx = X_pred(:,i) - x_pred;
    P_pred = P_pred + Wc(i) * (dx * dx');
end
P_pred = P_pred + Q;
end


function [x_upd, P_upd, y_tilde] = ukf_gnss_update(x_pred, P_pred, X_pred, Wm, Wc, z, H, R)
Z_sigma = H * X_pred;

% Expected
z_pred = Z_sigma * Wm';

% Innovationskovarians
P_zz = R;
for i = 1:size(X_pred,2)
    dz = Z_sigma(:,i) - z_pred;
    P_zz = P_zz + Wc(i) * (dz * dz');
end

%Cross covariance
P_xz = zeros(size(x_pred,1), size(z,1));
for i = 1:size(X_pred,2)
    dx = X_pred(:,i) - x_pred;
    dz = Z_sigma(:,i) - z_pred;
    P_xz = P_xz + Wc(i) * (dx * dz');
end

% Kalman gain
K = P_xz / P_zz;

% Innovation
y_tilde = z - z_pred;

% Opdater
x_upd = x_pred + K * y_tilde;
P_upd = P_pred - K * P_zz * K';
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


function [x_offset, y_offset, ax] = plot_ukf_result_with_arrows( ...
    x_ukf, y_ukf, xGNSS, yGNSS, xRef, yRef, ...
    arrow_step, arrow_length, arrow_width, rmse_val)

x_ukf = movmean(x_ukf, 5);
y_ukf = movmean(y_ukf, 5);

x_offset = min([x_ukf(:); xGNSS(:); xRef(:)]);
y_offset = min([y_ukf(:); yGNSS(:); yRef(:)]);
x_ukf = x_ukf - x_offset;
y_ukf = y_ukf - y_offset;
xGNSS = xGNSS - x_offset;
yGNSS = yGNSS - y_offset;
xRef = xRef - x_offset;
yRef = yRef - y_offset;

fig = figure;
ax = gca;

plot(ax, x_ukf, y_ukf, 'b-', 'LineWidth', 1.5, 'DisplayName', 'UKF position estimate'); hold(ax, 'on');
plot(ax, xGNSS, yGNSS, 'ro', 'DisplayName', 'GNSS position');
plot(ax, xRef, yRef, '-', 'Color', [0.5 0 0.5], 'LineWidth', 1.5, 'DisplayName', 'Geo Tracker Ref');

idx = 1:arrow_step:(length(x_ukf) - 1);
u = x_ukf(idx + 1) - x_ukf(idx);
v = y_ukf(idx + 1) - y_ukf(idx);
mag = sqrt(u.^2 + v.^2);
mag(mag == 0) = eps;
u = u ./ mag;
v = v ./ mag;

for i = 1:length(idx)
    x0 = x_ukf(idx(i));
    y0 = y_ukf(idx(i));
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
title(ax, 'UKF with GNSS Correction and Bias State');
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
%%

%% UKF ble and gnss

%UKF parameters
alpha = 1e-2;
kappa = 0;
beta  = 2;

%Initialise
x_ukf_all = zeros(6, N);
x_ukf_all(:,1) = [xGNSS(1); yGNSS(1); 0.01; h0; 0.5];
P_pred = eye(6);

% Q scaling constants
alpha_v = 1000;
alpha_h = 1000;
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

beaconPosList = [
    xGNSS_10014802, yGNSS_10014802;
    xGNSS_10014520, yGNSS_10014520
    xGNSS_10014565, yGNSS_10014565;

    ];

%BLE PArameters
sigma_shadowing = 6;
R_ble = sigma_shadowing^2;

last_rssi = containers.Map('KeyType', 'double', 'ValueType', 'double');
last_time = containers.Map('KeyType', 'double', 'ValueType', 'any');
n = 1.28;
A = -77.11;

% gnss Parameters
gnss_used = false(height(gnssTable), 1);
gnss_hits = 0;

R_gnss = diag([1.5 1.5])^2
H_gnss = [1 0 0 0 0 0;
    0 1 0 0 0 0];
update_type = zeros(1, N);  % 0 = prediction only, 1 = BLE, 2 = GNSS, 3 = both

for k = 2:N
    dt = dt_all(k-1);
    u = [accZ_corr(k); GyroX_corr(k)];

    Q = compute_Q_bias(x_ukf_all(:,k-1), alpha_v, alpha_h, Rpi2, acc_noise, gyro_noise, dt, alpha_bias);

    [X_sigma, Wm, Wc] = generate_sigma_points(x_ukf_all(:,k-1), P_pred, alpha, beta, kappa);

    [x_pred, P_pred, X_pred] = ukf_prediction_step(X_sigma, Wm, Wc, Q, dt, u, f_bias);
    x_pred(3) = max(x_pred(3),0);
    if k == 2
        gnss_used = false(height(gnssTable), 1);
    end
    gnss_updated = false;
    % GNSS-korrektion
    for idx_gnss = 1:height(gnssTable)
        if gnss_used(idx_gnss)
            continue;
        end

        deltaT = seconds(gnssTable.Timestamp(idx_gnss) - imuTable.Timestamp(k));
        if deltaT >= 0 && deltaT < 0.5
            z = [xGNSS(idx_gnss); yGNSS(idx_gnss)];
            [x_pred, P_pred, y_tilde] = ukf_gnss_update(x_pred, P_pred, X_pred, Wm, Wc, z, H_gnss, R_gnss);
            x_pred(3) = max(x_pred(3),0);
            gnss_used(idx_gnss) = true;
            gnss_updated = true;
            gnss_hits = gnss_hits + 1;

            fprintf('[k=%d] GNSS-korrektion med #%d, ændring = %.2f m\n', ...
                k, idx_gnss, norm(y_tilde));

            break;
        end
    end

    %BLE update
    allowed_indices = [1,2,3,4,5,6,7,8,9,10,11,12, 13, 14];
    ble_updated = false;
    for idx_ble = 1:height(ble_used)
        if ble_used_flags(idx_ble)
            continue;
        end

        deltaT_ble = seconds(ble_used.Timestamp(idx_ble) - imuTable.Timestamp(k));

        if deltaT_ble >= 0 && deltaT_ble < 0.1

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

            [x_pred, P_pred, y_ble, rssi_expected, d_measured] = ukf_ble_correction_rssi(x_pred, P_pred, X_pred, Wm, Wc, ref, rssi_now, R_ble, n, A);
            x_pred(3) = max(x_pred(3),0); % Have to cap velocity to keep direction

            fprintf('[k=%d] BLE-korrektion #%d fra beacon #%d, ændring = %.2f dB\n, expected = %.2f m,målt = %.2f m', ...
                k, idx_ble, beacon_id, abs(y_ble),rssi_expected, rssi);

            ble_used_flags(idx_ble) = true;
            ble_updated = true;
            ble_correction_indices(end+1) = k;
            break;
        end
    end

    x_ukf_all(:,k) = x_pred;
    if gnss_updated && ble_updated
        update_type(k) = 3;
    elseif gnss_updated
        update_type(k) = 2;
    elseif ble_updated
        update_type(k) = 1;
    else
        update_type(k) = 0;
    end
end

colors = {[0.7 0.7 0.7], [0 0.6 0], [1 0 0], [0.5 0 0.5]};  % gray, green, red, purple
labels = {'Prediction only', 'BLE correction', 'GNSS correction', 'Both corrections'};

figure; hold on;
shown = false(1, 4);  % For types 0 to 3

for i = 0:3
    idx = find(update_type(2:end) == i);
    for j = 1:length(idx)
        k = idx(j) + 1;
        if ~shown(i+1)
            h = plot(imuTable.Timestamp([k-1 k]), x_ukf_all(3, [k-1 k]), '-', ...
                'Color', colors{i+1}, 'LineWidth', 2, ...
                'DisplayName', labels{i+1});
            shown(i+1) = true;
        else
            plot(imuTable.Timestamp([k-1 k]), x_ukf_all(3, [k-1 k]), '-', ...
                'Color', colors{i+1}, 'LineWidth', 2, ...
                'HandleVisibility', 'off');  % no legend for repeats
        end
    end
end

xlabel('Time'); ylabel('Velocity [m/s]');
title('UKF with GNSS- and BLE corrections: Velocity over time with correction types');
legend('Location', 'best');
set(gca, 'FontSize', 20); grid on;


[mse_IMU_all, rmse_IMU_all] = compare_estimate_to_gnss( ...
    x_ukf_all(1:2,:), imuTable.Timestamp, xGNSSRef, yGNSSRef, RuteRef.TimeUTC, 'EKF');
