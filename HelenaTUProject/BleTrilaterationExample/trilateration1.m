
%beacons found from Device Data
lat = [57.04868; 57.04876; 57.0487333];
lon = [9.946495; 9.9465683; 9.9464817];

%GNSS location from TU700 data collector
GNSSLocation=[57.0487269, 9.946505];

%Calculated distances RSSI
r1 = 10^((-77.11 - (-89)) / (10 * 2));
r2 = 10^((-77.11 - (-93)) / (10 * 2));
r3 = 10^((-77.11 - (-90)) / (10 * 2));
distances = [r1; r2; r3]

% Convert to UTM 
[x, y, utmzone] = deg2utm(lat, lon);

%% Least squares

beacons = [x(1,:), y(1,:); x(2,:), y(2,:); x(3,:), y(3,:)]; 

positionLS = trilaterationLS(beacons, distances)

%fprintf('Estimated Position using least squares:\n Latitude: %.6f\n Longitude: %.6f\n', LatLS, LonLS);

%% file:///C:/Users/hts/Downloads/A_Comparison_Analysis_of_BLE-Based_Algorithms_for_.pdf

B1 = [x(1), y(1)];
B2 = [x(2), y(2)];
B3 = [x(3), y(3)];

positionEx = trilateration_explicit(B1, B2, B3, r1, r2, r3)

%fprintf('Estimated Position using explicit:\n Latitude: %.6f\n Longitude: %.6f\n', LatEx, LonEx);

%% Display in UTM

figure;
hold on;
axis equal;
grid on;

% Plot beacons
plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Transmitters');

% circles
theta = linspace(0, 2*pi, 360);
for i = 1:3
    cx = x(i) + distances(i) * cos(theta);
    cy = y(i) + distances(i) * sin(theta);
    plot(cx, cy, '--','LineWidth', 2, 'DisplayName', sprintf('Distance from Transmitter %d', i));
end

% Plot GNSS location
[xGNSS, yGNSS] = deg2utm(GNSSLocation(1), GNSSLocation(2));
plot(xGNSS, yGNSS, 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'GNSS Position');

% Plot estimated positions from the three methods
 plot(positionLS(1), positionLS(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Estimated Position Least Squares');
 plot(positionEx(1), positionEx(2), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Estimated Position Algebraic');

legend show;
xlabel('UTM X (meters)');
ylabel('UTM Y (meters)');
title('Trilateration: Transmitters and Estimated Position');
%%
[xGNSS, yGNSS] = deg2utm(GNSSLocation(1), GNSSLocation(2));
% Shift everything so the bottom-left is at (0, 0)
x_offset = min([x; xGNSS; positionLS(1); positionEx(1)]);
y_offset = min([y; yGNSS; positionLS(2); positionEx(2)]);

x_shifted = x - x_offset;
y_shifted = y - y_offset;
xGNSS_shifted = xGNSS - x_offset;
yGNSS_shifted = yGNSS - y_offset;
positionLS_shifted = positionLS - [x_offset; y_offset];
positionEx_shifted = positionEx' - [x_offset; y_offset];

figure;
hold on;
axis equal;
grid on;

plot(x_shifted, y_shifted, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Transmitters');
% circles
theta = linspace(0, 2*pi, 360);
for i = 1:3
    cx = x_shifted(i) + distances(i) * cos(theta);
    cy = y_shifted(i) + distances(i) * sin(theta);
    plot(cx, cy, '--','LineWidth', 2, 'DisplayName', sprintf('Distance from Transmitter %d', i));
end

% Plot GNSS location

plot(xGNSS_shifted, yGNSS_shifted, 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'GNSS Position');

% Plot estimated positions from the three methods
 plot(positionLS_shifted(1),positionLS_shifted(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Estimated Position Least Squares');
 plot(positionEx_shifted(1), positionEx_shifted(2), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Estimated Position Algebraic');

legend show;
xlabel('UTM X (meters)');
ylabel('UTM Y (meters)');
title('Trilateration: Transmitters and Estimated Position');
set(gca, 'FontSize', 18);
%% Error

[xGNSS, yGNSS] = deg2utm(GNSSLocation(1), GNSSLocation(2));

% Calculate errors
errorLS = norm([xGNSS - positionLS(1), yGNSS - positionLS(2)]);
errorEx = norm([xGNSS - positionEx(1), yGNSS - positionEx(2)]);

% Display errors
fprintf('\n--- Positioning Errors from GNSS [meters] ---\n');
fprintf('Least Squares Error: %.2f m\n', errorLS);
fprintf('Explicit Method Error: %.2f m\n', errorEx);


