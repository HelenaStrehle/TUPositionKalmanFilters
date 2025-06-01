% Clear workspace and command window
clear; clc;

% Define measured distances (in meters)
distances = [1, 2, 5, 7.5 , 10, 15]; 

% Meassured RSSI
rssi_values = {
    [-76, -82, -76, -76, -76, -73, -84, -76, -75, -77, -75, -75, -80, -76, -76, -74, -78, -77, -76, -76, -77, -76, -76, -76, -75, -79, -77, -74, -73, -78, -76, -75, -75, -75, -74, -77, -75, -75, -80, -76, -75, -75, -76, -76, -75, -66, -80, -76, -76], % RSSI at 1m
    [-81, -84, -83, -84, -80, -84, -82, -78, -82, -81, -83, -83, -79, -83, -92, -80, -82, -89, -80, -97, -82, -81, -84, -80], % RSSI at 2m
    [-88, -88, -92, -80, -80, -84], % RSSI at 5m
    [-90, -92, -84, -85, -85, -93, -85, -86, -84, -84, -85, -95, -85, -84, -84, -86, -93, -86, -84, -87, -84, -92, -86], % RSSI at 7.5m
    [-92, -89, -90, -88, -92, -92, -88, -90, -88, -95, -87, -92, -94, -88, -90, -93, -90, -89], % RSSI at 10m
    [-92, -92, -92, -94, -92, -92, -94] % RSSI at 15m
};

% Mean
mean_rssi = cellfun(@mean, rssi_values)

% Convert RSSI to path loss
PL = mean_rssi;
log_distances = log10(distances);

% PL = 10n * log10(d) + A
coeffs = polyfit(log_distances, PL, 1);

n = coeffs(1) / 10; 
A = coeffs(2);  

% Compute estimated path loss using the fitted model
PL_estimated = polyval(coeffs, log_distances);

% difference between measured and estimated path loss
residuals = PL - PL_estimated;

% standard diviation
sigma = std(residuals);

fprintf('Estimated Path Loss Exponent (n): %.2f\n', n);
fprintf('Offset Constant (A): %.2f dB\n', A);
fprintf('Shadowing Standard Deviation (σ): %.2f dB\n', sigma);

% Plot Path Loss vs. log10(Distance)

figure;
plot(log_distances, PL, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot(log_distances, PL_estimated, '-r', 'LineWidth', 1.5);
xlabel('log_{10}(Distance) [m]');
ylabel('Path Loss [dB]');
title('Path Loss vs. Distance');
grid on;
legend('Measured Data', 'Linear Fit');

figure;
plot(log_distances, residuals, 'o-', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('log_{10}(Distance) [m]');
ylabel('Residuals [dB]');
title('Residuals of Path Loss Model');
grid on;

%% Plot sammenligning
PL_est_polyfit = polyval(coeffs, log_distances);

% Beregn estimeret PL for originale afstande
PL_est_polyfit_normal = polyval(coeffs, log10(distances));

%
figure;
plot(distances, PL, 'ko', 'MarkerSize', 8, 'DisplayName', 'Meassured Data'); hold on;
plot(distances, PL_est_polyfit_normal, 'r-', 'LineWidth', 2, 'DisplayName', 'Pass Loss model');
xlabel('Distance [m]');
ylabel('Path Loss [dB]');
title('Path Loss as function of distance');
set(gca, 'FontSize', 20);  % for ticks og akse generelt
xlim([1 inf])
grid on;
legend show;