% Stationary GPS receiver 
% (c) Bj�rn �strand, 2012
close all; clear all; 

DATA = load('gps_ex1_window2012.txt');
Longitude = DATA(:,4); % read all rows in column 4
Latitude  = DATA(:,3); % read all rows in column 3
figure, plot(Longitude,Latitude);
title('Position in NMEA-0183 format');
xlabel('Longitude [dddmm.mmmm]');
ylabel('Latitude [ddmm.mmmm]');

% 3.1 Write a function that transform your longitude and latitude angles 
% from NMEA-0183 into meters
% 1. longitude and latitude angles from NMEA-0183 into degrees

[LongDeg, LatDeg] = nmea0183_to_degrees(Longitude, Latitude);

figure, plot(LongDeg,LatDeg);
title('Position in decimal degrees');
xlabel('Longitude [°]');
ylabel('Latitude [°]');

% 2. longitude and latitude angles from NMEA-0183 into degrees
[X, Y] = degrees_to_meters(LongDeg, LatDeg);

X = X - min(X);
Y = Y - min(Y);

figure, plot(X,Y);
title('Position in meters');
xlabel('X [m]');
ylabel('Y [m]');

% 3.2 Estimate the mean and variance of the position (in x and y)
% Matlab fuctions mean() and var()

mu_X = mean(X);
mu_Y = mean(Y);
e_X = X - mu_X;
e_Y = Y - mu_Y;
sigma_X = std(X);
sigma_Y = std(Y);
sigma2_X = var(e_X);
sigma2_Y = var(e_Y);


d = sqrt(e_X.^2 + e_Y.^2);
[a, b] = max(d);            % Maximum error
SIGMA = cov([e_X, e_Y]);  % Covariance matrix

figure, plot(e_X, e_Y, 'HandleVisibility', 'off')
hold on;
plot(e_X(1:10), e_Y(1:10), '-x', 'DisplayName', 'Error interval 1')
plot(e_X(101:110), e_Y(101:110), '-x', 'DisplayName', 'Error interval 2');
plot(e_X(b), e_Y(b), 'rx', 'DisplayName', 'Maximum error');
plot_uncertainty([0, 0].', SIGMA, 1, 2)
hold off;
title('Position error')
xlabel('X [m]')
ylabel('Y [m]')
legend

figure, histfit(e_X, 30, 'normal')
title('X error')
ylabel('Samples');
xlabel('Distance [m]')
legend('Error', 'Gaussian distribution')
figure, histfit(e_Y, 30, 'normal')
title('Y error')
ylabel('Samples');
xlabel('Distance [m]');
legend('Error', 'Gaussian distribution')
%%
% 3.3 Plot, with respect to time, the errors and the auto-correlation 
% in x and y separately.

R = randn(1, length(X));
cN = xcorr(R - mean(R), 'coeff');
cX = xcorr(e_X, 'coeff');
cY = xcorr(e_Y, 'coeff');

figure, subplot(4,1,1);
plot(e_X); 
title('X error')
xlabel('Time [s]')
ylabel('Distance [m]');
subplot(4,1,2);
plot(e_Y);
title('Y error')
xlabel('Time [s]')
ylabel('Distance [m]');
subplot(4,1,3);
plot(cX, 'HandleVisibility', 'off');
hold
plot(cN, 'DisplayName', 'Random signal')
legend
title('X error correlation')
xlabel('Time [s]')
ylabel('Correlation')
subplot(4,1,4);
plot(cY, 'HandleVisibility', 'off');
hold
plot(cN, 'DisplayName', 'Random signal')
legend
title('Y error correlation')
xlabel('Time [s]')
ylabel('Correlation')
%%
figure, subplot(2,1,1);
plot(cX, 'DisplayName', 'X error');
hold
plot(cN, 'DisplayName', 'Random signal')
xlim([9460 9500])
title('X error correlation')
xlabel('Time [s]')
ylabel('Correlation')
legend('Location', 'southeast')
subplot(2,1,2);
plot(cY, 'DisplayName', 'Y error');
hold
plot(cN, 'DisplayName', 'Random signal')
xlim([9460 9500])
title('Y error correlation')
xlabel('Time [s]')
ylabel('Correlation')
legend('Location', 'southeast')

%%
histfit(R, 100, 'normal')