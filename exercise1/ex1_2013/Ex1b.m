%% Part 2
% Mobile GPS receiver

close all; clear all;
%% Data conversion
DATA = load('Dokument/intelligent vehicles/exercises/exercise1/ex1_2013/gps_ex1_morningdrive2012.txt');
Longitude = DATA(:,4); % read all rows in column 4
Latitude  = DATA(:,3); % read all rows in column 3
Timestamp = DATA(:,2); % time in seconds

figure, plot(Longitude,Latitude);
title('Coordinates');
xlabel('Longitude [dddmm.mmmm]');
ylabel('Latitude [ddmm.mmmm]');

% 1. longitude and latitude angles from NMEA-0183 into degrees
[LongDeg, LatDeg] = nmea0183_to_degrees(Longitude, Latitude);

% 2. longitude and latitude angles from NMEA-0183 into degrees
[X, Y] = degrees_to_meters(LongDeg, LatDeg);

F_lon = 62393;
F_lat = 111342;

% Adjust map offset
X = X - min(X);
Y = Y - min(Y);

figure, plot(X, Y);
title('Position in meters');
xlabel('X [m]');
ylabel('Y [m]');



%% Error estimation
mu_X = mean(X);                     % Position mean
mu_Y = mean(Y);

e_X = X - mu_X;                     % Position error
e_Y = Y - mu_Y;

d = sqrt(e_X.^2 + e_Y.^2);
[a, b] = max(d);            % Maximum position error
SIGMA = cov([e_X, e_Y]);    % Covariance matrix

figure, scatter(X, Y, 5, d);
colorbar;
hold
plot(X(b), Y(b), 'rx', 'LineWidth', 2);
title("Position error")
xlabel("X [m]")
ylabel("Y [m]")

%% Maximum velocity
t0 = Timestamp(1);                  % Start time in seconds
v = sqrt(diff(X).^2 + diff(Y).^2) * 3.6;  % Speed in km/h
mu_v = mean(v);                     % Velocity mean
[v_max, v_max_index] = max(v);      % Maximum speed in km/h and corresponding index
    
figure, scatter(X(2:826), Y(2:826), 5, v);
title(colorbar,'km/h');
hold
plot(X(v_max_index), Y(v_max_index), 'rx', 'LineWidth', 2);
title("Velocity")
xlabel("X [m]")
ylabel("Y [m]")

%% Velocity error
e_v = v - mu_v;                     % Velocity error
[e_v_max, e_v_max_index] = max(e_v) % Maximum velocity error 
figure, plot((1:825)', e_v)
title('Velocity error')
xlabel('Time [s]')
ylabel('Error [km/h]')
hold
plot(e_v_max_index, e_v_max,'rx', 'LineWidth', 2);

%% Error correlation

R = randn(1, length(X));
cN = xcorr(R - mean(R), 'coeff');
cv = xcorr(e_v, 'coeff');
cd = xcorr(d, 'coeff');
cX = xcorr(e_X, 'coeff');
cY = xcorr(e_Y, 'coeff');

figure, subplot(4,1,1);
plot(d); 
title('Position error')
xlabel('Time [s]');
ylabel('Error [m]')
subplot(4,1,2);
plot(cd);
hold;
plot(cN);
title('Position error correlation')
xlabel('Time [s]');
ylabel('Correlation')
subplot(4,1,3);
plot(e_v);
title('Velocity error')
xlabel('Time [s]');
ylabel('Error [km/h]')
subplot(4,1,4);
plot(cv);
hold;
plot(cN)
title('Velocity error correlation')
xlabel('Time [s]');
ylabel('Correlation');
%% Heading
heading = atan2d(diff(Y), diff(X));
figure, subplot(2, 1, 1)
plot(X,Y,'HandleVisibility', 'off');
hold on;
plot(X(215:324), Y(215:324), 'LineWidth', 3, 'DisplayName', 'Laholmsvägen')
plot(X(233:324), Y(233:324), 'LineWidth', 2, 'DisplayName', 'Laholmsvägen shortened')
legend;
hold off;
title('Position in meters');
xlabel('X [m]');
ylabel('Y [m]');

subplot(2, 1, 2)
plot(heading)
hold on;
plot(215:324, heading(215:324),'LineWidth', 3)
plot(233:324, heading(233:324), 'LineWidth', 2)
title('Heading')
xlabel('Time [s]')
ylabel('Heading [°]')

mu_heading1 = mean(heading(215:324));
e_heading1 = heading(215:324) - mu_heading1;
mu_heading2 = mean(heading(233:324));
e_heading2 = heading(233:324) - mu_heading2;

figure, plot((215:324)', e_heading1);
hold
plot((233:324)', e_heading2);
xlim([215 324])
title('Heading error')
xlabel('Time [s]')
ylabel('Heading [°]')
legend('Laholmsvägen', 'Laholmsvägen shortened')


sigma_heading1 = var(heading(215:324));
sigma_heading2 = var(heading(233:324));
