% Stationary GPS receiver 
% (c) Björn Åstrand, 2012
close all; clear all; 

DATA = load('gps_ex1_morningdrive2012.txt');
Longitude = DATA(:,4); % read all rows in column 4
Latitude  = DATA(:,3); % read all rows in column 3

LongDeg = floor(Longitude/100) + (Longitude - floor(Longitude/100)*100)/60;
LatDeg =  floor(Latitude/100) + (Latitude - floor(Latitude/100)*100)/60;

F_lon = 62393; % from table
F_lat = 111342; % from table

X = F_lon * LongDeg;
Y = F_lat * LatDeg; 



% 3.2 Estimate the mean and variance of the position (in x and y)
% Matlab fuctions mean() and var()
meanValueX = mean(X);
meanValueY = mean(Y);
variance = var(X,Y);

%X-Led
ErrorX = X-meanValueX;

%Y-Led
ErrorY = Y-meanValueY;



figure
dMaxSpeed = sqrt((diff(X).^2) + (diff(Y).^2));
[a, b] = max(dMaxSpeed);
plot(3.6*dMaxSpeed)
hold on;
plot(b, 3.6*a,'rs','LineWidth',2)
hold on;
title('Speed km/h')


heading = atan2d(diff(Y), diff(X));
meanHeading = mean(heading);


errorHeading = heading - meanHeading;


figure
plot(errorHeading)


figure
plot(X, Y);
hold on;
plot(X(b), Y(b), 'rs','LineWidth',2)
hold on;


title('Position in meters');
xlabel('Longitude');
ylabel('Latitude');
grid on;



figure
subplot(2, 1, 1)
plot(X,Y);
hold on;
plot(X(240:310), Y(240:310), 'r.','LineWidth',1)
hold on;

title('Position in meters');
xlabel('Longitude');
ylabel('Latitude');


subplot(2, 1, 2)
plot(heading)
hold on;
plot(240:310,heading(240:310),'r','LineWidth',2.5)
