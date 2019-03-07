% Main script that reads controller data and laser data
clear all;
close all;
%warning('off','all')
format longG

% Coordinates of the ref. nodes
REF = [1920 9470;   % 01
       10012 8179;  % 02
       9770 7590;   % 03
       11405 7228;  % 04
       11275 6451;  % 05
       11628 6384.5;% 06
       11438 4948;  % 07
       8140 8274;   % 08
       8392 8486;   % 09
       3280 2750;   % 10
       7250 2085;   % 11
       9990 1620;   % 12
       7485 3225;   % 13
       9505 3893;   % 14
       9602 4278;   % 15
       10412 4150;  % 16
       4090 7920;   % 17
       8010 5290;   % 18
       8255 6099;   % 19
       7733 6151;   % 20
       7490 6136;   % 21
       7061 5420;   % 22
       7634 5342];  % 23

% LINE SEGMENT MODEL (Indeces in the REF-vector)
LINES = [1 8;       % L01
         9 2;       % L02
         2 3;       % L03
         3 4;       % L04
         4 5;       % L05
         5 6;       % L06
         6 7;       % L07
         17 10;     % L08
         10 12;     % L09
         11 13;     % L10
         12 16;     % L11
         16 15;     % L12
         15 14;     % L13
         19 21;     % L14
         22 18;     % L15
         20 23];    % L16
     
% Make the LINEMODEL matrix
LINEMODEL = [REF(LINES(:,1),1:2) REF(LINES(:,2),1:2)];
%disp(LINEMODEL)
         
% Control inputs (velocity and steering angle)
CONTROL = load('control_joy.txt');

% Laser data
LD_HEAD   = load('laser_header.txt');
LD_ANGLES = load('laser_angles.txt');
LD_MEASUR = load('laser_measurements.txt');
LD_FLAGS  = load('laser_flags.txt');

[no_inputs, co] = size(CONTROL);

% Robots initial position
X(1) = CONTROL(1,4);
Y(1) = CONTROL(1,5);
A(1) = CONTROL(1,6);
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];

% Position without Cox
XX(1) = CONTROL(1,4);
YY(1) = CONTROL(1,5);
AA(1) = CONTROL(1,6);

% Uncertainty settings for dead reckoning
SIGMAv = 18;
SIGMAalfa = 0.037;
SIGMAT = 0.025;

% Figures
fig_path = figure;
fig_env = figure;
fig_Cox = figure;
fig_error = figure;

% To open the figures at different locations to the right
movegui(fig_env,[-20,-20]);
movegui(fig_path,[-20,30]);
movegui(fig_Cox,[-590,-20]);
movegui(fig_error,[-590,30]);

% Plot the line model
figure(fig_env);
plot_line_segments(REF, LINES, 1);
axis([0 12000 1000 10000]);

% Scan index and Cox counter
scan_idx = 1;
cox_counter = 0;

%% 
for kk = 2:no_inputs
    % Check if we should get a position fix, i.e. if the time stamp of the
    % next laser scan is the same as the time stamp of the control input values
    if LD_HEAD(scan_idx,1) == CONTROL(kk,1)
        % Mark the position where the position fix is done - and the size
        % of the position fix to be found
        figure(fig_path);
        hold on
        plot(X(kk-1), Y(kk-1), 'ro');
        axis([5000 10000 4000 8000]);
        plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'ro');
        plot([X(kk-1) CONTROL(kk-1,4)], [Y(kk-1) CONTROL(kk-1,5)], 'r');
        hold off;
        
        % Get the position fix - Use data points that are ok, i.e. with
        % flags = 0
        DATAPOINTS = find(LD_FLAGS(scan_idx,:) == 0);
        angs = LD_ANGLES(scan_idx, DATAPOINTS);
        meas = LD_MEASUR(scan_idx, DATAPOINTS);
        
        % Plot schematic picture of the Snowhite robot
        alfa = 660;
        beta = 0;
        gamma = -90*pi/180;
        figure(fig_env);
        plot_line_segments(REF, LINES, 1);
        plot_threewheeled_laser([X(kk-1) Y(kk-1) A(kk-1)]', 100, 612, 2, CONTROL(kk-1,6), 150, 50, 680, alfa, beta, gamma, angs, meas, 1);
        axis([0 12000 1000 10000]);
        
        % Call the Cox algorithm
        figure(fig_Cox);
        cox_counter = cox_counter +1;
        [dx, dy, da, C] = Cox_LineFit_v(angs, meas, [X(kk-1) Y(kk-1) A(kk-1)]', LINEMODEL, cox_counter, REF, LINES);
        %dx = 0; dy = 0; da = 0; C = [0 0 0; 0 0 0; 0 0 0];
         
        %% From Ex 3 with no Kalman filter
        % Update the position, i.e. X(kk-1), Y(kk-1), A(kk-1) and C(kk-1)
        %X(kk-1) = X(kk-1) + dx;
        %Y(kk-1) = Y(kk-1) + dy;
        %A(kk-1) = mod(A(kk-1) + da, 2*pi);
        %P(kk-1,1:9) = [C(1,1:3),C(2,1:3),C(3,1:3)];
        
        %% New part with Kalman filter  
        % Dead reckoning position and covariance matrix
        XDR = X(kk-1);
        YDR = Y(kk-1);
        ADR = A(kk-1);
        XhatDR = [XDR YDR ADR]';
        PDR = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];
        
        % True position and covariance matrix of the added noise
        VarianceXTrue = 10; % 10 mm
        VarianceYTrue = 10; % 10 mm
        VarianceATrue = 1*pi/180; % 1 degree
        XTrue = CONTROL(kk-1,4) + randn*VarianceXTrue;
        YTrue = CONTROL(kk-1,5) + randn*VarianceYTrue;
        ATrue = CONTROL(kk-1,6) + randn*VarianceATrue;
        XhatTrue = [XTrue YTrue ATrue]';
        PTrue = [VarianceXTrue 0 0; 0 VarianceYTrue 0; 0 0 VarianceATrue];
       
        % Cox algorithm position and covariance matrix
        XCox = X(kk-1) + dx;
        YCox = Y(kk-1) + dy;
        ACox = A(kk-1) + da;
        XhatCox = [XCox YCox ACox]';
        PCox = C;
        
        % The Kalman filter DR+True. Use this or DR+Cox
%         Xhat = PTrue*(pinv(PTrue+PDR))*XhatDR + PDR*(pinv(PTrue+PDR))*XhatTrue;
%         Phat = pinv(pinv(PTrue)+pinv(PDR));
        
        % The Kalman filter DR+Cox
        Xhat = PCox*(pinv(PCox+PDR))*XhatDR + PDR*(pinv(PCox+PDR))*XhatCox;
        Phat = pinv(pinv(PCox)+pinv(PDR));
        
        % Update position and covariance matrix
        X(kk-1) = Xhat(1);
        Y(kk-1) = Xhat(2);
        A(kk-1) = mod(Xhat(3), 2*pi);
        P(kk-1,1:9) = [Phat(1,1:3),Phat(2,1:3),Phat(3,1:3)];
        %% End of the Kalman filter part
        
        % Next time use the next scan
        scan_idx = mod(scan_idx, max(size(LD_HEAD))) + 1;
    end
    %% End of the part running when there is a laser measurement
    
    % Mark the estimated position with Kalman filter
    figure(fig_path);
    hold on;
    plot(X(kk-1), Y(kk-1), 'b.');
    
    % Mark the position without Cox
    plot(XX(kk-1), YY(kk-1), 'r.');
    
    % Mark the true (from the LaserWay system) position
    plot(CONTROL(kk-1,4), CONTROL(kk-1,5), 'k.');
    axis([5000 10000 4000 8000]);
    hold off;
    
    % To plot the errors with and without Cox
    if mod(kk-1, 10)==0
        ErrCox = sqrt((X(kk-1)-CONTROL(kk-1,4))^2+(Y(kk-1)-CONTROL(kk-1,5))^2);
        ErrNoCox = sqrt((XX(kk-1)-CONTROL(kk-1,4))^2+(YY(kk-1)-CONTROL(kk-1,5))^2);
        figure(fig_error);
        hold on
        plot(kk-1, ErrCox, 'b.')
        plot(kk-1, ErrNoCox, 'r.')
        axis([0 4050 0 250]);
        title('Distance errors with Kalman filter (blue) and without (red)')
        hold off
    end
    
    % To estimate the new position based on the control inputs
    v = CONTROL(kk-1,2);
    a = CONTROL(kk-1,3);
    T = 0.050;
    L = 680;
    
    % Original
    %X(kk) = X(kk-1) + cos(a)*v*cos(A(kk-1))*T;
    %Y(kk) = Y(kk-1) + cos(a)*v*sin(A(kk-1))*T;
    %A(kk) = A(kk-1) + sin(a)*v*T/L;
    
    %% This part is from Exercise 2   
    % Predict the new state variables (World co-ordinates)
    X(kk) = X(kk-1) + cos(a)*v*T*cos(A(kk-1)+(v*sin(a)*T)/(2*L));
    Y(kk) = Y(kk-1) + cos(a)*v*T*sin(A(kk-1)+(v*sin(a)*T)/(2*L));
    A(kk) = mod(A(kk-1) + sin(a)*v*T/L, 2*pi);
   
    % To predict the new uncertainty in the state variables (Error prediction)
    dD = v*cos(a)*T;
    dA = v*sin(a)*T/L;
    
    % Uncertainty in state variables at time k-1 [3x3] 
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];   

    Jxya = [1 0 -dD*sin(A(kk-1)+dA/2);
            0 1 dD*cos(A(kk-1)+dA/2);
            0 0 1];
        
    CvalfaT = [SIGMAv^2 0 0;
               0 SIGMAalfa^2 0;
               0 0 SIGMAT^2];
    
    JvalfaT = [T*cos(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)-(T^2*v*sin(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)*sin(a))/(2*L), -T*v*cos(A(kk-1)+(T*v*sin(a))/(2*L))*sin(a)-(T^2*v^2*sin(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)^2)/(2*L), v*cos(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)-(T*v^2*sin(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)*sin(a))/(2*L);
               T*sin(A(kk-1) + (T*v*sin(a))/(2*L))*cos(a) + (T^2*v*cos(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)*sin(a))/(2*L), (T^2*v^2*cos(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)^2)/(2*L)-T*v*sin(A(kk-1)+(T*v*sin(a))/(2*L))*sin(a), v*sin(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)+(T*v^2*cos(A(kk-1)+(T*v*sin(a))/(2*L))*cos(a)*sin(a))/(2*L);
              (T*sin(a))/L, (T*v*cos(a))/L, (v*sin(a))/L];
        
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Jxya*Cxya_old*Jxya' + JvalfaT*CvalfaT*JvalfaT';
    
    % Store the new co-variance matrix
    P(kk,1:9) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];
    % End of part from Exercise 2
    %%
    
    % Update position without Cox
    XX(kk) = XX(kk-1) + cos(a)*v*cos(AA(kk-1))*T;
    YY(kk) = YY(kk-1) + cos(a)*v*sin(AA(kk-1))*T;
    AA(kk) = AA(kk-1) + sin(a)*v*T/L;
    
end

% Subplot showing errors
ERROR = [X' Y' A'] - CONTROL(:,4:6);
ERROR(:,3) = AngDifference(A',CONTROL(:,6));
ERROR = abs(ERROR);

figure;
subplot(3,1,1); hold on;
plot(ERROR(:,1),'b');
plot(sqrt(P(:,1)),'r'); % one sigma
title('Error X [mm] and uncertainty [std] (red)');
axis([0 4050 0 100]);
hold off;

subplot(3,1,2); hold on;
plot(ERROR(:,2),'b');
plot(sqrt(P(:,5)),'r'); % one sigma
title('Error Y [mm] and uncertainty [std] (red)');
axis([0 4050 0 100]);
hold off;

subplot(3,1,3); hold on;
plot(ERROR(:,3)*180/pi,'b');
plot(sqrt(P(:,9))*180/pi,'r'); % one sigma
title('Error A [degree] and uncertainty [std] (red)');
axis([0 4050 0 5]);
hold off;