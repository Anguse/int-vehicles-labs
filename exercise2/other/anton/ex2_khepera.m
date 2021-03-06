%
% Dead Reckoning with Khepera Mini Robot
%
% Ola Bengtsson, 2004.02.06
%

clear all;
close all;

% %%% Khepera settings 
WHEEL_BASE = 53;                % [mm]
WHEEL_DIAMETER = 15.3;          % [mm]
PULSES_PER_REVOLUTION = 600;    %
WHEELOMKR = pi*WHEEL_DIAMETER;
MM_PER_PULSE = WHEELOMKR / 600;               % [mm / pulse]
%MM_PER_PULSE = ??; % You should write the correct one, which replaces the one above!


% %%% Uncertainty settings, which are be the same for the left and right encoders
SIGMA_WHEEL_ENCODER = 0.5/12;   % The error in the encoder is 0.5mm / 12mm travelled
% Use the same uncertainty in both of the wheel encoders
SIGMAl = SIGMA_WHEEL_ENCODER;
SIGMAr = SIGMA_WHEEL_ENCODER;


% Load encoder values
ENC = load('khepera.txt');

% Transform encoder values (pulses) into distance travelled by the wheels (mm)

Dr = ENC(1:20:end,2) * MM_PER_PULSE;
Dl = ENC(1:20:end,1) * MM_PER_PULSE;
N = max(size(Dr));
F = 2; %Frequency
% Init Robot Position, i.e. (0, 0, 90*pi/180) and the Robots Uncertainty
X(1) = 0;
Y(1) = 0;
A(1) = 90*pi/180;
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];


%Calcs for variances
    vdD = (SIGMAl+SIGMAr)/4;
    vdTheta = (SIGMAl+SIGMAr)/WHEEL_BASE^2;


% Run until no more encoder values are available
disp('Calculating ...');
CT(1) = 0;
CT(2) = 0;
for kk=2:N,
    % Change of wheel displacements, i.e displacement of left and right wheels
    dDr = Dr(kk) - Dr(kk-1);
    dDl = Dl(kk) - Dl(kk-1);
    
    % Change of relative movements
    %dD = 0;
    %dA = 0.017;
    dD = (dDr + dDl)/2;   % You should write the correct one, which replaces the one above!
    dA = (dDr - dDl)/WHEEL_BASE;   % You should write the correct one, which replaces the one above!
    
    %Loop to calc the compensation term
    % Have to turn large angles for compensation term to matter, for "slow"
    % freq
    if (dA==0),
        CT(1,kk) = 1;
        CT(2,kk) = dA/(pi/180);
    else
        CT(1,kk) = sin(dA/2)/(dA/2);
        CT(2,kk) = dA/(pi/180);
    end
    
    dDc = CT(kk)*dD;
    
    % Calculate the change in X and Y (World co-ordinates)
    dX = dD*cos(A(kk-1)+dA/2);
    dY = dD*sin(A(kk-1)+dA/2);
    %dX = ??;   % You should write the correct one, which replaces the one above!
    %dY = ??;   % You should write the correct one, which replaces the one above!
    dXc =  dDc*cos(A(kk-1)+dA/2)
    dYc =  dDc*sin(A(kk-1)+dA/2);
    % Predict the new state variables (World co-ordinates) with
    % compensation
    Xc(kk) = X(kk-1) + dXc;
    Yc(kk) = Y(kk-1) + dYc;
    Ac(kk) = mod(A(kk-1) + dA, 2*pi);
    
    % Predict the new state variables (World co-ordinates) without
    % compensation
    X(kk) = X(kk-1) + dX;
    Y(kk) = Y(kk-1) + dY;
    A(kk) = mod(A(kk-1) + dA, 2*pi);
    
    
    % Predict the new uncertainty in the state variables (Error prediction)
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];   % Uncertainty in state variables at time k-1 [3x3]    

    Cu =   [vdD 0;0 vdTheta];               % Uncertainty in the input variables [2x2]
    Axya = [1 0 -dD*sin(A(kk-1)+dA/2) ;0 1 dD*cos(A(kk-1)+dA/2);0 0 1];     % Jacobian matrix w.r.t. X, Y and A [3x3]
    Au =   [cos(A(kk-1)+dA/2) (-dD/2)*sin(A(kk-1)+dA/2);sin(A(kk-1)+dA/2) (dD/2)*cos(A(kk-1)+dA/2);0 1];           % Jacobian matrix w.r.t. dD and dA [3x2]
    
    Cu_c =   [vdD 0;0 vdTheta];               % Uncertainty in the input variables [2x2]
    Axya_c = [1 0 -dDc*sin(A(kk-1)+dA/2) ;0 1 dDc*cos(A(kk-1)+dA/2);0 0 1];     % Jacobian matrix w.r.t. X, Y and A [3x3]
    Au_c =   [cos(A(kk-1)+dA/2) (-dDc/2)*sin(A(kk-1)+dA/2);sin(A(kk-1)+dA/2) (dDc/2)*cos(A(kk-1)+dA/2);0 1];           % Jacobian matrix w.r.t. dD and dA [3x2]
    
    
    
    
   
    
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Axya*Cxya_old*Axya' + Au*Cu*Au';
    Cxya_newc = Axya_c*Cxya_old*Axya_c' + Au_c*Cu_c*Au_c';

    % Store the new co-variance matrix
    P(kk,1:9) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];
    
    P_c(kk,1:9) = [Cxya_newc(1,1:3) Cxya_newc(2,1:3) Cxya_newc(3,1:3)];
    
    all_dA(kk) = dA;
    all_dD(kk) = dD;
end;


disp('Plotting ...');

% Plot the path taken by the robot, by plotting the uncertainty in the current position
figure; 
    %plot(X, Y, 'b.');
    title('Path taken by the robot [Wang]');
    xlabel('X [mm] World co-ordinates');
    ylabel('Y [mm] World co-ordinates');
    hold on;
        for kk = 1:2:N,
            C = [P(kk,1:3);P(kk,4:6);P(kk,7:9)];
            plot_uncertainty([X(kk) Y(kk) A(kk)]', C, 1, 2, 'r');
        end;
    hold off;
    axis('equal');
    
%figure; 
    %plot(X, Y, 'b.');
 %   title('Path taken by the robot [Wang] with C');
  %  xlabel('Xc [mm] World co-ordinates');
   % ylabel('Yc [mm] World co-ordinates');
    %hold on;
     %   for kk = 1:2:N,
      %      C_t = [P_c(kk,1:3);P_c(kk,4:6);P_c(kk,7:9)];
       %     plot_uncertainty([Xc(kk) Yc(kk) Ac(kk)]', C_t, 1, 2);
       % end;
    %hold off;
    %axis('equal');    
%____________________________________________________________________
% After the run, plot the results (X,Y,A), i.e. the estimated positions 
%figure;
 %   subplot(3,1,1); plot(X, 'b'); title('X [mm]');
  %  subplot(3,1,2); plot(Y, 'b'); title('Y [mm]');
   % subplot(3,1,3); plot(A*180/pi, 'b'); title('A [deg]');
% After the run, plot the results (Xc,Yc,Ac), i.e. the estimated positions 
figure;
    subplot(3,1,1); plot(Xc, 'g','DisplayName', 'X with compensation'); title('Xc [mm]');
    xlabel('Samples');ylabel('[mm]');
    hold;
    plot(X, 'r','DisplayName', 'X without compensation'); title('X[mm]');
    legend;
    subplot(3,1,2); plot(Yc, 'g','DisplayName','Y with compensation'); title('Yc[mm]');
    xlabel('Samples');ylabel('[mm]');
    hold;
    plot(Y, 'r','DisplayName','Y without compensation'); title('Y[mm]');
    legend;
    subplot(3,1,3); plot(Ac*180/pi, 'g','DisplayName','A with compensation'); title('Ac [deg]');
    xlabel('Samples');ylabel('Angle');
    hold;
    plot(A*180/pi, 'r','DisplayName','A without compensation'); title('A[deg]');
    legend;

% Plot the estimated variances (in X, Y and A) - 1 standard deviation
subplot(3,1,1); hold on;
    plot(X'+sqrt(P(:,1)), 'g:','HandleVisibility','off');
    plot(X'-sqrt(P(:,1)), 'g:','HandleVisibility','off');
hold off;
subplot(3,1,2); hold on;
    plot(Y'+sqrt(P(:,5)), 'g:','HandleVisibility','off');
    plot(Y'-sqrt(P(:,5)), 'g:','HandleVisibility','off');
hold off;
subplot(3,1,3); hold on;
    plot((A'+sqrt(P(:,9)))*180/pi, 'g:','HandleVisibility','off');
    plot((A'-sqrt(P(:,9)))*180/pi, 'g:','HandleVisibility','off');
hold off;

% Plot the estimated variances (in Xc, Yc and Ac) - 1 standard deviation
subplot(3,1,1); hold on;
    plot(Xc'+sqrt(P(:,1)), 'g:','HandleVisibility','off');
    plot(Xc'-sqrt(P(:,1)), 'g:','HandleVisibility','off');
hold off;
subplot(3,1,2); hold on;
    plot(Yc'+sqrt(P(:,5)), 'g:','HandleVisibility','off');
    plot(Yc'-sqrt(P(:,5)), 'g:','HandleVisibility','off');
hold off;
subplot(3,1,3); hold on;
    plot((Ac'+sqrt(P(:,9)))*180/pi, 'g:','HandleVisibility','off');
    plot((Ac'-sqrt(P(:,9)))*180/pi, 'g:','HandleVisibility','off');
hold off;


