% Intelligent vehicles
% Exercise 2 part 2:5-6
clear all, close all;
warning ('off', 'MATLAB:plot:IgnoreImaginaryXYPart');

% %%% Khepera settings 
WHEEL_BASE = 53;                % 53 [mm]
b = WHEEL_BASE;
WHEEL_DIAMETER = 15.3;          % 15.3 [mm]
PULSES_PER_REVOLUTION = 600;    % 600
WHEEL_CIRCUMFERENCE = 15.3*pi;
MM_PER_PULSE = WHEEL_CIRCUMFERENCE/PULSES_PER_REVOLUTION; % [mm / pulse]

% %%% Uncertainty settings, which are be the same for the left and right encoders
SIGMA_WHEEL_ENCODER = 0.5/12;   % The error in the encoder is 0.5mm / 12mm travelled
% Use the same uncertainty in both of the wheel encoders
SIGMAl = SIGMA_WHEEL_ENCODER;
SIGMAr = SIGMA_WHEEL_ENCODER;
% uncertainty of the wheelbase
SIGMAb = 4;

% Load encoder values
%ENC = load('khepera_circle.txt');
ENC = load('khepera_circle.txt');

% Transform encoder values (pulses) into distance travelled by the wheels (mm)
Dr = ENC(1:1:end,2) * MM_PER_PULSE; %Change here not to read all (1:10:end,1)
Dl = ENC(1:1:end,1) * MM_PER_PULSE;
N = max(size(Dr));

% Init Robot Position, i.e. (0, 0, 90*pi/180) and the Robots Uncertainty
X(1) = 0;
Y(1) = 0;
A(1) = 90*pi/180;
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];

% Run until no more encoder values are available
disp('Calculating ...');
for kk=2:N
    % Change of wheel displacements, i.e displacement of left and right wheels
    dDr = Dr(kk) - Dr(kk-1);
    dDl = Dl(kk) - Dl(kk-1);
    
    % Change of relative movements
    dD = (dDr+dDl)/2;
    dA = (dDr-dDl)/b;
    
    % Calculate the change in X and Y (World co-ordinates)
    %dX = dD*cos(A(kk-1)+(dA/2));
    %dY = dD*sin(A(kk-1)+(dA/2));
   
    % Predict the new state variables (World co-ordinates)
    X(kk) = X(kk-1) + (dDr+dDl)*0.5*cos(A(kk-1)+(dDr-dDl)/(2*b));
    Y(kk) = Y(kk-1) + (dDr+dDl)*0.5*sin(A(kk-1)+(dDr-dDl)/(2*b));
    A(kk) = mod(A(kk-1) + (dDr-dDl)/b, 2*pi);
    
    % Part with the compensation
    %X(kk) = X(kk-1) + sinc(dA/2)*dX;
    %Y(kk) = Y(kk-1) + sinc(dA/2)*dY;
    %A(kk) = mod(A(kk-1) + dA, 2*pi);
    
    % Predict the new uncertainty in the state variables (Error prediction)
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];   % Uncertainty in state variables at time k-1 [3x3]    

    Jxya = [1 0 -dD*sin(A(kk-1)+dA/2);
            0 1 dD*cos(A(kk-1)+dA/2);
            0 0 1];
    
    Jrl = [0.5*cos(A(kk-1)+dA/2)-dD*0.5*sin(A(kk-1)+dA/2)/b 0.5*cos(A(kk-1)+dA/2)+dD*0.5*sin(A(kk-1)+dA/2)/b;
            0.5*sin(A(kk-1)+dA/2)+dD*0.5*cos(A(kk-1)+dA/2)/b 0.5*sin(A(kk-1)+dA/2)-dD*0.5*cos(A(kk-1)+dA/2)/b;
            1/b 1/b];
    
    kr = 0.001;
    kl = 0.001;
        
    Crl = [kr*abs(dDr) 0;
            0 kl*abs(dDl)];
    
    Jb = [dD*((dDr-dDl)/(2*b^2))*sin(A(kk-1)+dA/2);
            -dD*((dDr-dDl)/(2*b^2))*cos(A(kk-1)+dA/2);
            -(dDr-dDl)/b^2];
    
    Cb = [SIGMAb];
        
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Jxya*Cxya_old*Jxya' + Jrl*Crl*Jrl' + Jb*Cb*Jb';
    
    % Store the new co-variance matrix
    P(kk,1:9) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];
end


disp('Plotting ...');

% Plot the path taken by the robot, by plotting the uncertainty in the current position
figure; 
    %plot(X, Y, 'b.');
    title('Path taken by the robot [Wang]');
    xlabel('X [mm] World co-ordinates');
    ylabel('Y [mm] World co-ordinates');
    hold on;
        for kk = 1:10:N % to plot every 10:th write 1:10:N
            C = [P(kk,1:3);P(kk,4:6);P(kk,7:9)];
            plot_uncertainty([X(kk) Y(kk) A(kk)]', C, 1, 2);
            %plot_khepera([X(kk) Y(kk) A(kk)]', 10, 10, 2);
        end
    hold off;
    axis('equal');

% After the run, plot the results (X,Y,A), i.e. the estimated positions 
figure;
    subplot(3,1,1); plot(X, 'b'); title('X [mm]');
    subplot(3,1,2); plot(Y, 'b'); title('Y [mm]');
    subplot(3,1,3); plot(A*180/pi, 'b'); title('A [deg]');

% Plot the estimated variances (in X, Y and A) - 1 standard deviation
subplot(3,1,1); hold on;
    plot(X'+sqrt(P(:,1)), 'b:');
    plot(X'-sqrt(P(:,1)), 'b:');
hold off;
subplot(3,1,2); hold on;
    plot(Y'+sqrt(P(:,5)), 'b:');
    plot(Y'-sqrt(P(:,5)), 'b:');
hold off;
subplot(3,1,3); hold on;
    plot((A'+sqrt(P(:,9)))*180/pi, 'b:');
    plot((A'-sqrt(P(:,9)))*180/pi, 'b:');
hold off;
