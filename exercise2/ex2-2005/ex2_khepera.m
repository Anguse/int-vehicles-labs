%
% Dead Reckoning with Khepera Mini Robot
%
% Ola Bengtsson, 2004.02.06
%

clear all; close all;

% %%% Khepera settings 1
WHEEL_BASE = 53;                % [mm]
WHEEL_DIAMETER = 15.3;          % [mm]
WHEEL_CIRCUMFERENCE = pi*WHEEL_DIAMETER;    % [mm]
PULSES_PER_REVOLUTION = 600;
MM_PER_PULSE = WHEEL_CIRCUMFERENCE/PULSES_PER_REVOLUTION; % [mm / pulse]

% %%% Uncertainty settings, which are be the same for the left and right encoders
SIGMA_WHEEL_ENCODER = 0.5/12;   % The error in the encoder is 0.5mm / 12mm travelled
% Use the same uncertainty in both of the wheel encoders
SIGMAl = SIGMA_WHEEL_ENCODER;
SIGMAr = SIGMA_WHEEL_ENCODER;

% Load encoder values
ENC = load('khepera.txt');

% Transform encoder values (pulses) into distance travelled by the wheels (mm)
Dr = ENC(1:1:end,2) * MM_PER_PULSE;
Dl = ENC(1:1:end,1) * MM_PER_PULSE;
N = max(size(Dr));

% Init Robot Position, i.e. (0, 0, 90*pi/180) and the Robots Uncertainty
X(1) = 0;
Y(1) = 0;
A(1) = 90*pi/180;
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];
AF(1) = 0;

% Variance of Direction and Angle
sigma2_dD = (SIGMAl+SIGMAr)/4; % [mm]
sigma2_dA = (SIGMAl+SIGMAr)/WHEEL_BASE^2; % [degrees]

% Run until no more encoder values are available
disp('Calculating ...');
for kk=2:N,
    % Change of wheel displacements, i.e displacement of left and right wheels
    dDr = Dr(kk) - Dr(kk-1);
    dDl = Dl(kk) - Dl(kk-1);
    
    % Change of relative movements
    dD = (dDr + dDl)/2;
    dA = (dDr - dDl)/WHEEL_BASE; 
    
    % Adjustment factor
    if(dA == 0)
        AF(kk) = 1;
    else
        AF(kk) = (sin(dA/2))/(dA/2);     
    end 
    
    % Calculate the change in X and Y (World co-ordinates)
    dX = AF(kk)*dD*cos(A(kk-1)+ dA/2);
    dY = AF(kk)*dD*sin(A(kk-1)+dA/2);
    
    % Predict the new state variables (World co-ordinates)
    X(kk) = X(kk-1) + dX;
    Y(kk) = Y(kk-1) + dY;
    A(kk) = mod(A(kk-1) + dA, 2*pi);
    
    % Predict the new uncertainty in the state variables (Error prediction)
    Cxya_old = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];   % Uncertainty in state variables at time k-1 [3x3]    

    Cu =   [sigma2_dD 0;0 sigma2_dA];               % Uncertainty in the input variables [2x2]
    Axya = [1, 0, -dD*sin(A(kk-1)+dA/2);            % Jacobian matrix w.r.t. X, Y and A [3x3]
            0, 1, dD*cos(A(kk-1)+dA/2);
            0, 0, 1];     
    Au =   [cos(A(kk-1)+dA/2), -dD/2*sin(A(kk-1)+dA/2); % Jacobian matrix w.r.t. dD and dA [3x2]
            sin(A(kk-1)+dA/2), (dD/2)*cos(A(kk-1)+dA/2);
            0, 1];           
      
    % Use the law of error predictions, which gives the new uncertainty
    Cxya_new = Axya*Cxya_old*Axya' + Au*Cu*Au';
    
    % Store the new co-variance matrix
    P(kk,:) = [Cxya_new(1,1:3) Cxya_new(2,1:3) Cxya_new(3,1:3)];
    
    % Store all values from dA and dD 
    all_dA(kk) = dA;
    all_dD(kk) = dD;
end

% Khepera settings 2
WHEEL_BASE = 45;                % [mm]
WHEEL_DIAMETER = 14;          % [mm]
WHEEL_CIRCUMFERENCE = pi*WHEEL_DIAMETER;    % [mm]
MM_PER_PULSE = WHEEL_CIRCUMFERENCE/PULSES_PER_REVOLUTION; % [mm / pulse]

% %%% Uncertainty settings, which are be the same for the left and right encoders
SIGMA_WHEEL_ENCODER = 0.5/12;   % The error in the encoder is 0.5mm / 12mm travelled
SIGMA_WHEEL_BASE = 4; % The error in wheelbase is 13mm

% Use the same uncertainty in both of the wheel encoders
kl = 1.3*pi/WHEEL_CIRCUMFERENCE   % The error on each wheel
kr = 1.3*pi/WHEEL_CIRCUMFERENCE
kb = SIGMA_WHEEL_BASE;

% Load encoder values
ENC = load('khepera.txt');


% Transform encoder values (pulses) into distance travelled by the wheels (mm)
Dr = ENC(1:1:end,2) * MM_PER_PULSE;
Dl = ENC(1:1:end,1) * MM_PER_PULSE;
N = max(size(Dr));

% Init Robot Position, i.e. (0, 0, 90*pi/180) and the Robots Uncertainty
X2(1) = 0;
Y2(1) = 0;
A2(1) = 90*pi/180;
P2(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];
AF(1) = 0;

syms x y t dSr dSl b
% Function for the next state variables
f = [x + (dSr+dSl)/2*cos(t+(dSr-dSl)/(2*b));y + (dSr+dSl)/2*sin(t+(dSr-dSl)/(2*b)); t+(dSr-dSl)/b] % (5.7)
% Calculate the jacobian with respect to state, distance traveled with left and right wheel and wheelbase
v1 = [x y t]
v2 = [dSr dSl]
v3 = b;
Rs = jacobian(f, v1) 
Rm = jacobian(f, v2)
Rb = jacobian(f, v3);

% Variance of Direction and Angle
sigma2_dD = (SIGMAl+SIGMAr)/4; % [mm]
sigma2_dA = (SIGMAl+SIGMAr)/WHEEL_BASE^2; % [degrees]

% Run until no more encoder values are available
disp('Calculating ...');

for kk=2:N,
    
    % Change of wheel displacements, i.e displacement of left and right
    % wheels
    dDr = Dr(kk) - Dr(kk-1);
    dDl = Dl(kk) - Dl(kk-1);
    
    % Change of relative movements
    dD = (dDr + dDl)/2;
    dA = (dDr - dDl)/WHEEL_BASE; 
    
    % Predict the new state variables (World co-ordinates)
    Xn = eval(subs(f,[x, y, t, dSr, dSl, b],[X(kk-1), Y(kk-1), A(kk-1), dDr, dDl, WHEEL_BASE]));
    X2(kk) = Xn(1);
    Y2(kk) = Xn(2);
    A2(kk) = Xn(3);
    
    % Jacobians for state variables and measurement variables
    %Js = eval(subs(Rs,[t,dSr,dSl,b],[A(kk-1), dDr, dDl, WHEEL_BASE]));
    %Jm = eval(subs(Rm,[t,dSr,dSl,b],[A(kk-1), dDr, dDl, WHEEL_BASE]));
    %Jb = eval(subs(Rb,[t,dSr,dSl,b],[A(kk-1), dDr, dDl, WHEEL_BASE]));
    
    Js = [ 1, 0, -sin(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl/2 + dDr/2)
          0, 1,  cos(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl/2 + dDr/2)
          0, 0,                                           1];
    Jm = [ cos(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))/2 - (sin(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl/2 + dDr/2))/(2*WHEEL_BASE), cos(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))/2 + (sin(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl/2 + dDr/2))/(2*WHEEL_BASE)
          sin(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))/2 + (cos(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl/2 + dDr/2))/(2*WHEEL_BASE), sin(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))/2 - (cos(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl/2 + dDr/2))/(2*WHEEL_BASE)
                                                                                        1/WHEEL_BASE,                                                                              -1/WHEEL_BASE];
    Jb =  [-(sin(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl - dDr)*(dDl/2 + dDr/2))/(2*WHEEL_BASE^2);
            (cos(A(kk-1) - (dDl - dDr)/(2*WHEEL_BASE))*(dDl - dDr)*(dDl/2 + dDr/2))/(2*WHEEL_BASE^2);
                                                   (dDl - dDr)/WHEEL_BASE^2]
    
    % Old state uncertainty and measurement uncertainty
    Sx = [P2(kk-1,1:3);P2(kk-1,4:6);P2(kk-1,7:9)];   % Uncertainty in state variables at time k-1 [3x3]
    Su = [sigma2_dD 0; 0 sigma2_dA];
    Sm = [kr*abs(dDr) 0;0 kl*abs(dDl)];
    Sb = kb;
    
    % New state uncertainty
    Sxn = Js*Sx*Js' + Jm*Sm*Jm' + Jb*Sb*Jb';
    P2(kk,1:9) =[Sxn(1,1:3) Sxn(2,1:3) Sxn(3,1:3)];
    
    % Store all values from dA and dD 
    all_dA(kk) = dA;
    all_dD(kk) = dD;
end

figure, plot(all_dA, 'DisplayName', '\delta\theta');
hold;
plot(AF, 'DisplayName', 'Adjustment factor');
legend;

figure, subplot(2,1,1);
plot(all_dA, 'DisplayName', 'dA');
legend;
subplot(2,1,2);
plot(all_dD, 'DisplayName', 'dD');
legend;

figure, subplot(3,1,1);
plot(X, 'DisplayName', 'X');
legend;
subplot(3,1,2);
plot(Y, 'DisplayName', 'Y');
subplot(3,1,3);
plot(A, 'DisplayName', 'A');
legend;
%%

disp('Plotting ...');

% Plot the path taken by the robot, by plotting the uncertainty in the current position
figure; 
    %plot(X, Y, 'b.');
    title('Path taken by the robot [Wang]');
    xlabel('X [mm] World co-ordinates');
    ylabel('Y [mm] World co-ordinates');
    hold on;
        for kk = 1:1:N,
            C = [P(kk,1:3);P(kk,4:6);P(kk,7:9)];
            C2 = [P2(kk,1:3);P2(kk,4:6);P2(kk,7:9)];
            plot_uncertainty([X(kk) Y(kk) A(kk)]', C, 1, 2, 'r');
            plot_uncertainty([X2(kk) Y2(kk) A2(kk)]', C2, 1, 2, 'g');
        end;
        legend("With measurement uncertainty", "With uncertainty in wheelbase and wheel diameter")
    hold off;
    axis('equal');

% After the run, plot the results (X,Y,A), i.e. the estimated positions 
figure;
    subplot(3,1,1); plot(X, 'r'); title('X');
    hold;
    plot(X2, 'g');
    legend("With measurement uncertainty", "With uncertainty in wheelbase and wheel diameter");
    xlabel('Samples');
    ylabel('[mm]');
    subplot(3,1,2); plot(Y, 'r'); title('Y');
    hold;
    plot(Y2, 'g');
    legend("With measurement uncertainty", "With uncertainty in wheelbase and wheel diameter");
    xlabel('Samples');
    ylabel('[mm]');
    subplot(3,1,3); plot(A*180/pi, 'r'); title('A');
    hold;
    plot(A2*180/pi, 'g');
    legend("With measurement uncertainty", "With uncertainty in wheelbase and wheel diameter");
    xlabel('Samples');
    ylabel("Angle [°]");

% Plot the estimated variances (in X, Y and A) - 1 standard deviation
subplot(3,1,1); hold on;
    plot(X'+sqrt(P(:,1)), 'r:', 'HandleVisibility', 'off');
    plot(X'-sqrt(P(:,1)), 'r:', 'HandleVisibility', 'off');
    plot(X2'+sqrt(P2(:,1)), 'g:', 'HandleVisibility', 'off');
    plot(X2'-sqrt(P2(:,1)), 'g:', 'HandleVisibility', 'off');
hold off;
subplot(3,1,2); hold on;
    plot(Y'+sqrt(P(:,5)), 'r:', 'HandleVisibility', 'off');
    plot(Y'-sqrt(P(:,5)), 'r:', 'HandleVisibility', 'off');
    plot(Y2'+sqrt(P2(:,5)), 'g:', 'HandleVisibility', 'off');
    plot(Y2'-sqrt(P2(:,5)), 'g:', 'HandleVisibility', 'off');
hold off;
subplot(3,1,3); hold on;
    plot((A'+sqrt(P(:,9)))*180/pi, 'r:', 'HandleVisibility', 'off');
    plot((A'-sqrt(P(:,9)))*180/pi, 'r:', 'HandleVisibility', 'off');
    plot(A2'+sqrt(P2(:,9)), 'g:', 'HandleVisibility', 'off');
    plot(A2'-sqrt(P2(:,9)), 'g:', 'HandleVisibility', 'off');
hold off;
