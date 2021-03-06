%% Snowhite robot data
% Summary of example objective

%% Section 1 : Setup
clear all; close all;
L = 680; % Wheel base[mm]
T = 0.05 % Sample time [s]
data = load("snowhite.txt");

% Load input data
v = data(:,1);
a = data(:,2);
true_X = data(:,3);
true_Y = data(:,4);
true_A = data(:,5);
N = max(size(data));

% Initial position and uncertainty
X(1) = true_X(1);
Y(1) = true_Y(1);
A(1) = true_A(1);
P(1,1:9) = [1 0 0 0 1 0 0 0 (1*pi/180)^2];

% Estimated measurement uncertainty
SIGMAV = var(v);
SIGMAa = var(a);
SIGMAT = 0.00001;

%% Computation
for kk=2:N,
    
    dD = v(kk)*cos(a(kk))*T
    dA = (v(kk)*sin(a(kk))*T)/L
    
    % Calculate the change in X and Y (World co-ordinates)
    dX = dD*cos(A(kk-1)+ dA/2);
    dY = dD*sin(A(kk-1)+dA/2);
    
    % Predict the new state variables (World co-ordinates)
    X(kk) = X(kk-1) + dX;
    Y(kk) = Y(kk-1) + dY;
    A(kk) = mod(A(kk-1) + dA, 2*pi);
    
    % Old state uncertainty and measurement uncertainty
    Sx = [P(kk-1,1:3);P(kk-1,4:6);P(kk-1,7:9)];   % Uncertainty in state variables at time k-1 [3x3]
    Su = [SIGMAV 0 0;
          0 SIGMAa 0;
          0 0 SIGMAT];
    
    % Jacobians for state variables and measurement variables
    Js = [ 1, 0, -T*v(kk)*cos(a(kk))*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L);
          0, 1,  T*v(kk)*cos(a(kk))*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L);
          0, 0,                                        1];
    % w.r.t [v,a,T]
    Ju = [ T*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk)) - (T^2*v(kk)*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk))*sin(a(kk)))/L, - T*v(kk)*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*sin(a(kk)) - (T^2*v(kk)^2*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk))^2)/L, v(kk)*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk)) - (T*v(kk)^2*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk))*sin(a(kk)))/L
          T*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk)) + (T^2*v(kk)*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk))*sin(a(kk)))/L,   (T^2*v(kk)^2*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk))^2)/L - T*v(kk)*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*sin(a(kk)), v(kk)*sin(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk)) + (T*v(kk)^2*cos(A(kk-1) + (T*v(kk)*sin(a(kk)))/L)*cos(a(kk))*sin(a(kk)))/L
                                                                                (T*sin(a(kk)))/L,                                                                      (T*v(kk)*cos(a(kk)))/L,                                                                       (v(kk)*sin(a(kk)))/L];
                                                                                     
    Sxn = Js*Sx*Js' + Ju*Su*Ju';
    P(kk,1:9) =[Sxn(1,1:3) Sxn(2,1:3) Sxn(3,1:3)];
    
    
    
    % Actual error in current state
    Es(kk,1) = abs(X(kk)-true_X(kk));
    Es(kk,2) = abs(Y(kk)-true_Y(kk));
    Es(kk,3) = abs(A(kk)-true_A(kk));
    
    % Estimated error in current state
    eEs(kk,1) = sqrt(Sxn(1,1));
    eEs(kk,2) = sqrt(Sxn(2,2));
    eEs(kk,3) = sqrt(Sxn(3,3));
    
    % Store all values from dA and dD 
    all_dA(kk) = dA;
    all_dD(kk) = dD;
    all_dX(kk) = dX;
    all_dY(kk) = dY;
    
end;

%% Section 2 Plot
% Plot the path taken by the robot, by plotting the uncertainty in the current position
figure; 
    %plot(X, Y, 'b.');
    title('Path taken by the robot [Wang]');
    xlabel('X [mm] World co-ordinates');
    ylabel('Y [mm] World co-ordinates');
    plot(true_X, true_Y, 'r', 'LineWidth', 5);
    hold on;
        for kk = 1:10:N,
            C = [P(kk,1:3);P(kk,4:6);P(kk,7:9)];
            plot_uncertainty([X(kk) Y(kk) A(kk)]', C, 1, 2, 'g');
        end;
        legend("True value", "Estimated value with variance")
        xlabel('X[mm]');
        ylabel('Y[mm]');
    hold off;
    axis('equal');

figure, plot(true_X, true_Y, 'r');
hold;
plot(X,Y, 'g');
legend('True value', 'Estimated value');
xlabel('X[mm]');
ylabel('Y[mm]');
    
% After the run, plot the results (X,Y,A), i.e. the estimated positions 
figure;
    subplot(3,1,1); plot(X, 'g'); title('X');
    hold;
    plot(true_X, 'r');
    legend("True value", "Estimated value");
    xlabel('Samples');
    ylabel('[mm]');
    subplot(3,1,2); plot(Y, 'g'); title('Y');
    hold;
    plot(true_Y, 'r');
    legend("True value", "Estimated value");
    xlabel('Samples');
    ylabel('[mm]');
    subplot(3,1,3); plot(A*180/pi, 'g'); title('A');
    hold;
    plot(true_A*180/pi, 'r');
    legend("True value", "Estimated value");
    xlabel('Samples');
    ylabel("Angle [°]");

% Plot the estimated variances (in X, Y and A) - 1 standard deviation
subplot(3,1,1); hold on;
    plot(X'+sqrt(P(:,1)), 'r:', 'HandleVisibility', 'off');
    plot(X'-sqrt(P(:,1)), 'r:', 'HandleVisibility', 'off');
hold off;
subplot(3,1,2); hold on;
    plot(Y'+sqrt(P(:,5)), 'r:', 'HandleVisibility', 'off');
    plot(Y'-sqrt(P(:,5)), 'r:', 'HandleVisibility', 'off');
hold off;
subplot(3,1,3); hold on;
    plot((A'+sqrt(P(:,9)))*180/pi, 'r:', 'HandleVisibility', 'off');
    plot((A'-sqrt(P(:,9)))*180/pi, 'r:', 'HandleVisibility', 'off');
hold off;

    figure, subplot(3,1,1);
    plot(Es(:,1), 'r')
    hold;
    plot(eEs(:,1), 'g')
    title('X error')
    legend('True error', 'Estimated error')
    xlabel('Sample')
    ylabel('Distance [mm]')
    subplot(3,1,2);
    plot(Es(:,2), 'r');
    hold;
    plot(eEs(:,2), 'g')
    title('Y error')
    legend('True error', 'Estimated error')
    xlabel('Sample')
    ylabel('Distance [mm]')
    subplot(3,1,3);
    plot(Es(:,3), 'r');
    hold;
    plot(eEs(:,3), 'g')
    title('A error')
    legend('True error', 'Estimated error')
    xlabel('Sample')
    ylabel('°')
    