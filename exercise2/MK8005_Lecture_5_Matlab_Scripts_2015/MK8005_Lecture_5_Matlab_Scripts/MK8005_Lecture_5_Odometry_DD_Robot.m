% (C), Bj�rn �strand 2013, 2015
% To be used in Intelligent Vehicles, Odometry example
% 
close all; clear all;
% definde symbolic values
% http://www.mathworks.se/help/symbolic/jacobian.html
% http://www.mathworks.se/help/symbolic/creating-symbolic-variables-and-expressions.html
% http://www.mathworks.se/help/symbolic/performing-symbolic-computations.html

% Init
X = [0 0 0]';
Sx = [0.1 0 0;0 0.1 0;0 0 0]'; % Starting state unsertainty
% Uncertainty parameters (symbolic values)
syms x y t dSr dSl b
% define the function
f = [x + (dSr+dSl)/2*cos(t+(dSr-dSl)/(2*b));y + (dSr+dSl)/2*sin(t+(dSr-dSl)/(2*b)); t+(dSr-dSl)/b] % (5.7)
% Calculate the jacobian with respect to STATE
v1 = [x y t]
% claculate the Jacobian R = jacobian(f, v)
Rs = jacobian(f, v1) 
% Calculate the jacobian with respect to MEASUREMENT
v2 = [dSr dSl]
% claculate the Jacobian R = jacobian(f, v)
Rm = jacobian(f, v2)
v3 = [b];
Rb = jacobian(f, v3);

figure; hold on;
plot(X(1), X(2));
for p=1:60,
    % New state (moving constant speed d = 0.1 and turning 3 degrees)
    % http://www.mathworks.se/help/symbolic/subs.html
    %x = X(1); y = X(2); t = X(3); d = 0.01; a = 3/180*pi;
    ddSr = 1.05; % equal = drive straight
    ddSl = 1.00;
    wb = 0.5;
    Xn = eval(subs(f,[x, y, t, dSr, dSl, b],[X(1), X(2), X(3), ddSr, ddSl,wb]));
    %Xn = double(Xn)
    % New state unsertainty
    Js = eval(subs(Rs,[t,dSr,dSl,b],[X(3), ddSr, ddSl,wb]));
    Jm = eval(subs(Rm,[t,dSr,dSl,b],[X(3), ddSr, ddSl,wb]));
    Jb = eval(subs(Rb,[t,dSr,dSl,b],[X(3), ddSr, ddSl,wb]));
    % Measurement covariance
    kr = 0.01^2;
    kl = 0.01^2;
    kb = 0.01;
    Sm = [kr*abs(ddSr) 0;0 kl*abs(ddSl)];
    Sb = kb;
    % New state uncertainty
    Sxn = Js*Sx*Js' + Jm*Sm*Jm' + Jb*Sb*Jb';
    % print
    plot(Xn(1), Xn(2));
    plot_uncertainty([Xn(1) Xn(2)]', Sxn, 1, 2);
    hold on;
    axis equal;
    grid on;
    %
    X = Xn; Sx = Sxn;
    drawnow();
end


