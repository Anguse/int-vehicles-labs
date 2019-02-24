function [dx, dy, da, C] = Cox_LineFit(angs, meas, p_state, sensor_offset, LINEMODEL)
% Setup
alpha = sensor_offset(1);  % [Scalar] Offset between sensors centre and robots x-axis
beta = sensor_offset(2);   % [Scalar] Offset between sensors centre and robots y-axis
gamma = sensor_offset(3);  % [Scalar] Heading offset between sensors and robots x-axis(-pi/2 radians)

% Get unit vectors of all lines
rot = [0 -1;1 0] % 90 degree rotation matrix
for kk = 1:size(LINEMODEL, 1)
    L1 = [LINEMODEL(kk,1:2); LINEMODEL(kk,3:4)]; % [X1, X2; Y1, Y2]
    L2 = [LINEMODEL(kk,1), LINEMODEL(kk,1); LINEMODEL(kk,3), LINEMODEL(kk,3)]; % [X1, X1; Y1, Y1]
    V = rot*(L1-L2);      % Rotate line 90 degrees and move it to origo (rot*[0,  X2-X1; 0, Y2-Y1])
    Ui = V/norm(V);        % Normalize for unitvector
    Ri = dot(Ui, L1);
    U(kk,1:2) = Ui(1,:);
    U(kk,3:4) = Ui(2,:);
    R(kk,1:2) = Ri;
end;

U(1,:)

% The loop
for kk = 1:size(angs, 1)
    
    % 1.Translate and rotate data points
    x = meas(kk)*cos(angs(kk));
    y = meas(kk)*sin(angs(kk));
    R = [cos(gamma) -sin(gamma) alpha; sin(gamma) cos(gamma) beta; 0 0 1]
    Xs = R*[x y 1]'
    R = [cos(p_state(3)) -sin(p_state(3)) p_state(1); sin(p_state(3)) cos(p_state(3)) p_state(2); 0 0 1]
    Xw = R*[Xs(1) Xs(2) 1]'
    
    % 2.Find the target for data points
    vi = [Xs(1), Xs(2)];
    for ii = 1:size(LINEMODEL, 1)
        yi(ii) = R(ii,:) - dot([U(ii,1:2); U(ii,3:4)], vi)
    end
    y = min(yi);
    
    % 3.Set up linear equation system
    
end

dx = 1;
dy = 2;
da = 3;
C = 4;

