function [ddx, ddy, dda, C] = Cox_LineFit(angs, meas, p_state, sensor_offset, LINEMODEL)

% Setup
alpha = sensor_offset(1);  % [Scalar] Offset between sensors centre and robots x-axis
beta = sensor_offset(2);   % [Scalar] Offset between sensors centre and robots y-axis
gamma = sensor_offset(3);  % [Scalar] Heading offset between sensors and robots x-axis(-pi/2 radians)
vm = mean([meas.*cos(angs) meas.*sin(angs)]) % Mean of all points
ddx = 0
ddy = 0
dda = 0

% Get unit vectors of all lines
rot = [0 -1;1 0] % 90 degree rotation matrix
for kk = 1:size(LINEMODEL, 1)
    L1 = [LINEMODEL(kk,1:2); LINEMODEL(kk,3:4)]; % [X1, X2; Y1, Y2]
    L2 = [LINEMODEL(kk,1), LINEMODEL(kk,1);      % [X1, X1; Y1, Y1]
          LINEMODEL(kk,3), LINEMODEL(kk,3)]; 
    V = rot*(L1-L2);       % rot*[0,  X2-X1; 0, Y2-Y1] Rotate line 90 degrees and move it to origo 
    Ui = V/norm(V)         % Normalize for unitvector
    Ri = dot(Ui, L1);      % Project unitvector on L1
    U(kk,1:2) = [Ui(3) Ui(4)]
    R(kk,1:2) = Ri
end;


% The loop
for kk = 1:size(angs, 1)
    % 1.Translate and rotate data points
    x = meas(kk)*cos(angs(kk));
    y = meas(kk)*sin(angs(kk));
    rot = [cos(gamma) -sin(gamma) alpha; sin(gamma) cos(gamma) beta; 0 0 1]
    Xs = rot*[x y 1]'
    rot = [cos(p_state(3)) -sin(p_state(3)) p_state(1); sin(p_state(3)) cos(p_state(3)) p_state(2); 0 0 1]
    Xw = rot*[Xs(1) Xs(2) 1]'
    
    % 2.Find the target for data points
    vi = [Xw(1), Xw(2)]'
    for ii = 1:size(LINEMODEL, 1)
        yi(ii) = R(ii,2) - dot(U(ii,:), vi)
    end
    y = min(yi);
    
    % 3.Set up linear equation system
    X1 = U(ii,1)
    X2 = U(ii,2)
    X3 = U(ii,:)*[0 -1;1 0]*(vi-vm)
    
    A = [X1 X2 X3]
    B = inv(A'*A)*A'*y;
    S2 = ((y-A*B)^2)/(max(size(A))-4);
    
    % 4.Add latest contribution to the overall congruence
    dx = B(1)
    dy = B(2)
    da = B(3)
    
    ddx = ddx + dx;
    ddy = ddy + dy;
    dda = dda + da;
    
    if(sqrt(dx^2+dy^2+da^2) < 5)&&(abs(da<0.1*pi/180))
        break;
    else
        continue;
    end
    
end

