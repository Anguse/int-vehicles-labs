function [dx, dy, da, C] = Cox_LineFit(angs, meas, p_state, sensor_offset, LINEMODEL)

% Setup
alpha = sensor_offset(1);  % [Scalar] Offset between sensors centre and robots x-axis
beta = sensor_offset(2);   % [Scalar] Offset between sensors centre and robots y-axis
gamma = sensor_offset(3);  % [Scalar] Heading offset between sensors and robots x-axis(-pi/2 radians)
vm = mean([meas.*cos(angs) meas.*sin(angs)]); % Mean of all points
dx = 0;
dy = 0;
da = 0;

% Get unit vectors of all lines
rot = [0 -1;1 0]; % 90 degree rotation matrix
for kk = 1:size(LINEMODEL, 1)
    L1 = [LINEMODEL(kk,1) LINEMODEL(kk,3); LINEMODEL(kk,2) LINEMODEL(kk,4)]; % [X1, Y1; X2, Y2]
    L2 = [LINEMODEL(kk,1), LINEMODEL(kk,1);      % [X1, Y1; X1, Y1]
          LINEMODEL(kk,2), LINEMODEL(kk,2)]; 
    V = rot*(L1-L2);       % rot*[0,  0; X2-X1, Y2-Y1] Rotate line 90 degrees and move it to origo 
    Ui = V/norm(V);         % Normalize for unitvector
    pU(kk,1:2) = Ui(1,1:2);
    pU(kk,3:4) = Ui(2,1:2);
    Ui = [Ui(3) Ui(4)];     % [Xu Yu]
    z = [L1(3) L1(4)];     % A point on L1
    Ri = dot(Ui, z);       % Project unitvector on L1
    U(kk,1:2) = Ui;        
    R(kk,1) = Ri;
end

% Get meadian value of all target distances
for ii = 1:size(angs, 2)
    x = meas(kk)*cos(angs(kk));
    y = meas(kk)*sin(angs(kk));
    rot = [cos(gamma) -sin(gamma) alpha; sin(gamma) cos(gamma) beta; 0 0 1];
    Xs = rot*[x y 1]';
    rot = [cos(p_state(3)) -sin(p_state(3)) p_state(1); sin(p_state(3)) cos(p_state(3)) p_state(2); 0 0 1];
    Xw = rot*[Xs(1) Xs(2) 1]';
    v = [Xw(1), Xw(2)]';
    for kk = 1:size(U, 1)
        yy(kk) = abs(R(kk) - dot(U(kk,:), v))
    end
    [targ(ii), targind] = min(yy)
end
target_med = median(targ);
threshold = target_med;


%for kj = 1:size(pU, 1)
%    line(pU(kj,1:2), pU(kj,3:4));
%end

% The loop
finished = false;
while ~finished  
    
    for kk = 1:size(angs, 2)
        
        % 1.Translate and rotate data points
        x = meas(kk)*cos(angs(kk));
        y = meas(kk)*sin(angs(kk));
        rot = [cos(gamma) -sin(gamma) alpha; sin(gamma) cos(gamma) beta; 0 0 1];
        Xs = rot*[x y 1]';
        rot = [cos(p_state(3)) -sin(p_state(3)) p_state(1); sin(p_state(3)) cos(p_state(3)) p_state(2); 0 0 1];
        Xw = rot*[Xs(1) Xs(2) 1]';
        v = [Xw(1), Xw(2)]';
        
        % 2.Find the target line
        for ii = 1:size(U, 1)
            yi(ii) = abs(R(ii) - dot(U(ii,:), v));
        end
        [target, target_index] = min(yi);
        if target^2 < threshold % Reject all outliers
            X1(kk,1) = U(target_index, 1);
            X2(kk,1) = U(target_index, 2);
            X3(kk,1) = U(target_index, :)*[0 -1;1 0]*(v-vm);
            targets(kk) = target;
        end
    end
    
    targets = nonzeros(targets);
    if size(targets,1) == 0
        if threshold > 300
            dx = 0
            dy = 0
            da = 0
            C = [1 0 0;
                 0 1 0;
                 0 0 1]
            return;
        end
        threshold = threshold + 1
        continue;
    end
    % 3.Set up linear equation system
    A = [X1 X2 X3];
    A = A(any(A,2),:) % Remove empty rows

    B =  A\targets;%inv(A'*A)*A'*target;A\targets;

    S2 = ((targets-A*B)'*(targets-A*B))/(max(size(A))-4);
    C = S2*(inv(A'*A));

    % 4.Add latest contribution to the overall congruence
    dx = B(1);
    dy = B(2);
    da = B(3);
    
    % Change the position
    p_state(1) = p_state(1)+dx;
    p_state(2) = p_state(2)+dy;
    p_state(3) = mod(p_state(3)+da, 2*pi);

    % 5.Check if process has converged
    if(sqrt(dx^2+dy^2+da^2) < 5)&&(abs(da<0.1*pi/180))
        finished = true;
        break;
    end
   
    clear X1;
    clear X2;
    clear X3;
    targets = 0;
end

dx
dy
da

