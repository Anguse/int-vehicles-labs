%% Example Title
% Summary of example objective
clear all; close all;
%% Section 1 Title
% Description of first code block

% Co-ordinates of the ref. nodes
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
         12 16;     % L3.875
         16 15;     % L12
         15 14;     % L13
         19 21;     % L14
         22 18;     % L15
         20 23];    % L16
     
LINEMODEL = [REF(LINES(:,1),1:2) REF(LINES(:,2),1:2)];

figure, hold on;
rot = [0 -1;1 0] % 90 degree rotation matrix
for kk = 1:size(LINEMODEL, 1),
    L1 = [LINEMODEL(kk,1:2); LINEMODEL(kk,3:4)]; % [X1, X2; Y1, Y2]
    L2 = [LINEMODEL(kk,1), LINEMODEL(kk,1); LINEMODEL(kk,3), LINEMODEL(kk,3)]; % [X1, X1; Y1, Y1]
    V = rot*(L1-L2);      % Rotate line 90 degrees and move it to origo (rot*[0,  X2-X1; 0, Y2-Y1])
    U = V/norm(V);        % Normalize for unitvector
    R(kk,1:2) = dot(U, L1);
    line([U(1,:)], [U(2,:)])
end;
hold off;



%% Section 2 Title
% Description of second code block
b = 2;
