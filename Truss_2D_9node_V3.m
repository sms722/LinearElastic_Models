%% 2D 9-node problem

clc;    
close all; 
clear;

sel = 1; % unit cube side length of 1 m/ truss length of 0.5 m
A = 0.0005; % Cross-sectional area of 0.0005 m^2 / 5 cm^2
r = sqrt(A/pi); % Radius for (assumed circular) member c/s area 
                % listed above 
E = 10000; % Young's Modulus of 10 kPa (polymeric material)

NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1];
CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];

% Forming Elemental Stiffness Matrices
Kbasket = [];
for i = 1:size(CA,1)
    x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
    y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
    L = sqrt(((x2-x1)^2)+((y2-y1)^2));
    c=(x2-x1)/L; c2 = c^2;
    s=(y2-y1)/L; s2 = s^2;
    ktemp = [c2,   c*s,  -c2,  -c*s;  
             c*s,   s2,  -c*s,  -s2;    
             -c2, -c*s,  c2,    c*s; 
             -c*s, -s2,  c*s,   s2];
    ke = ((A.*E)./L).*ktemp;
    Kbasket(:,:,i) = ke;
end

% Global-to-local-coordinate-system Coordination
GlobToLoc=zeros(size(CA,1),4);
for n=1:2  
    GN=CA(:,n); 
    for d=1:2
        GlobToLoc(:,(n-1)*2+d)=(GN-1)*2+d;
    end
end

% Forming Global Truss Stiffness Matrix
K = zeros(2*size(NC,1));
for e=1:size(CA,1) 
    ke = Kbasket(:,:,e);
    for lr = 1:4
        gr = GlobToLoc(e,lr); 
        for lc = 1:4
            gc = GlobToLoc(e,lc); 
            K(gr,gc) = K(gr,gc) + ke(lr,lc);
        end
    end
end

% Solving for C matrix components (master pseudocode):
% Define strain vector: [e11, e22, e12]'
% for each strainvec component:
%   set that component equal to a dummy value, set all other values to zero
%   use strain relations and BCs to develop all displacement values, except
%       for center node
%   develop force vector for all nodes (except center node) using sized-down K
%       matrix
%   use force/displacement relations to solve at center node
%   use system-wide F matrix and force relations to solve for stress vector
%       [s11,s22,s12]'
%   use strain and stress vectors to solve for the corresponding row of the
%       C matrix

C = [];
for i = 1:3
%   Define strain vector: [e11, e22, e12]'
    strainvec = [0;0;0];
    
%   set that component equal to a dummy value, set all other values to zero
    strainvec(i) = 0.01; 
    
%   use strain relations and BCs to develop all displacement values, except
%       for center node
    e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
    u = zeros(16,1);
    u(5) = e12*sel; u(6) = e11*sel;
    u(15) = sel.*(e11+e12); u(16) = sel.*(e12+e22);
    u(11) = e11*sel;
    if e11 ~= 0 % when testing non-zero e11
        u(3) = 0; u(13) = e11*sel;
        u(4) = 0; u(14) = 0;
        u(8) = 0; u(10) = 0;
        u(7) = u(9); % INCOMPLETE
    elseif e22 ~= 0 % when testing non-zero e22
        u(3) = 0; u(13) = 0;
        u(7) = 0; u(9) = 0;
        u(8) = 0; u(10) = e22*sel;
        u(4) = u(14); % INCOMPLETE     
    else % when testing non-zero e12
        u(7) = 0; u(8) = 0;
        u(9) = e12*sel; u(10) = 0;
        u(4) = 0; u(14) = 0;
        u(3) = []; u(13) = []; % INCOMPLETE  
    end
    
    
%   develop force vector for all nodes (except center node) using sized-down K
%       matrix
    F = K([1:8,11:18],[1:8,11:18])*u; 
    
%   use force/displacement relations to solve at center node
    F_5x = F(1)+F(3)+F(5)+F(7)+F(9)+F(11)+F(13)+F(15); 
    F_5y = F(2)+F(4)+F(6)+F(8)+F(10)+F(12)+F(14)+F(16); 
    u_5 = K(9:10,9:10)\[F_5x;F_5y];
    F = [F(1:8);F_5x;F_5y;F(9:16)];
    u = [u(1:8);u_5(1);u_5(2);u(9:16)];
    
%   use system-wide F matrix and force relations to solve for stress vector
%       [s11,s22,s12]'
    F_x = F(13)+F(15)+F(17); s11 = F_x/(3*A);
    F_y = F(6)+F(12)+F(18); s22 = F_y/(3*A);
    F_xy = F(5)+F(11)+F(17); s12 = F_xy/(3*sel*r);
    stressvec = [s11;s22;s12];
    
%   use strain and stress vectors to solve for the corresponding row of the
%       C matrix
    Cdummy = stressvec/strainvec;
    C(i,:) = Cdummy(:,i)';
end


