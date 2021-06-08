%% 2D 9-node problem

clc;    
close all; 
clear;

sel = 1; % unit cube side length of 1 m/ truss length of 0.5 m
A = 0.0005; % Cross-sectional area of 0.0005 m^2 / 5 cm^2
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

% Assigning values in known-values global u vector
    ux = 0.01; uy = 0.01;
    %u_actual = [0;0;0;[u2y];0;uy;
    %            [u4x];0;[u5x];[u5y];[u6x];uy;
    %            ux;0;ux;[u8y];ux;uy];
    u_u = [0;0;0;0;uy;0;uy;ux;0;ux;ux;uy];

    % Downsizing K to match known-values u vector
    K_u = K([1:3,5,6,8,12,13:15,17,18],[1:3,5,6,8,12,13:15,17,18]);

    % Solving unknown F's with K_u
    F_u = K_u*u_u;

% Assigning values in known-values global F vector
    %F_actual = [[F1x];[F1y];[F2x];0;[F3x];[F3y];
    %            0;[F4y];[F5x];[F5y];0;[F6y];
    %            [F7x];[F7y];[F8x];0;[F9x];[F9y]];
    F_F = [0;0;0;0];

    % Downsizing K to match known-values F vector
    K_F = K([4,7,11,16],[4,7,11,16]);
    
    % Solving unknown u's with K_F
    u_F = (K_F)\F_F;

%Solving for results at the central node (node 5)
    % Calculating forces at node 5
    F_5x = F_u(1)+F_u(3)+F_u(4)+F_u(8)+F_u(10)+F_u(11);
    F_5y = F_u(2)+F_u(5)+F_u(6)+F_u(7)+F_u(9)+F_u(12);

    % Calculating displacements at node 5
    u_5 = K(9:10,9:10)\[F_5x;F_5y];
