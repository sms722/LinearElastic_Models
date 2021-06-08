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

C = [];
for i = 1:3
%   Define strain vector: [e11, e22, e12]'
    strainvec = [0;0;0];
    
%   set that component equal to a dummy value (0.01 strain), 
%       set all other values to zero
    strainvec(i) = 0.01; 
    
%   use strain relations and BCs to develop all displacement values, 
%       except for center node
    e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
    if e11 ~= 0 % when testing non-zero e11
        u_u = [0;0;0;0;0;0;0;e11*sel;0;e11*sel;e11*sel;0];
        K_u = K([1:3,5,6,8,12,13:15,17,18],[1:3,5,6,8,12,13:15,17,18]);
        F_u = K_u*u_u;
        F_F = [0;0;0;0];
        K_F = K([4,7,11,16],[4,7,11,16]);
        u_F = (K_F)\F_F;
        F_5x = F_u(1)+F_u(3)+F_u(4)+F_u(8)+F_u(10)+F_u(11);
        F_5y = F_u(2)+F_u(5)+F_u(6)+F_u(7)+F_u(9)+F_u(12);
        u_5 = K(9:10,9:10)\[F_5x;F_5y];
        F = [F_u(1:3);0;F_u(4:5);0;F_u(6);F_5x;...
             F_5y;0;F_u(7:10);0;F_u(11:12)];
        u = [u_u(1:3);u_F(1);u_u(4:5);u_F(2);u_u(6);u_5(1);...
             u_5(2);u_F(3);u_u(7:10);u_F(4);u_u(11:12)];
    elseif e22 ~= 0 % when testing non-zero e22
        u_u = [0;0;0;0;e22*sel;0;e22*sel;0;0;0;0;e22*sel];
        K_u = K([1:3,5,6,8,12,13:15,17,18],[1:3,5,6,8,12,13:15,17,18]);
        F_u = K_u*u_u;
        F_F = [0;0;0;0];
        K_F = K([4,7,11,16],[4,7,11,16]);
        u_F = (K_F)\F_F;
        F_5x = F_u(1)+F_u(3)+F_u(4)+F_u(8)+F_u(10)+F_u(11);
        F_5y = F_u(2)+F_u(5)+F_u(6)+F_u(7)+F_u(9)+F_u(12);
        u_5 = K(9:10,9:10)\[F_5x;F_5y];
        F = [F_u(1:3);0;F_u(4:5);0;F_u(6);F_5x;...
             F_5y;0;F_u(7:10);0;F_u(11:12)];
        u = [u_u(1:3);u_F(1);u_u(4:5);u_F(2);u_u(6);u_5(1);...
             u_5(2);u_F(3);u_u(7:10);u_F(4);u_u(11:12)];
    else % when testing non-zero e12
        u_u = [0;0;0.5*e12*sel;e12*sel;0;0;e12*sel;0;0;0;0.5*e12*sel;e12*sel;0];
        K_u = K([1:3,5,6,8,11:15,17,18],[1:3,5,6,8,11:15,17,18]);
        F_u = K_u*u_u;
        F_F = [0;0;0];
        K_F = K([4,7,16],[4,7,16]);
        u_F = (K_F)\F_F;
        F_5x = F_u(1)+F_u(3)+F_u(4)+F_u(9)+F_u(11)+F_u(12);
        F_5y = F_u(2)+F_u(5)+F_u(6)+F_u(8)+F_u(10)+F_u(13);
        u_5 = K(9:10,9:10)\[F_5x;F_5y];
        F = [F_u(1:3);0;F_u(4:5);0;F_u(6);F_5x;...
             F_5y;0;F_u(7:10);0;F_u(11:12)];
        u = [u_u(1:3);u_F(1);u_u(4:5);u_F(2);u_u(6);u_5(1);...
             u_5(2);u_u(7:11);u_F(3);u_u(12:13)];
    end

%   use system-wide F matrix and force relations to solve for stress vector
%       [s11,s22,s12]'
    F_x = F(13)+F(15)+F(17); s11 = F_x/(2*sel*2*r);
    F_y = F(6)+F(12)+F(18); s22 = F_y/(2*sel*2*r);
    F_xy = F(5)+F(11)+F(17); s12 = F_xy/(2*sel*2*r);
    stressvec = [s11;s22;s12];
    
%   use strain and stress vectors to solve for the corresponding row of the
%       C matrix
    Cdummy = stressvec/strainvec;
    C(i,:) = Cdummy(:,i)';
end

%{
% Plotting nodal displacement
sf = 5; % scale factor (for magnifying displacement plotting)
newNC = [];
n = 1;
for j = 1:2:size(u)
    newrow = [NC(n,1)+(sf*u(j)),NC(n,2)+(sf*u(j+1))];
    newNC = [newNC;newrow];  
    n = n+1;
end
plot(NC(:,1),NC(:,2),'b.',newNC(:,1),newNC(:,2),'r*');
%}

