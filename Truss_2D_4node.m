%% 2D 4-node problem

clc;    
close all; 
clear;

sel = 1; % unit cube side length of 1 m/ truss length of 0.5 m
A = 0.0005; % Cross-sectional area of 0.0005 m^2 / 5 cm^2
E = 10000; % Young's Modulus of 10 kPa (polymeric material)

NC = sel.*[0,0;0,1;1,1;1,0];
CA = [1,2;2,3;3,4;4,1;1,3];

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

% Assigning values in u vector
ux = 0.01; uy = 0;
u = [0;0;0;uy;ux;uy;ux;0];

% Solving for F values
F = K*u;