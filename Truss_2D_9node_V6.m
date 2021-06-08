%% 2D 9-node problem

% Initialization
clc;    
close all; 
clear;

% Input Constants
sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                %   area of (assumed circular) truss members
A = pi*(r^2); % Cross-sectional area of truss member
E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

% Input Connectivity Array
%yes-stable
CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
%CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;5,9;3,5;5,7];
%CA = [1,2;2,3;1,4;2,4;2,5;2,6;3,6;4,5;5,6;4,7;4,8;5,8;6,8;6,9;7,8;8,9];

%yes-unstable
%CA = [1,2;2,3;1,4;2,4;3,6;2,6;4,7;4,8;6,9;7,8;8,9;6,8];
%CA = [1,2;2,3;1,4;2,5;3,6;4,5;4,7;5,6;6,9;7,8;8,9;5,8];

%no...
%CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,9;5,7];
%CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;4,8;5,8;6,8;6,9;7,8;8,9];
%CA = [1,2;2,3;1,4;1,5;3,5;3,6;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];

% Nodal Coordinate Vector (Standard 3x3, 2D Grid)
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1];

% Check stability pre-FEA (function below)
[N,stability_pre] = stabilityTester(CA);
if stability_pre == 0
    disp('This truss design is not mechanically stable');
    return
end

% Develop C-matrix from K-matrix (functions below)
C = [];
[C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C);

% Plot nodal displacement (function below)
%plotNDisp(NC,uBasket);

% Print C-matrix as output
disp('The C-matrix is: '); disp(C);

% Post-FEA C-symmetry check (function below)
symm_post = symmCCheck(C);
if symm_post == 0
    disp('This truss design is not eligible due to a lack of homogenizability');
end

%----------%
% FUNCTION TO TEST TRUSS STABILITY
function [N,stability] = stabilityTester(CA)
    stability = [];

    % Add up counters based on connectivities
    [N,~] = histcounts(CA,9);
    
    % Determine stability based on counter values
    if N(5)>=3
        if (N(1)>=3)||(N(7)>=3)||(N(3)>=3)||(N(9)>=3)
            if (N(1)>=3)&&(N(3)>=3)
                stability = true;
            elseif (N(7)>=3)&&(N(3)>=3)
                stability = true;
            elseif (N(9)>=3)&&(N(3)>=3)
                stability = true;
            elseif (N(1)>=3)&&(N(7)>=3)
                stability = true;
            elseif (N(1)>=3)&&(N(9)>=3)
                stability = true;
            elseif (N(7)>=3)&&(N(9)>=3)
                stability = true;
            elseif (N(2)>=3)||(N(4)>=3)||(N(6)>=3)||(N(8)>=3)
                if (N(2)>=3)&&(N(4)>=3)
                    stability = true;
                elseif (N(6)>=3)&&(N(4)>=3)
                    stability = true;
                elseif (N(8)>=3)&&(N(4)>=3)
                    stability = true;
                elseif (N(8)>=3)&&(N(2)>=3)
                    stability = true;
                elseif (N(8)>=3)&&(N(6)>=3)
                    stability = true;
                elseif (N(6)>=3)&&(N(2)>=3)
                    stability = true;
                else
                    stability = false;
                end
            else
                stability = false;
            end    
        elseif (N(2)>=4)||(N(4)>=4)||(N(6)>=4)||(N(8)>=4)
            if (N(2)>=4)&&(N(4)>=4)&&(N(6)>=4)
                stability = true;
            elseif (N(4)>=4)&&(N(6)>=4)&&(N(8)>=4)
                stability = true;
            elseif (N(2)>=4)&&(N(6)>=4)&&(N(8)>=4)
                stability = true;
            elseif (N(2)>=4)&&(N(4)>=4)&&(N(8)>=4)
                stability = true;
            else
                stability = false;
            end
        else
            stability = false;
        end
    else
        stability = false;
    end
end

% FUNCTION TO CALCULATE C-MATRIX
function [C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C)
    % Initialize outputs
    uBasket = []; FBasket = [];
    
    % Iterate through once for each strain component
    for y = 1:1:3
    %   Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

    %   set that component equal to a dummy value (0.01 strain), 
    %       set all other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;

    %   use strain relations, BCs, and partitioned K-matrix to 
    %       solve for all unknowns
        e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
        if e11 ~= 0 % when testing non-zero e11
            qrvec = [4,7,9,10,11,16,1,2,3,5,6,8,12,13,14,15,17,18];
            K = formK(NC,CA,A,E); % function for this below
            newK = [K([4,7,9,10,11,16],:);...
                    K([1,2,3,5,6,8,12,13,14,15,17,18],:)];
            newK = [newK(:,[4,7,9,10,11,16]),...
                    newK(:,[1,2,3,5,6,8,12,13,14,15,17,18])];
            K = newK;
            u_r = [0;0;0;0;0;0;0;e11*sel;0;e11*sel;e11*sel;0];
            F_q = [0;0;0;0;0;0];
            K_qq = K(1:6,1:6);
            K_rq = K(7:18,1:6);
            K_qr = K(1:6,7:18);
            K_rr = K(7:18,7:18);
            u_q = K_qq\(F_q-(K_qr*u_r)); 
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end            
        elseif e22 ~= 0 % when testing non-zero e22
            K = formK(NC,CA,A,E); % function for this below
            newK = [K([4,7,9,10,11,16],:);...
                    K([1,2,3,5,6,8,12,13,14,15,17,18],:)];
            newK = [newK(:,[4,7,9,10,11,16]),...
                    newK(:,[1,2,3,5,6,8,12,13,14,15,17,18])];
            K = newK;
            u_r = [0;0;0;0;e22*sel;0;e22*sel;0;0;0;0;e22*sel];
            F_q = [0;0;0;0;0;0];
            K_qq = K(1:6,1:6);
            K_rq = K(7:18,1:6);
            K_qr = K(1:6,7:18);
            K_rr = K(7:18,7:18);
            u_q = K_qq\(F_q-(K_qr*u_r));
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end 
        else % when testing non-zero e12
            qrvec = [4,9,10,16,1,2,3,5,6,7,8,11,12,13,14,15,17,18];
            K = formK(NC,CA,A,E);            
            newK = [K([4,9,10,16],:);K([1:3,5:8,11:15,17,18],:)];
            newK = [newK(:,[4,9,10,16]),newK(:,[1:3,5:8,11:15,17,18])];
            K = newK;
            u_r = [0;0;0.5*e12*sel;e12*sel;0;0;0;e12*sel;0;0;...
                   0;0.5*e12*sel;e12*sel;0];
            F_q = [0;0;0;0];
            K_qq = K(1:4,1:4);
            K_rq = K(5:18,1:4);
            K_qr = K(1:4,5:18);
            K_rr = K(5:18,5:18);
            u_q = K_qq\(F_q-(K_qr*u_r));
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end 
        end

    %   use system-wide F matrix and force relations to solve for 
    %       stress vector [s11,s22,s12]'
        F_x = F(13)+F(15)+F(17);
        F_y = F(6)+F(12)+F(18);
        F_xy = F(5)+F(11)+F(17);
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

    %   use strain and stress vectors to solve for the corresponding
    %       row of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
        FBasket(:,y) = F;
        uBasket(:,y) = u;
    end
end

% FUNCTION TO FORM GLOBAL STRUCTURAL STIFFNESS MATRIX
function K = formK(NC,CA,A,E)
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
end

% FUNCTION TO PLOT NODAL DISPLACEMENT
function plotNDisp(NC,uBasket)
    sf = 5; % scale factor (for magnifying displacement plotting)
    newNC = [];
    n = 1;
    for i = 1:1:size(uBasket,2)
        u = uBasket(:,i);
        nnc = [];
        n = 1;
        for j = 1:2:size(u)
            newrow = [NC(n,1)+(sf*u(j)),NC(n,2)+(sf*u(j+1))];
            nnc = [nnc;newrow];
            n = n+1;
        end
        newNC(:,:,i) = nnc;
    end
    subplot(2,2,1);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,1),newNC(:,2,1),'r*');
    axis([-0.01 0.06 -0.01 0.06]);
    title('X-Direction Axial Strain (e11)');
    
    subplot(2,2,2);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,2),newNC(:,2,2),'g*');
    axis([-0.01 0.06 -0.01 0.06]);
    title('Y-Direction Axial Strain (e22)');
    
    subplot(2,2,3);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,3),newNC(:,2,3),'k*');
    axis([-0.01 0.06 -0.01 0.06]);
    title('Shear Strain (e12)');
end

% FUNCTION TO TEST C-MATRIX SYMMETRY
function symm_post = symmCCheck(C)
    if C(1,2) == C(2,1)
        if C(3,1) == C(1,3)
            if C(2,3) == C(3,2)
                symm_post = true;
            else
                symm_post = false;
            end
        else
            symm_post = false;
        end
    else
        symm_post = false;
    end
end
