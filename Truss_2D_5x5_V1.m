%% 2D 9-node problem

% Initialization
clc;    
close all; 
clear;

% Input Constants
nucFac = 2; % factor to indicate how many unit cells are in this model
            %   (so nucFac = 1 indicates the entire 5x5 is one unit cell,
            %    and nucFac = 2 indicates that the 5x5 contains four
            %    3x3 unit cells inside of it
sel = nucFac*0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                %   area of (assumed circular) truss members
A = pi*(r^2); % Cross-sectional area of truss member
E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

% Input Connectivity Array
%{
% Case 1, 4 unit cells in 2x2
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      1,7;2,7;3,7;6,7;7,8;7,11;7,12;7,13;
      3,9;4,9;5,9;8,9;9,10;9,13;9,14;9,15;
      11,17;12,17;13,17;16,17;17,18;17,21;17,22;17,23;
      13,19;14,19;15,19;18,19;19,20;19,23;19,24;19,25];

% Case 2, 4 unit cells in 2x2
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      2,6;2,7;2,8;6,7;7,8;6,12;7,12;8,12;
      4,8;4,9;4,10;8,9;9,10;8,14;9,14;10,14;
      12,16;12,17;12,18;16,17;17,18;16,22;17,22;18,22;
      14,18;14,19;14,20;18,19;19,20;18,24;19,24;20,24];

% Case 3, 4 unit cells in 2x2 (DNR)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      1,7;3,7;7,11;7,13;3,9;5,9;9,13;9,15;
      11,17;13,17;17,21;17,23;13,19;15,19;19,23;19,25];
%}

% Nodal Coordinate Vector (Standard 5x5, 2D Grid)
NC = sel.*[0,0;0,0.25;0,0.5;0,0.75;0,1;
           0.25,0;0.25,0.25;0.25,0.5;0.25,0.75;0.25,1;
           0.5,0;0.5,0.25;0.5,0.5;0.5,0.75;0.5,1;
           0.75,0;0.75,0.25;0.75,0.5;0.75,0.75;0.75,1;
           1,0;1,0.25;1,0.5;1,0.75;1,1];


% Check stability pre-FEA (function below)
[N,stability_pre] = stabilityTester(CA);
%{
if stability_pre == 0
    disp('This truss design is not mechanically stable');
    return
end
%}


% Develop C-matrix from K-matrix (functions below)
C = [];
[C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C);

% Plot nodal displacement (function below)
plotNDisp(NC,uBasket);

% Print C-matrix as output
disp('The C-matrix is: '); disp(C);

%{
% Post-FEA C-symmetry check (function below)
symm_post = symmCCheck(C);
if symm_post == 0
    disp('This truss design is not eligible due to a lack of homogenizability');
end
%}

%----------%
% FUNCTION TO TEST TRUSS STABILITY
function [N,stability] = stabilityTester(CA)
    stability = [];

    % Add up counters based on connectivities
    [N,~] = histcounts(CA,25);
    %{
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
    %}
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
            qrvec = [4,6,8,11,13,14,15,16,17,18,19,21,23,24,25,26,27,28,...
                     29,31,33,34,35,36,37,38,39,44,46,48,1,2,3,5,7,9,10,...
                     12,20,22,30,32,40,41,42,43,45,47,49,50];
            K = formK(NC,CA,A,E); % function for this below
            newK = [K([4,6,8,11,13:19,21,23:29,31,33:39,44,46,48],:);...
                    K([1:3,5,7,9,10,12,20,22,30,32,40:43,45,47,49,50],:)];
            newK = [newK(:,[4,6,8,11,13:19,21,23:29,31,33:39,44,46,48]),...
                    newK(:,[1:3,5,7,9,10,12,20,22,30,32,...
                    40:43,45,47,49,50])];
            K = newK;
            u_r = [0;0;0;0;0;0;0;0;0;0;0;0;0;e11*sel;0;e11*sel;e11*sel;...
                   e11*sel;e11*sel;0];
            F_q = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;...
                   0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
            K_qq = K(1:30,1:30);
            K_rq = K(31:50,1:30);
            K_qr = K(1:30,31:50);
            K_rr = K(31:50,31:50);
            u_q = K_qq\(F_q-(K_qr*u_r)); 
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end            
        elseif e22 ~= 0 % when testing non-zero e22
            qrvec = [4,6,8,11,13,14,15,16,17,18,19,21,23,24,25,26,27,28,...
                     29,31,33,34,35,36,37,38,39,44,46,48,1,2,3,5,7,9,10,...
                     12,20,22,30,32,40,41,42,43,45,47,49,50];
            K = formK(NC,CA,A,E); % function for this below
            newK = [K([4,6,8,11,13:19,21,23:29,31,33:39,44,46,48],:);...
                    K([1:3,5,7,9,10,12,20,22,30,32,40:43,45,47,49,50],:)];
            newK = [newK(:,[4,6,8,11,13:19,21,23:29,31,33:39,44,46,48]),...
                    newK(:,[1:3,5,7,9,10,12,20,22,30,32,...
                    40:43,45,47,49,50])];
            K = newK;
            u_r = [0;0;0;0;0;0;e22*sel;0;e22*sel;0;e22*sel;0;e22*sel;0;...
                   0;0;0;0;0;e22*sel];
            F_q = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;...
                   0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
            K_qq = K(1:30,1:30); % DOFs 16 23 27 36
            K_rq = K(31:50,1:30);
            K_qr = K(1:30,31:50);
            K_rr = K(31:50,31:50);
            u_q = K_qq\(F_q-(K_qr*u_r));
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end 
        else % when testing non-zero e12
            qrvec = [4,6,8,13,14,15,16,17,18,23,24,25,26,27,28,...
                     33,34,35,36,37,38,44,46,48,1,2,3,5,7,9,10,...
                     11,12,19,20,21,22,29,30,31,32,39,40,41,42,...
                     43,45,47,49,50];
            K = formK(NC,CA,A,E); % function for this below
            newK = [K([4,6,8,13:18,23:28,33:38,44,46,48],:);...
                    K([1:3,5,7,9:12,19:22,29:32,39:43,45,47,49,50],:)];
            newK = [newK(:,[4,6,8,13:18,23:28,33:38,44,46,48]),...
                    newK(:,[1:3,5,7,9:12,19:22,29:32,39:43,45,47,49,50])];
            K = newK;
            u_r = [0;0;0.25*sel*e12;0.5*sel*e12;0.75*sel*e12;sel*e12;0;...
                   0;0;sel*e12;0;0;0;sel*e12;0;0;0;sel*e12;0;0;0;...
                   0.25*sel*e12;0.5*sel*e12;0.75*sel*e12;sel*e12;0];
            F_q = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;...
                   0;0;0;0;0;0;0;0;0;];
            K_qq = K(1:24,1:24);
            K_rq = K(25:50,1:24);
            K_qr = K(1:24,25:50);
            K_rr = K(25:50,25:50);
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
        F_x = F(41)+F(43)+F(45)+F(47)+F(49);
        F_y = F(10)+F(20)+F(30)+F(40)+F(50);
        F_xy = F(9)+F(19)+F(29)+F(39)+F(49);
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
    axis([-0.01 0.12 -0.01 0.12]);
    title('X-Direction Axial Strain (e11)');
    
    subplot(2,2,2);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,2),newNC(:,2,2),'g*');
    axis([-0.01 0.12 -0.01 0.12]);
    title('Y-Direction Axial Strain (e22)');
    
    subplot(2,2,3);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,3),newNC(:,2,3),'k*');
    axis([-0.01 0.12 -0.01 0.12]);
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
