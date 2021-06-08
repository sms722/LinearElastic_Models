% 2D NxN Square Truss Problem

% Initialization
clc;    
close all; 
clear;

% Input Constants
sidenum = 5; % Number of nodes on each size of square grid
nucFac = 1; % factor to indicate how many unit cells are in this model
            %   (so nucFac = 1 indicates the entire 5x5 is one unit cell,
            %    and nucFac = 2 indicates that the 5x5 contains four
            %    3x3 unit cells inside of it
sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                %   area of (assumed circular) truss members
A = pi*(r^2); % Cross-sectional area of truss member
E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

% Input Connectivity Array
%{
% Case 1, 4 unit cells in 2x2 (5x5 grid)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      1,7;2,7;3,7;6,7;7,8;7,11;7,12;7,13;
      3,9;4,9;5,9;8,9;9,10;9,13;9,14;9,15;
      11,17;12,17;13,17;16,17;17,18;17,21;17,22;17,23;
      13,19;14,19;15,19;18,19;19,20;19,23;19,24;19,25];


% Case 2, 4 unit cells in 2x2 (5x5 grid)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      2,6;2,7;2,8;6,7;7,8;6,12;7,12;8,12;
      4,8;4,9;4,10;8,9;9,10;8,14;9,14;10,14;
      12,16;12,17;12,18;16,17;17,18;16,22;17,22;18,22;
      14,18;14,19;14,20;18,19;19,20;18,24;19,24;20,24];

% Case 3, 4 unit cells in 2x2 (5x5 grid)(unstable, doesn't run)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      1,7;3,7;7,11;7,13;3,9;5,9;9,13;9,15;
      11,17;13,17;17,21;17,23;13,19;15,19;19,23;19,25];

% Case 1, 9 unit cells in 3x3 (7x7 grid)
CA = [1,2;2,3;3,4;4,5;5,6;6,7;7,14;14,21;21,28;28,35;35,42;42,49;
      49,48;48,47;47,46;46,45;45,44;44,43;36,43;29,36;22,29;15,22;8,15;1,8;
      1,9;2,9;3,9;8,9;9,10;9,15;9,16;9,17;
      3,11;4,11;5,11;10,11;11,12;11,17;11,18;11,19;
      5,13;6,13;7,13;12,13;13,14;13,19;13,20;13,21;
      15,23;16,23;17,23;22,23;23,24;23,29;23,30;23,31;
      17,25;18,25;19,25;24,25;25,26;25,31;25,32;25,33;
      19,27;20,27;21,27;26,27;27,28;27,33;27,34;27,35;
      29,37;30,37;31,37;36,37;37,38;37,43;37,44;37,45;
      31,39;32,39;33,39;38,39;39,40;39,45;39,46;39,47;
      33,41;34,41;35,41;40,41;41,42;41,47;41,48;41,49;
      3,10;10,17;17,24;24,31;31,38;38,45;
      5,12;12,19;19,26;26,33;33,40;40,47;
      15,16;16,17;17,18;18,19;19,20;20,21;
      29,30;30,31;31,32;32,33;33,34;34,35];

% Case 6, 1 unit cell in 4x4 (4x4 grid)(nonsymmetric)
CA = [1,2;2,3;3,4;4,8;8,12;12,16;15,16;14,15;13,14;9,13;5,9;1,5;
      2,6;6,10;10,14;3,7;7,11;11,15;5,6;6,7;7,8;9,10;10,11;11,12;
      1,6;5,10;9,14;2,7;6,11;10,15;3,8;7,12;11,16];

% Case 7, 1 unit cell in 5x5 (5x5 grid)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      1,7;2,7;3,7;6,7;7,8;7,11;7,12;7,13;
      3,9;4,9;5,9;8,9;9,10;9,13;9,14;9,15;
      12,16;12,17;12,18;16,17;17,18;16,22;17,22;18,22;
      14,18;14,19;14,20;18,19;19,20;18,24;19,24;20,24];

% Case 8, 1 unit cell in 5x5 (5x5 grid)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      1,7;2,7;3,7;6,7;7,8;7,11;7,12;7,13;
      4,8;4,9;4,10;8,9;9,10;8,14;9,14;10,14;
      11,17;12,17;13,17;16,17;17,18;17,21;17,22;17,23;
      14,18;14,19;14,20;18,19;19,20;18,24;19,24;20,24];
%}

% Case 9, 1 unit cell in 5x5 (5x5 grid)
CA = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
      25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
      3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
      2,6;2,7;2,8;6,7;7,8;6,12;7,12;8,12;
      4,8;4,9;4,10;8,9;9,10;8,14;9,14;10,14;
      12,16;12,17;16,17;17,18;16,22;17,22;18,22;
      14,19;14,20;18,19;19,20;18,24;19,24;20,24];
%}
% Generate vector with nodal coordinates
NC = generateNC(sel,sidenum);

% Check stability pre-FEA (function below)
[N,stability_pre] = stabilityTester(sidenum,CA,NC);
if stability_pre == 0
    disp('This truss design is not mechanically stable');
    return
end

% Develop C-matrix from K-matrix (functions below)
C = [];
[C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C);
C = C./nucFac;

%-% Calculate volume fraction (function below)
volFrac = calcVF(NC,CA,r,sel);

% Print C-matrix as output
disp('The C-matrix is: '); disp(C);

% Plot nodal displacement (function below)
%plotNDisp(sel,NC,uBasket);

% Post-FEA C-symmetry check (function below)
[C,symm_post] = symmCCheck(C);
if symm_post == 0
    disp('Truss design not eligible due to a lack of homogenizability');
end


%----------%
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

% FUNCTION TO TEST TRUSS STABILITY (INCOMPLETE/INVALID!!!)
function [N,stability] = stabilityTester(sidenum,CA,NC)
    stability = true;

    % Add up counters based on connectivities
    [N,~] = histcounts(CA,size(NC,1));
    
    % First stability check: number of "holes" (unconnected nodes) in truss
    %   should be less than or equal to [(number of side nodes) - 2]
    zeros = find(~N);
    if length(zeros) >= (sidenum-2)
        stability = false;
        return
    end
    
    % Second stability check: nodes with connections are connected to at
    %   least three other nodes apiece (except for the corner nodes)
    Ns = N([2:(sidenum-1),(sidenum+1):((sidenum^2)-sidenum),...
         ((sidenum^2)-(sidenum-2)):(sidenum^2)-1]);
    Nnz = Ns(Ns>0);
    for a = 1:1:length(Nnz)
       if Nnz(a) < 3
           stability = false;
           break
       end
    end
end

% FUNCTION TO CALCULATE C-MATRIX
function [C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C)
    % Initialize outputs
    uBasket = []; FBasket = [];
    
    % Iterate through once for each strain component
    for y = 1:1:3
    %   Define vectors to hold indexes for output forces
        Fi_x = []; Fi_y = []; Fi_xy = [];
    
    %   Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

    %   set that component equal to a dummy value (0.01 strain), 
    %       set all other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;

    %   use strain relations, BCs, and partitioned K-matrix to 
    %       solve for all unknowns
        e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
        if (e11 ~= 0) || (e22 ~= 0) % when testing non-zero e11 or e22
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                ND = NC./sel;
                % Separating conditions for exterior nodes
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((2*x)-1)];
                    elseif ND(x,1) == 1
                        % displacement in x = e11*sel
                        u_r = [u_r;(e11*sel)];
                        rvec = [rvec,((2*x)-1)];
                        Fi_x = [Fi_x,((2*x)-1)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((2*x)-1)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(2*x)];
                    elseif ND(x,2) == 1
                        % displacement in y = e22*sel
                        u_r = [u_r;(e22*sel)];
                        rvec = [rvec,(2*x)];
                        Fi_y = [Fi_y,(2*x)];
                        Fi_xy = [Fi_xy,((2*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(2*x)];
                    end
                else % Condition for all interior nodes
                    % both x- and y-DOFs are force BCs
                    F_q = [F_q;0;0];
                    qvec = [qvec,((2*x)-1),(2*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(2*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(2*size(NC,1)));
            K_rr = newK((length(qvec)+1):(2*size(NC,1)),...
                   (length(qvec)+1):(2*size(NC,1)));
            u_q = K_qq\(F_q-(K_qr*u_r)); 
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        else % when testing non-zero e12
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes 
                ND = NC./sel;
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e12*sel*ND(x,2))];
                        rvec = [rvec,((2*x)-1)];
                    elseif ND(x,1) == 1
                        % displacement in x is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e12*sel*ND(x,2))];
                        rvec = [rvec,((2*x)-1)];
                        Fi_x = [Fi_x,((2*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in x = e12*sel
                        u_r = [u_r;(e12*sel)];
                        rvec = [rvec,((2*x)-1)];
                    elseif ND(x,2) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((2*x)-1)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((2*x)-1)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(2*x)];
                    elseif ND(x,2) == 1
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(2*x)];
                        Fi_y = [Fi_y,(2*x)];
                        Fi_xy = [Fi_xy,((2*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(2*x)];
                    end
                else % Blanket condition for all interior nodes
                    % both x- and y-DOFs are force BCs
                    F_q = [F_q;0;0];
                    qvec = [qvec,((2*x)-1),(2*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(2*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(2*size(NC,1)));
            K_rr = newK((length(qvec)+1):(2*size(NC,1)),...
                   (length(qvec)+1):(2*size(NC,1)));
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
        F_x = 0; F_y = 0; F_xy = 0;
        for n = 1:1:size(Fi_xy,2)
            F_x = F_x + F(Fi_x(n));
            F_y = F_y + F(Fi_y(n));
            F_xy = F_xy + F(Fi_xy(n));
        end
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

% FUNCTION TO CALCULATE VOLUME FRACTION
function volFrac = calcVF(NC,CA,r,sel)
    totalTrussVol = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*pi*(r^2));
    end
    % Calculating volume fraction (using a solid square with 2*r thickness
    %   as a baseline)
    volFrac = totalTrussVol/(2*r*(sel^2));
end

% FUNCTION TO PLOT NODAL DISPLACEMENT
function plotNDisp(sel,NC,uBasket)
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
    
    min = NC(1,1)-(sel*0.1);
    max = NC(size(NC,1),1)+(sel*0.2);
    
    subplot(2,2,1);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,1),newNC(:,2,1),'r*');
    axis([min max min max]);
    title('X-Direction Axial Strain (e11)');
    
    subplot(2,2,2);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,2),newNC(:,2,2),'g*');
    axis([min max min max]);
    title('Y-Direction Axial Strain (e22)');
    
    subplot(2,2,3);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,3),newNC(:,2,3),'k*');
    axis([min max min max]);
    title('Shear Strain (e12)');
end

% FUNCTION TO TEST C-MATRIX SYMMETRY
function [C,symm_post] = symmCCheck(C)
    for x = 1:1:(size(C,1)*size(C,2))
        if C(x) < 10^(-10)
            C(x) = 0;
        end
    end
    if abs(C(1,2)-C(2,1)) < 10^(-10)
        if abs(C(1,3)-C(3,1)) < 10^(-10)
            if abs(C(2,3)-C(3,2)) < 10^(-10)
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
