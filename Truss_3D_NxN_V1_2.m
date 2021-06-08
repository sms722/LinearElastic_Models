% 3D NxN Square Truss Problem

% Initialization
clc;    
close all; 
clear;

% Input Constants
sidenum = 3; % Number of nodes on each edge of cubic grid
nucFac = 1; % factor to indicate how many unit cells are in this model
            %   (so nucFac = 1 indicates the entire 5x5x5 is one unit cell,
            %    and nucFac = 2 indicates that the 5x5x5 contains eight
            %    3x3x3 unit cells inside of it
sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                %   area of (assumed circular) truss members
A = pi*(r^2); % Cross-sectional area of truss member
E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

% Input Connectivity Array

% Case 1-3D: 1 unit cell in 3x3x3, with Case 1-2D unit cells in all planes

CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;1,5;5,9;3,5;5,7;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;11,14;14,17;13,14;
      14,15;10,14;14,18;12,14;14,16;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;20,23;23,26;22,23;
      23,24;19,23;23,27;21,23;23,25;
      1,10;10,19;2,11;11,20;3,12;12,21;4,13;13,22;5,14;14,23;6,15;15,24;
      7,16;16,25;8,17;17,26;9,18;18,27;
      7,17;17,27;9,17;17,25;4,14;14,24;6,14;14,22;1,11;11,21;3,11;11,19;
      7,13;13,19;1,13;13,25;8,14;14,20;2,14;14,26;9,15;15,21;3,15;15,27;
      1,14;3,14;7,14;9,14;14,19;14,21;14,25;14,27];
%{
% Case 2-3D: 1 unit cell in 3x3x3, with Case 2-2D unit cells in all planes
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;2,4;2,6;6,8;4,8;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;11,14;14,17;13,14;
      14,15;11,13;11,15;15,17;13,17;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;20,23;23,26;22,23;
      23,24;20,22;20,24;24,26;22,26;
      1,10;10,19;2,11;11,20;3,12;12,21;4,13;13,22;5,14;14,23;6,15;15,24;
      7,16;16,25;8,17;17,26;9,18;18,27;
      8,16;8,18;16,26;18,26;5,13;5,15;15,23;13,23;2,10;2,12;10,20;12,20;
      4,10;10,22;16,22;4,16;5,11;11,23;17,23;5,17;6,12;12,24;18,24;6,18;
      2,13;13,26;15,26;2,15;8,13;13,20;15,20;8,15];

% Case 3-3D: 1 unit cell in 3x3x3, with Case 1-2D unit cells along
% diagonals
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;
      11,14;14,17;13,14;14,15;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;
      20,23;23,26;22,23;23,24;
      1,10;10,19;2,11;11,20;3,12;12,21;6,15;15,24;
      9,18;18,27;8,17;17,26;7,16;16,25;4,13;13,22;5,14;14,23;
      1,14;3,14;7,14;9,14;14,19;14,21;14,25;14,27];

% Case 4-3D: 1 unit cell in 3x3x3, with only 4 diagonals
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;
      11,14;14,17;13,14;14,15;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;
      20,23;23,26;22,23;23,24;
      1,10;10,19;2,11;11,20;3,12;12,21;6,15;15,24;
      9,18;18,27;8,17;17,26;7,16;16,25;4,13;13,22;5,14;14,23;
      1,14;3,14;14,25;14,27];

% Case 5-3D: 1 unit cell in 3x3x3, consisting of 8 2x2x2 cubes and an
% octahedron in the middle
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;11,14;14,17;13,14;
      14,15;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;20,23;23,26;22,23;
      23,24;
      1,10;10,19;2,11;11,20;3,12;12,21;4,13;13,22;5,14;14,23;6,15;15,24;
      7,16;16,25;8,17;17,26;9,18;18,27;
      2,13;13,26;15,26;2,15;8,13;13,20;15,20;8,15];
  
% Case 6-3D: 1 unit cell in 3x3x3, essentially Case 1-3D minus some
% 3D diagonals
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;1,5;5,9;3,5;5,7;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;11,14;14,17;13,14;
      14,15;10,14;14,18;12,14;14,16;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;20,23;23,26;22,23;
      23,24;19,23;23,27;21,23;23,25;
      1,10;10,19;2,11;11,20;3,12;12,21;4,13;13,22;5,14;14,23;6,15;15,24;
      7,16;16,25;8,17;17,26;9,18;18,27;
      7,17;17,27;9,17;17,25;4,14;14,24;6,14;14,22;1,11;11,21;3,11;11,19;
      7,13;13,19;1,13;13,25;8,14;14,20;2,14;14,26;9,15;15,21;3,15;15,27;
      1,14;3,14;7,14;14,21;14,25;14,27];
%}

% Generate vector with nodal coordinates
NC = generateNC(sel,sidenum);

% Check stability pre-FEA (function below)
%{
[N,stability_pre] = stabilityTester(sidenum,CA,NC);
if stability_pre == 0
    disp('This truss design is not mechanically stable');
    return
end
%}

% Develop C-matrix from K-matrix (functions below)
C = [];
[C,uBasket,FBasket] = generateC(sel,NC,CA,A,E,C);
C = C./nucFac;

% Print C-matrix as output
disp('The C-matrix is: '); disp(C);

% Plot nodal displacement (function below)
%plotNDisp(sel,NC,uBasket);

% Post-FEA C-symmetry check (function below)
[C,symm_post] = symmCCheck(C);
if symm_post == 0
    disp('Truss design not eligible due to a lack of homogenizability');
end

% Calculate and print volume fraction (function below)
volFrac = calcVF(NC,CA,r,sel);
%disp('The volume fraction is: '); disp(volFrac);


%----------%
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            for k = 1:1:sidenum
                NC = [NC;notchvec(i),notchvec(j),notchvec(k)];
            end
        end
    end
    NC = sel.*NC;
end

% FUNCTION TO TEST TRUSS STABILITY (NOTCONVYET)
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
    
    % Third stability check: at least one diagonal member present
    nodiags = false;
    for i = 1:1:size(CA,1)
        if CA(i,1)+sidenum+1 == CA(i,2)
            nodiags = true;
        elseif CA(i,1)+sidenum-1 == CA(i,2)
            nodiags = true;
        end
    end
    if nodiags == false
        stability = false;
        return
    end
end

% FUNCTION TO CALCULATE C-MATRIX 
function [C,uBasket,FBasket] = generateC(sel,NC,CA,A,E,C)
    % Initialize variables
    uBasket = []; FBasket = []; ND = NC./sel;
    
    % Iterate through once for each strain component
    for y = 1:1:6
    %   Define vectors to hold indexes for output forces
        Fi_x = []; Fi_y = []; Fi_z = []; 
        Fi_xy = []; Fi_yz = []; Fi_xz = [];
    
    %   Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0;0;0;0];

    %   set that component equal to a dummy value (0.01 strain), 
    %       set all other values to zero
        strainvec(y) = 0.001; 
        strainvec(4) = strainvec(4)*2;
        strainvec(5) = strainvec(5)*2;
        strainvec(6) = strainvec(6)*2;

    %   use strain relations, BCs, and partitioned K-matrix to 
    %       solve for all unknowns
        e11 = strainvec(1); e22 = strainvec(2); e33 = strainvec(3);
        e12 = strainvec(4); e23 = strainvec(5); e31 = strainvec(6);
        % when testing non-zero e11 or e22 or e33:
        if e11 ~= 0
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true) || ...
                (ismember(ND(x,3),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,1) == 1
                        % displacement in x = e11*sel
                        u_r = [u_r;(e11*sel)];
                        rvec = [rvec,((3*x)-2)];
                        Fi_x = [Fi_x,((3*x)-2)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif (ND(x,3) == 0) && (e33 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-2)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in y = e22*sel
                        u_r = [u_r;(e22*sel)];
                        rvec = [rvec,((3*x)-1)];
                        Fi_y = [Fi_y,((3*x)-1)];
                        Fi_xy = [Fi_xy,((3*x)-2)];
                        Fi_yz = [Fi_yz,(3*x)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif (ND(x,3) == 0) && (e33 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-1)];
                    end
                    
                    % Finding z-DOF
                    if ND(x,3) == 0
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif ND(x,3) == 1
                        % displacement in z = e33*sel
                        u_r = [u_r;(e33*sel)];
                        rvec = [rvec,(3*x)];
                        Fi_z = [Fi_z,(3*x)];
                        Fi_xz = [Fi_xz,((3*x)-2)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    else
                        % z-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(3*x)];
                    end
                    
                else % Condition for all interior nodes
                    % x-, y-, and z-DOFs are all force BCs
                    F_q = [F_q;0;0;0];
                    qvec = [qvec,((3*x)-2),((3*x)-1),(3*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(3*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(3*size(NC,1)));
            K_rr = newK((length(qvec)+1):(3*size(NC,1)),...
                   (length(qvec)+1):(3*size(NC,1)));
            u_q = full(sparseinv(K_qq)*(F_q-(K_qr*u_r)));
            F_r = full((K_rq*u_q)+(K_rr*u_r));
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        elseif e22 ~= 0
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true) || ...
                (ismember(ND(x,3),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,1) == 1
                        % displacement in x = e11*sel
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                        Fi_x = [Fi_x,((3*x)-2)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif (ND(x,3) == 0) && (e33 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-2)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in y = e22*sel
                        u_r = [u_r;(e22*sel)];
                        rvec = [rvec,((3*x)-1)];
                        Fi_y = [Fi_y,((3*x)-1)];
                        Fi_xy = [Fi_xy,((3*x)-2)];
                        Fi_yz = [Fi_yz,(3*x)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif (ND(x,3) == 0) && (e33 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-1)];
                    end
                    
                    % Finding z-DOF
                    if ND(x,3) == 0
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif ND(x,3) == 1
                        % displacement in z = e33*sel
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                        Fi_z = [Fi_z,((3*x))];
                        Fi_xz = [Fi_xz,((3*x)-2)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    else
                        % z-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(3*x)];
                    end
                    
                else % Condition for all interior nodes
                    % x-, y-, and z-DOFs are all force BCs
                    F_q = [F_q;0;0;0];
                    qvec = [qvec,((3*x)-2),((3*x)-1),(3*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(3*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(3*size(NC,1)));
            K_rr = newK((length(qvec)+1):(3*size(NC,1)),...
                   (length(qvec)+1):(3*size(NC,1)));
            u_q = full(sparseinv(K_qq)*(F_q-(K_qr*u_r)));
            F_r = full((K_rq*u_q)+(K_rr*u_r));
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        elseif e33 ~= 0
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true) || ...
                (ismember(ND(x,3),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,1) == 1
                        % displacement in x = e11*sel
                        u_r = [u_r;(e11*sel)];
                        rvec = [rvec,((3*x)-2)];
                        Fi_x = [Fi_x,((3*x)-2)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif (ND(x,3) == 0) && (e33 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-2)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in y = e22*sel
                        u_r = [u_r;(e22*sel)];
                        rvec = [rvec,((3*x)-1)];
                        Fi_y = [Fi_y,((3*x)-1)];
                        Fi_xy = [Fi_xy,((3*x)-2)];
                        Fi_yz = [Fi_yz,(3*x)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif (ND(x,3) == 0) && (e33 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-1)];
                    end
                    
                    % Finding z-DOF
                    if ND(x,3) == 0
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif ND(x,3) == 1
                        % displacement in z = e33*sel
                        u_r = [u_r;(e33*sel)];
                        rvec = [rvec,(3*x)];
                        Fi_z = [Fi_z,(3*x)];
                        Fi_xz = [Fi_xz,((3*x)-2)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    else
                        % z-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(3*x)];
                    end
                    
                else % Condition for all interior nodes
                    % x-, y-, and z-DOFs are all force BCs
                    F_q = [F_q;0;0;0];
                    qvec = [qvec,((3*x)-2),((3*x)-1),(3*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(3*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(3*size(NC,1)));
            K_rr = newK((length(qvec)+1):(3*size(NC,1)),...
                   (length(qvec)+1):(3*size(NC,1)));
            u_q = full(sparseinv(K_qq)*(F_q-(K_qr*u_r)));
            F_r = full((K_rq*u_q)+(K_rr*u_r));
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        elseif (e12 ~= 0) % when testing non-zero e12
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes 
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true) || ...
                (ismember(ND(x,3),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e12*sel*ND(x,2))];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,1) == 1
                        % displacement in x is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e12*sel*ND(x,2))];
                        rvec = [rvec,((3*x)-2)];
                        Fi_x = [Fi_x,((3*x)-2)];
                    elseif ND(x,2) == 1
                        % displacement in x = e12*sel
                        u_r = [u_r;(e12*sel)];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,2) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-2)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                        Fi_y = [Fi_y,((3*x)-1)];
                        Fi_xy = [Fi_xy,((3*x)-2)];
                        Fi_yz = [Fi_yz,(3*x)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-1)];
                    end
                    
                    % Finding z-DOF
                    if ND(x,3) == 0
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x))];
                    elseif ND(x,3) == 1
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x))];
                        Fi_z = [Fi_z,(3*x)];
                        Fi_xz = [Fi_xz,((3*x)-2)];
                    else
                        % z-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x))];
                    end
                else % Condition for all interior nodes
                    % x-, y-, and z-DOFs are all force BCs
                    F_q = [F_q;0;0;0];
                    qvec = [qvec,((3*x)-2),((3*x)-1),(3*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(3*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(3*size(NC,1)));
            K_rr = newK((length(qvec)+1):(3*size(NC,1)),...
                   (length(qvec)+1):(3*size(NC,1)));
            u_q = full(sparseinv(K_qq)*(F_q-(K_qr*u_r)));
            F_r = full((K_rq*u_q)+(K_rr*u_r));
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        elseif (e23 ~= 0) % when testing non-zero e23 
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes 
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true) || ...
                (ismember(ND(x,3),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,1) == 1
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                        Fi_x = [Fi_x,((3*x)-2)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-2)];
                    end
            
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                        Fi_y = [Fi_y,((3*x)-1)];
                        Fi_xy = [Fi_xy,((3*x)-2)];
                        Fi_yz = [Fi_yz,(3*x)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-1)];
                    end
                    
                    % Finding z-DOF
                    if ND(x,3) == 0
                        % displacement in z is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e23*sel*ND(x,2))];
                        rvec = [rvec,(3*x)];
                    elseif ND(x,3) == 1
                        % displacement in z is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e23*sel*ND(x,2))];
                        rvec = [rvec,(3*x)];
                        Fi_z = [Fi_z,(3*x)];
                        Fi_xz = [Fi_xz,((3*x)-2)];
                    elseif ND(x,2) == 1
                        % displacement in z = e23*sel
                        u_r = [u_r;(e23*sel)];
                        rvec = [rvec,(3*x)];
                    elseif ND(x,2) == 0
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    else
                        % z-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(3*x)];
                    end
                    
                else % Condition for all interior nodes
                    % x-, y-, and z-DOFs are all force BCs
                    F_q = [F_q;0;0;0];
                    qvec = [qvec,((3*x)-2),((3*x)-1),(3*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(3*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(3*size(NC,1)));
            K_rr = newK((length(qvec)+1):(3*size(NC,1)),...
                   (length(qvec)+1):(3*size(NC,1)));
            u_q = full(sparseinv(K_qq)*(F_q-(K_qr*u_r)));
            F_r = full((K_rq*u_q)+(K_rr*u_r));
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        else % when testing non-zero e31
            K = formK(NC,CA,A,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes 
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true) || ...
                (ismember(ND(x,3),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x is proportional to z-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e31*sel*ND(x,3))];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,1) == 1
                        % displacement in x is proportional to z-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e31*sel*ND(x,3))];
                        rvec = [rvec,((3*x)-2)];
                        Fi_x = [Fi_x,((3*x)-2)];
                    elseif ND(x,3) == 1
                        % displacement in x = e31*sel
                        u_r = [u_r;(e31*sel)];
                        rvec = [rvec,((3*x)-2)];
                    elseif ND(x,3) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-2)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-2)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((3*x)-1)];
                        Fi_y = [Fi_y,((3*x)-1)];
                        Fi_xy = [Fi_xy,((3*x)-2)];
                        Fi_yz = [Fi_yz,(3*x)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((3*x)-1)];
                    end
                    
                    % Finding z-DOF
                    if ND(x,3) == 0
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                    elseif ND(x,3) == 1
                        % displacement in z = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(3*x)];
                        Fi_z = [Fi_z,(3*x)];
                        Fi_xz = [Fi_xz,((3*x)-2)];
                    else
                        % z-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(3*x)];
                    end
                else % Condition for all interior nodes
                    % x-, y-, and z-DOFs are all force BCs
                    F_q = [F_q;0;0;0];
                    qvec = [qvec,((3*x)-2),((3*x)-1),(3*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(3*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(3*size(NC,1)));
            K_rr = newK((length(qvec)+1):(3*size(NC,1)),...
                   (length(qvec)+1):(3*size(NC,1)));
            u_q = full(sparseinv(K_qq)*(F_q-(K_qr*u_r)));
            F_r = full((K_rq*u_q)+(K_rr*u_r));
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        end

    %   use system-wide F matrix and force relations to solve for 
    %       stress vector [s11,s22,s12]'
        F_x = 0; F_y = 0; F_z = 0; 
        F_xy = 0; F_yz = 0; F_xz = 0;
        for n = 1:1:size(Fi_z,2)
            F_x = F_x + F(Fi_x(n));
            F_y = F_y + F(Fi_y(n));
            F_z = F_z + F(Fi_z(n));
            F_xy = F_xy + F(Fi_xy(n));
            F_yz = F_yz + F(Fi_yz(n));
            F_xz = F_xz + F(Fi_xz(n));
        end
        %{
        %ruh-roh
        if (e22 ~= 0)
            F_z = F_x;
        end
        %}
        stressvec = (1/(sel^2)).*[(F_x);(F_y);(F_z);...
                    (F_xy);(F_yz);(F_xz)];

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
        z1 = NC(CA(i,1),3); z2 = NC(CA(i,2),3);
        
        % Plotting elements (comment out if undesired)
        line([x1,x2],[y1,y2],[z1,z2]);
        
        L = sqrt(((x2-x1)^2)+((y2-y1)^2)+((z2-z1)^2));
        c=(x2-x1)/L; c2 = c^2;
        s=(y2-y1)/L; s2 = s^2;
        n=(z2-z1)/L; n2 = n^2;
        ktemp = [c2,   c*s,  c*n,  -c2,  -c*s, -c*n;  
                 c*s,  s2,   s*n,  -c*s, -s2,  -s*n;  
                 c*n,  s*n,  n2,   -c*n, -s*n, -n2 ;   
                 -c2,  -c*s, -c*n, c2,   c*s,  c*n ; 
                 -c*s, -s2,  -s*n, c*s,  s2,   s*n ;
                 -c*n, -s*n, -n2,  c*n,  c*n,  n2  ];
        ke = ((A.*E)./L).*ktemp;
        Kbasket(:,:,i) = ke;
    end

    % Global-to-local-coordinate-system Coordination
    GlobToLoc=zeros(size(CA,1),6);
    for n=1:2  
        GN=CA(:,n); 
        for d=1:3 
            GlobToLoc(:,(n-1)*3+d)=(GN-1)*3+d;
        end
    end

    % Forming Global Truss Stiffness Matrix
    K = zeros(3*size(NC,1));
    for e=1:size(CA,1) 
        ke = Kbasket(:,:,e);
        for lr = 1:6
            gr = GlobToLoc(e,lr); 
            for lc = 1:6
                gc = GlobToLoc(e,lc); 
                K(gr,gc) = K(gr,gc) + ke(lr,lc);
            end
        end
    end
    
    % Sparsify Global Truss Stiffness Matrix
    K = sparse(K);
end

% FUNCTION TO CALCULATE VOLUME FRACTION
function volFrac = calcVF(NC,CA,r,sel)
    totalTrussVol = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        z1 = NC(CA(i,1),3); z2 = NC(CA(i,2),3);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2)+((z2-z1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*pi*(r^2));
    end
    % Calculating volume fraction (using a solid cube (sel) as a baseline) 
    volFrac = totalTrussVol/(sel^3);
end

% FUNCTION TO PLOT NODAL DISPLACEMENT (NOTCONVYET)
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
        if abs(C(x)) < 10^(-10)
            C(x) = 0;
        end
    end
    if abs(C(1,2)-C(2,1)) < 10^(-10)
        if abs(C(1,3)-C(3,1)) < 10^(-10)
            if abs(C(2,3)-C(3,2)) < 10^(-10)
                if (C(1,1) == C(2,2)) && (C(2,2) == C(3,3))
                    if (C(4,4) == C(5,5)) && (C(5,5) == C(6,6))
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
        else
            symm_post = false;
        end
    else
        symm_post = false;
    end
end
