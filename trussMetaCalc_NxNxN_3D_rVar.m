% -----------------------------------------------------------------
% Function for calculating metamaterial properties for a 3D NxNxN truss 
% with individually-adjustable element radii
% -----------------------------------------------------------------
% All lengths are in [m], all stresses and moduli are in [Pa] 
% Each node has 3 degrees of freedom (DOF), x,y, and z.  The 
% (node number*3) gives the number for that node's zDOF, (zDOF - 1) is the
% yDOF, and (zDOF - 2) is the xDOF.
% -----------------------------------------------------------------
% Sample values to test code (copy-paste into command window):
%{
clc;    
close all; 
clear;
nucFac = 1; 
sel = 0.05; 
E = 10000;
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;
      10,11;11,12;12,15;15,18;17,18;16,17;13,16;10,13;
      11,14;14,17;13,14;14,15;
      19,20;20,21;21,24;24,27;26,27;25,26;22,25;19,22;
      20,23;23,26;22,23;23,24;
      1,10;10,19;2,11;11,20;3,12;12,21;6,15;15,24;
      9,18;18,27;8,17;17,26;7,16;16,25;4,13;13,22;5,14;14,23;
      1,14;3,14;7,14;9,14;14,19;14,21;14,25;14,27];
rvar = (50*(10^-6)).*ones(1,size(CA,1)); 
%}

function [C,volFrac] = trussMetaCalc_NxNxN_3D_rVar(nucFac,sel,rvar,E,CA)
    
    % Calculated Inputs
    sidenum = (2*nucFac) + 1; % Number of nodes on each size of square grid
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members
    
    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Develop C-matrix from K-matrix (functions below)
    C = [];
    [C,uBasket,FBasket] = generateC(sel,NC,CA,Avar,E,C);
    C = C./nucFac;

    % Print C-matrix as output
    %disp('The C-matrix is: '); disp(C);
    
    % Calculate and print volume fraction (function below)
    volFrac = calcVF(NC,CA,rvar,sel);
    %disp('The volume fraction is: '); disp(volFrac);
end

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

% FUNCTION TO CALCULATE C-MATRIX 
function [C,uBasket,FBasket] = generateC(sel,NC,CA,Avar,E,C)
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
            K = formK(NC,CA,Avar,E); % function for this below
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
            K = formK(NC,CA,Avar,E); % function for this below
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
            K = formK(NC,CA,Avar,E); % function for this below
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
            K = formK(NC,CA,Avar,E); % function for this below
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
            K = formK(NC,CA,Avar,E); % function for this below
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
            K = formK(NC,CA,Avar,E); % function for this below
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
function K = formK(NC,CA,Avar,E)
    % Forming Elemental Stiffness Matrices
    Kbasket = [];
    for i = 1:size(CA,1)
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        z1 = NC(CA(i,1),3); z2 = NC(CA(i,2),3);
        
        % Plotting elements (comment out if undesired)
        %line([x1,x2],[y1,y2],[z1,z2]);
        
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
        ke = ((Avar(i).*E)./L).*ktemp;
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
function volFrac = calcVF(NC,CA,rvar,sel)
    totalTrussVol = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        z1 = NC(CA(i,1),3); z2 = NC(CA(i,2),3);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2)+((z2-z1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*pi*((rvar(i))^2));
    end
    % Calculating volume fraction (using a solid cube (sel) as a baseline) 
    volFrac = totalTrussVol/(sel^3);
end


