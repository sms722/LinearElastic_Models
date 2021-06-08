% 3D NxNxN Cubic Truss Problem-- Specialized for "Case 1-3D" Unit Cell
% "Case 1-3D" is a 3x3x3 unit cell, symmetric in x, y, and z

% Initialization
clc;    
close all; 
clear;

UC = [1,2,3,4,5,6,7,8];
C11 = []; C12 = []; C44 = []; C55 = []; C66 = []; 
runtimes = []; CAsize = []; Ksize = [];
for nuc = 1:1:size(UC,2)
    % Start time counter
    tic
    
    % Input Constants
    nucFac = UC(nuc); % factor to indicate how many unit cells are in this 
                %   model (nucFac^3 gives the actual number of unit cells)
    sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
    r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                    %   area of (assumed circular) truss members
    E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

    % Calculated Inputs
    sidenum = (2*nucFac) + 1; % Number of nodes on each size of cube grid
    A = pi*(r^2); % Constant cross-sectional area of truss members

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum,nucFac);
    if nuc == 1
        NC1 = NC;
    elseif nuc == 2
        NC2 = NC;
    end

    % Generate "Case 2-3D" Connectivity Array
    CA = generateCA(nucFac,sidenum);
    if nuc == 1
        CA1 = CA;
    elseif nuc == 2
        CA2 = CA;
    end

    % Develop C-matrix from K-matrix (functions below)
    C = [];
    [C,uBasket,FBasket,K] = generateC(sel,NC,CA,A,E,C,nucFac);
    C = C./nucFac;
    if nuc == 1
        C_1 = C;
    elseif nuc == 2
        C_2 = C;
    end
    
    % Record outputs for plotting
    C11 = [C11,C(1,1)];
    C12 = [C12,C(1,2)];
    C44 = [C44,C(4,4)];
    C55 = [C55,C(5,5)];
    C66 = [C66,C(6,6)];
    
    % Record runtime, CA-matrix size, K-matrix size
    runtimes = [runtimes,toc]; 
    CAsize = [CAsize,size(CA,1)]; 
    Ksize = [Ksize,size(K,1)];
end

% Plot C-matrix values for increasing # of unit cells
subplot(2,2,1)
plot(UC.^3,C11,'b*',UC.^3,C12,'g*',...
     UC.^3,C44,'m*',UC.^3,C55,'k*',UC.^3,C66,'c*');
xlabel('Number of Unit Cells');
ylabel('C-Value (Pa)');
%legend('C11,C22','C12/C21','C33');

% Plot runtimes
subplot(2,2,2)
f = fit((UC.^3)',(runtimes)','exp1');
fitvals = coeffvalues(f); a = fitvals(1); b = fitvals(2);
new_UC3 = (1:1:500).^3;
fit_runtimes = a.*exp(b.*(new_UC3));
plot(UC.^3,runtimes,'b*',new_UC3,fit_runtimes,'r-');
xlabel('Number of Unit Cells');
ylabel('Runtime (Seconds)');
axis([0 ((UC(end).^3)+250) 0 (runtimes(end)+500)]);

% Plot CA sizes
subplot(2,2,3)
plot(UC.^3,CAsize,'k*');
xlabel('Number of Unit Cells');
ylabel('Size of CA');
axis([0 ((UC(end).^3)+200) 0 (CAsize(end)+100)]);

% Plot K-matrix sizes
subplot(2,2,4)
plot(UC.^3,Ksize,'m*');
xlabel('Number of Unit Cells');
ylabel('Size of K-matrix');
axis([0 ((UC(end).^3)+100) 0 (Ksize(end)+500)]);

%----------%
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum,nucFac)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            for k = 1:1:sidenum
                NC = [NC;notchvec(i),notchvec(j),notchvec(k)];
            end
        end
    end
    NC = (sel*nucFac).*NC;
end

% FUNCTION TO GENERATE "CASE 2-3D"-SPECIFIC ELEMENT CONNECTIVITY ARRAY
function CA = generateCA(nucFac,sidenum)
    CA = [];
    
    % All z-axis-parallel elements
    for i = 0:1:(sidenum-1)
        for j = 0:1:(sidenum-1)
            for k = 1:1:(sidenum-1)
                CA = [CA;((i*sidenum^2)+(j*sidenum)+k),...
                         ((i*sidenum^2)+(j*sidenum)+(k+1))];
            end
        end
    end
    
    % All y-axis-parallel elements
    for i = 0:1:(sidenum-1)
        for j = 0:1:(sidenum-2)
            for k = 1:1:(sidenum)
                CA = [CA;((i*sidenum^2)+(j*sidenum)+k),...
                         ((i*sidenum^2)+(j*sidenum)+(k+sidenum))];
            end
        end
    end
    
    % All x-axis-parallel elements
    for i = 0:1:(sidenum-2)
        for j = 0:1:(sidenum-1)
            for k = 1:1:(sidenum)
                CA = [CA;((i*sidenum^2)+(j*sidenum)+k),...
                         ((i*sidenum^2)+(j*sidenum)+k+(sidenum^2))];
            end
        end
    end
    
    % Identify face-center nodes on yz-parallel planes
    fcn_yz = [];
    for i = 0:1:(sidenum-1)
        for j = 1:1:nucFac
            for k = 1:1:nucFac
                fcn_yz = [fcn_yz,...
                    (((sidenum^2)*i)+(((2*j)-1)*sidenum)+(2*k))];
            end
        end
    end
    
    % Identify face-center nodes on xz-parallel planes
    fcn_xz = [];
    for i = 1:1:nucFac
        for j = 0:1:(sidenum-1)
            for k = 1:1:nucFac
                fcn_xz = [fcn_xz,...
                    (((sidenum^2)*((2*i)-1))+(j*sidenum)+(2*k))];
            end
        end
    end
    
    % Identify face-center nodes on xy-parallel planes
    fcn_xy = [];
    for i = 1:1:nucFac
        for j = 1:1:nucFac
            for k = 1:1:sidenum
                fcn_xy = [fcn_xy,...
                    (((sidenum^2)*((2*i)-1))+(((2*j)-1)*sidenum)+k)];
            end
        end
    end
    
    % Identify center-center nodes
    ccn = [];
    for i = 1:1:nucFac
        for j = 1:1:nucFac
            for k = 1:1:nucFac
                ccn = [ccn,...
                    ((((2*i)-1)*(sidenum^2))+(((2*j)-1)*sidenum)+(2*k))];
            end
        end
    end
    
    % Create crosses around yz-plane face-center nodes
    for q = 1:1:size(fcn_yz,2)
        CA = [CA;fcn_yz(q),fcn_yz(q)-sidenum-1;
                 fcn_yz(q),fcn_yz(q)+sidenum-1;
                 fcn_yz(q),fcn_yz(q)-sidenum+1;
                 fcn_yz(q),fcn_yz(q)+sidenum+1];
    end
    
    % Create crosses around xz-plane face-center nodes
    for q = 1:1:size(fcn_xz,2)
        CA = [CA;fcn_xz(q),fcn_xz(q)-(sidenum^2)-1;
                 fcn_xz(q),fcn_xz(q)+(sidenum^2)-1;
                 fcn_xz(q),fcn_xz(q)-(sidenum^2)+1;
                 fcn_xz(q),fcn_xz(q)+(sidenum^2)+1];
    end
    
    % Create crosses around xy-plane face-center nodes
    for q = 1:1:size(fcn_xy,2)
        CA = [CA;fcn_xy(q),fcn_xy(q)-sidenum-(sidenum^2);
                 fcn_xy(q),fcn_xy(q)-sidenum+(sidenum^2);
                 fcn_xy(q),fcn_xy(q)+sidenum-(sidenum^2);
                 fcn_xy(q),fcn_xy(q)+sidenum+(sidenum^2)];
    end
    
    % Create spikeballs around center-center nodes
    for q = 1:1:size(ccn,2)
        CA = [CA;ccn(q),ccn(q)-(sidenum^2)-sidenum-1;
                 ccn(q),ccn(q)-(sidenum^2)+sidenum-1;
                 ccn(q),ccn(q)+(sidenum^2)+sidenum-1;
                 ccn(q),ccn(q)+(sidenum^2)-sidenum-1;
                 ccn(q),ccn(q)-(sidenum^2)-sidenum+1;
                 ccn(q),ccn(q)-(sidenum^2)+sidenum+1;
                 ccn(q),ccn(q)+(sidenum^2)+sidenum+1;
                 ccn(q),ccn(q)+(sidenum^2)-sidenum+1];
    end
end

% FUNCTION TO CALCULATE C-MATRIX 
function [C,uBasket,FBasket,Kout] = generateC(sel,NC,CA,A,E,C,nucFac)
    % Initialize variables
    uBasket = []; FBasket = []; ND = NC./(nucFac*sel);
    
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
            newK = sparse([K(qvec,:);K(rvec,:)]);
            newK = sparse([newK(:,qvec),newK(:,rvec)]);
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
            newK = sparse([K(qvec,:);K(rvec,:)]);
            newK = sparse([newK(:,qvec),newK(:,rvec)]);
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
            newK = sparse([K(qvec,:);K(rvec,:)]);
            newK = sparse([newK(:,qvec),newK(:,rvec)]);
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
            newK = sparse([K(qvec,:);K(rvec,:)]);
            newK = sparse([newK(:,qvec),newK(:,rvec)]);
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
            newK = sparse([K(qvec,:);K(rvec,:)]);
            newK = sparse([newK(:,qvec),newK(:,rvec)]);
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
            newK = sparse([K(qvec,:);K(rvec,:)]);
            newK = sparse([newK(:,qvec),newK(:,rvec)]);
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
            Kout = K;
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
        if (e22 ~= 0)
            F_z = F_x;
        end
        %}
        stressvec = (1/((sel*nucFac)^2)).*[(F_x);(F_y);(F_z);...
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
