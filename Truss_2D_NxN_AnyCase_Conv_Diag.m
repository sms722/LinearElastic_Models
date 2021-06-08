% 2D NxN Square Truss Problem-- Specialized for "Case 1" Unit Cell
% "Case 1" is a 3x3 unit cell, symmetric in both x and y, with all possible
% positions for elements filled in, and a unit cell connectivity of 8 at
% the central node (and unit cell connectivity of 3 at all other nodes)

% Initialization
clc;   
close all; 
clear;

UC = [1,2,3,4,6,8,10]; UC2 = ([1,2,3,4,6,8,10].^2);
C11basket = []; C12basket = []; C44basket = [];
bC11 = []; bC12 = []; bC44 = []; dC11 = []; dC12 = []; dC44 = [];
fitObjs = {};

% Unit Cell CAs-- Change as needed
CABasket(1) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,6;2,7;5,7;6,7]};

CABasket(2) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;...
                5,7;5,8]};
CABasket(3) = {[2,3;3,6;8,9;1,4;1,5;1,6;2,6;4,5;5,6;5,7;6,8]};
CABasket(4) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;4,5;5,6;2,5;5,8;5,7;2,6]};
CABasket(5) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;4,5;5,6;3,4;3,5;6,7;...
                6,8]};

CABasket(6) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;4,5;5,6;2,5;1,5;4,9]};
            %{
CABasket(7) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,5;2,6;5,6;5,9;6,8]};
CABasket(8) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;2,6;2,4;...
                4,8]};
CABasket(9) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;4,5;5,7;5,9]};
CABasket(10) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,6;4,5;5,7;6,8]};
CABasket(11) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,6;5,6;5,7;5,9]};
CABasket(12) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;3,5;5,6;5,7;5,8;5,9]};
CABasket(13) = {[1,2;2,3;3,6;8,9;7,8;1,4;1,5;3,5;5,6;5,7;5,8;5,9]};
CABasket(14) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,6;6,8;1,5;4,5;5,7]};
CABasket(15) = {[1,2;6,9;7,8;4,7;1,5;2,5;3,5;5,6;5,7;5,9]};
CABasket(16) = {[2,3;3,6;6,9;7,8;4,7;1,4;1,5;2,6;4,5;5,6;4,8;5,9]};
CABasket(17) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;5,8;4,8]};
CABasket(18) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;...
                 5,7;5,8;5,9]};
CABasket(19) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,5;2,6;4,5;5,6;...
                 4,8;5,8;6,8]};
CABasket(20) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,5;2,6;4,5;5,7;6,8]};
CABasket(21) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,5;3,5;4,5;5,6;6,7]};
CABasket(22) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;3,4;3,5;4,5;5,6;6,7]};
CABasket(23) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;3,4;3,5;4,5;5,6;...
                 5,7;6,7;6,8]};
CABasket(24) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,7;2,9;5,7;5,9;5,8;...
                 2,4;2,6]};
CABasket(25) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;1,6;5,6;5,7;5,8;5,9]};
CABasket(26) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,6;2,9;4,8;4,9;...
                 2,5;4,5;5,9]};
CABasket(27) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,7;2,5;3,5;5,6;5,7;...
                 5,8;5,9]};
CABasket(28) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;3,4;3,5;5,6;4,8;5,8;6,8]};
CABasket(29) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,7;2,6;5,6;5,8;5,9]};
CABasket(30) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,7;2,5;3,5;3,8;...
                 5,7;5,8;6,8]};
%}
% Loop through each design to solve for, plot C-matrix values
for cac = 1:1:size(CABasket,2)
    tinyCA = cell2mat(CABasket(cac));
    %disp(tinyCA);
    C11 = []; C12 = []; C44 = [];

    % Loop through each # of unit cells
    for nuc = 1:1:size(UC,2)
        % Input Constants
        nucFac = UC(nuc); % factor to indicate how many unit cells are in  
              % this model (nucFac^2 gives the actual number of unit cells)
        sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
        r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                         %  area of (assumed circular) truss members
        E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

        % Calculate C, volFrac for each # of unit cells case

        [C,volFrac] = calcCTruss(nucFac,sel,r,E,tinyCA);

        % Record outputs for plotting/reference
        C11 = [C11,C(1,1)];
        C12 = [C12,C(1,2)];
        C44 = [C44,C(3,3)];
    end
    C11basket = [C11basket;C11]; 
    %C12basket = [C12basket;C12]; 
    %C44basket = [C44basket;C44]; 

    fC11 = fit(UC',C11','exp2'); fitObjs(cac) = {fC11};
    %fC12 = fit(UC',C12','exp2');
    %fC44 = fit(UC',C44','exp2');
    %{
    C11vals = coeffvalues(fC11);
    %C12vals = coeffvalues(fC12);
    %C44vals = coeffvalues(fC44);

    bC11 = [bC11;(-1)/C11vals(2)]; 
    dC11 = [dC11;(-1)/C11vals(4)];
    %bC12 = [bC12;(-1)/C12vals(2)]; dC12 = [dC12;(-1)/C12vals(4)];
    %bC44 = [bC44;(-1)/C44vals(2)]; dC44 = [dC44;(-1)/C44vals(4)];

    % Plot results (fits vs. points)

    %figure(1);
    subplot(3,2,cac)
    plot(fC11,UC,C11)
    title(['Design ',num2str(cac)])
    xlabel('Number of unit cells per side');
    ylabel('C11');
    hold on;
    %}
    %{
    figure(2);
    subplot(3,2,cac)
    plot(fC12,UC2,C12)
    ttl = ['Design ',cac];
    title(ttl)
    hold on;

    figure(3);
    subplot(3,2,cac)
    plot(fC44,UC2,C44)
    title(['Design ',cac])
    hold on;
    %}
end
%{
% Plot results (boxplot of C11, C12, C44)
bvals = [bC11, bC12, bC44];
groups = ['C11';'C12';'C44'];
boxplot(bvals,groups);
xlabel('C-Matrix Value');
ylabel('Convergence Constant');
%}


%----------%
% FUNCTION TO CALCULATE FULL CA, THEN C, FOR GIVEN # OF UNIT CELLS
function [C,volFrac] = calcCTruss(nucFac,sel,r,E,tinyCA)
    % Calculated Inputs
    sidenum = (2*nucFac) + 1; % Number of nodes on each size of square grid
    A = pi*(r^2); % Constant cross-sectional area of truss members

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);

    % Generate Connectivity Array for given # of UC's
    CA = generateCA(sidenum,tinyCA);

    % Develop C-matrix from K-matrix (functions below)
    C = [];
    [C,~,~] = generateC(sel,r,NC,CA,A,E,C);
    C = C./nucFac;

    % Calculate volume fraction (function below)
    volFrac = calcVF(NC,CA,r,sel);
end

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

% FUNCTION TO GENERATE CONNECTIVITY ARRAY FOR 3x3 UC'S, nucFac x nucFac 
function CA = generateCA(sidenum,tinyCA)
    CA = [];

    % Outside frame
    for i = 1:1:(sidenum-1)
        CA = [CA;i,(i+1)];
    end
    for i = sidenum:sidenum:((sidenum^2)-sidenum)
        CA = [CA;i,(i+sidenum)];
    end
    for i = ((sidenum^2)-sidenum+1):1:((sidenum^2)-1)
        CA = [CA;i,(i+1)];
    end
    for i = 1:sidenum:((sidenum^2)-(2*sidenum)+1)
        CA = [CA;i,(i+sidenum)];
    end
    
    % Identify central node for all unit cells
    cnvec = [];
    for j = (sidenum+2):(2*sidenum):((sidenum^2)-(2*sidenum)+2)
        cnvec = [cnvec;j];
        for i = 2:2:(sidenum-2)
            cnvec = [cnvec;(j+i)];
        end
    end
    
    % Populate CA using template of tinyCA and locations of cn's
    for w = 1:1:size(cnvec,1)
        cn = cnvec(w);
        for q = 1:1:size(tinyCA,1)
            row = [];
            for r = 1:1:2
                if tinyCA(q,r) == 1
                    row = [row,(cn-sidenum-1)];
                elseif tinyCA(q,r) == 2
                    row = [row,(cn-sidenum)];
                elseif tinyCA(q,r) == 3
                    row = [row,(cn-sidenum+1)];
                elseif tinyCA(q,r) == 4
                    row = [row,(cn-1)];
                elseif tinyCA(q,r) == 5
                    row = [row,cn];
                elseif tinyCA(q,r) == 6
                    row = [row,(cn+1)];
                elseif tinyCA(q,r) == 7
                    row = [row,(cn+sidenum-1)];
                elseif tinyCA(q,r) == 8
                    row = [row,(cn+sidenum)];
                elseif tinyCA(q,r) == 9
                    row = [row,(cn+sidenum+1)];
                end
            end
            CA = [CA;row];
        end
    end
    
    % Eliminate duplicate members
    CA = unique(CA,'rows');
end

% FUNCTION TO CALCULATE C-MATRIX
function [C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C)
    % Initialize outputs
    uBasket = []; FBasket = [];

    % Iterate through once for each strain component
    for y = 1:1:3
      % Define vectors to hold indexes for output forces
        Fi_x = []; Fi_y = []; Fi_xy = [];

      % Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

      % set that component equal to a dummy value (0.01 strain), 
      %     set all other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;

      % use strain relations, BCs, and partitioned K-matrix to 
      %     solve for all unknowns
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

      % use system-wide F matrix and force relations to solve for 
      %     stress vector [s11,s22,s12]'
        F_x = 0; F_y = 0; F_xy = 0;
        for n = 1:1:size(Fi_xy,2)
            F_x = F_x + F(Fi_x(n));
            F_y = F_y + F(Fi_y(n));
            F_xy = F_xy + F(Fi_xy(n));
        end
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

      % use strain and stress vectors to solve for the corresponding
      %     row of the C matrix
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



