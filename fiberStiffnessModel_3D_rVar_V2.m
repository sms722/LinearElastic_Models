% Fiber Stiffness-based Model ("Simple Model")in 3D
% -----------------------------------------------------------------
% This model is intended to find C-matrix (material stiffness) values based
%    approximations of individual truss stiffnesses from treating them as
%    fibers at various orientations (w/ volume correction)
% -----------------------------------------------------------------
% All lengths are in [m], all stresses and moduli are in [Pa] 
% The model uses any 3D NxNxN nodal grid
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
function [C11,C22,C33,volFrac] = ...
                 fiberStiffnessModel_3D_rVar_V2(sel,rvar,E,CA,nucFac)
    % Calculated Inputs
    sidenum = (2*nucFac) + 1;
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members
             
    % Generate nodal grid
    NC = generateNC(sel,sidenum);
    
    % Find volume fraction
    volFrac = calcVF(NC,CA,rvar,sel);
    
    % Calculating C-matrix values
    C11 = fiberCalc(volFrac,NC,CA,E,1,Avar)/nucFac;
    C22 = fiberCalc(volFrac,NC,CA,E,2,Avar)/nucFac;
    C33 = fiberCalc(volFrac,NC,CA,E,3,Avar)/nucFac;

end

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

% FUNCTION TO CALCULATE C-MATRIX VALUES VIA FIBER METHOD
function Cval = fiberCalc(volFrac,NC,CA,E,dir,Avar)
    % Find effective structural stiffness
    K = E*volFrac;
    
    % Find length-corrected sum of cosines for all fibers 
    cVsum = 0;
    Vsum = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        z1 = NC(CA(i,1),3); z2 = NC(CA(i,2),3);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2)+((z2-z1)^2));
        
        % Considering directionality for length-correction
        if dir == 1
            c=(x2-x1)/L; 
        elseif dir == 2
            c=(y2-y1)/L; 
        elseif dir == 3
            c=(z2-z1)/L;
        end
        cVsum = cVsum + (L*(Avar(i))*(c^4));
        Vsum = Vsum + (L*Avar(i));
    end
    
    % Find desired C-value
    Cval = (K*cVsum)/Vsum;
    
end