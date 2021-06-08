% Volume Fraction-based model ("Really Dumb Model")
% This model is intended to find effective material constants (only C11,
%    C22) based on the volume fraction of a given unit cell design
% All lengths are in [m], all stresses and moduli are in [Pa] 
% The model uses any 2D NxN nodal grid
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
function [Cval,volFrac]=volFracModel(sel,r,E,CA,sidenum)
    % Generate nodal grid
    NC = generateNC(sel,sidenum);
    
    % Find volume fraction
    volFrac = calcVF(NC,CA,r,sel);
    
    % Find C11/C22
    Cval = volFrac*E;
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