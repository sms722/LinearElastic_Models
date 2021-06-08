%% 2D Truss Completion Script-- FEA/ANSYS Verification (1 UC)
%   This script is used to process raw results from the ANSYS verification 
%   simulations.  This script is specified for Case 1, 3x3 UC, 3x3 grid.

% Initialization
clc;    
close all; 
clear;

% Values from ANSYS
Fx = [1.0623e-6,7.8479e-7,1.0623e-6;
      2.7746e-7,-1.8424e-12,2.7746e-7;
      2.7747e-7,-7.6268e-22,-2.7747e-7];
Fy = [2.7746e-7,-1.8424e-12,2.7746e-7;
      1.0623e-6,7.8479e-7,1.0623e-6;
      -2.7747e-7,2.0667e-22,2.7747e-7];
Fxy = [-1.0623e-6,0,1.0623e-6;
       -2.7746e-7,-1.0296e-22,2.7746e-7;
       2.7747e-7,5.878e-12,2.7747e-7];
    
% Input Constants
sel = 0.05;
r = 50*(10^-6);

% Calculating C-matrix
C = [];
for y = 1:1:3
    strainvec = [0;0;0];
    strainvec(y) = 0.01; 
    strainvec(3) = strainvec(3)*2;
    F_x = sum(Fx(y,:));
    F_y = sum(Fy(y,:));
    F_xy = sum(Fxy(y,:));
    stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];
    Cdummy = stressvec/strainvec;
    C(:,y) = Cdummy(:,y);
end

% Print C-matrix as output
disp('The C-matrix is: '); disp(C);

%% 2D Truss Completion Script-- ANSYS M-APDL Verification (1 UC)
%   This script is used to process raw results from the APDL verification 
%   simulations.  This script is specified for Case 1, 3x3 UC, 3x3 grid.

% Initialization
clc;    
close all; 
clear;

% Values from ANSYS
Fx = [1.0623e-6,7.8479e-7,1.0623e-6;
      2.7746e-7,-1.8424e-12,2.7746e-7;
      2.7747e-7,-7.6268e-22,-2.7747e-7];
Fy = [2.7746e-7,-1.8424e-12,2.7746e-7;
      1.0623e-6,7.8479e-7,1.0623e-6;
      -2.7747e-7,2.0667e-22,2.7747e-7];
Fxy = [-1.0623e-6,0,1.0623e-6;
       -2.7746e-7,-1.0296e-22,2.7746e-7;
       2.7747e-7,5.878e-12,2.7747e-7];
    
% Input Constants
sel = 0.05;
r = 50*(10^-6);

% Calculating C-matrix
C = [];
for y = 1:1:3
    strainvec = [0;0;0];
    strainvec(y) = 0.01; 
    strainvec(3) = strainvec(3)*2;
    F_x = sum(Fx(y,:));
    F_y = sum(Fy(y,:));
    F_xy = sum(Fxy(y,:));
    stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];
    Cdummy = stressvec/strainvec;
    C(:,y) = Cdummy(:,y);
end

% Print C-matrix as output
disp('The C-matrix is: '); disp(C);

