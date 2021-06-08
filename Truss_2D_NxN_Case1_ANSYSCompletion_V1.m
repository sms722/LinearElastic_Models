%% 2D Truss Completion Script-- FEA/ANSYS Verification (1 UC)
%   This script is used to process raw results from the ANSYS verification 
%   simulations.  This script is specified for Case 1, 3x3 UC, 3x3 grid.

% Initialization
clc;    
close all; 
clear;

% Values from ANSYS
Fx = [1.0623e-6,7.4797e-7,1.0362e-6;
       -2.5143e-7,8.1476e-14,2.7746e-7;
       2.7747e-7,-8.5166e-23,-2.7747e-7];
Fy = [0,3.6822e-8,2.7746e-7;
       1.0362e-6,7.4797e-7,1.0623e-6;
       -2.7747e-7,1.939e-22,2.7747e-7];
Fxy = [-9.6257e-7,0,1.0623e-6;
        -2.5143e-7,8.1476e-14,2.7746e-7;
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

%% 2D Truss Completion Script-- FEA/ANSYS Verification (4 UC)
%   This script is used to process raw results from the ANSYS verification 
%   simulations.  This script is specified for Case 1, 3x3 UC, 5x5 grid.

% Initialization
clc;    
close all; 
clear;

% Values from ANSYS
Fx = [5.3211e-7,3.8932e-7,6.624e-7,3.8698e-7,5.2512e-7;
      1.3762e-7,-6.5552e-11,2.9239e-7,4.4823e-8,0;
      1.3873e-7,-6.3089e-13,-7.9372e-22,6.3089e-13,-1.3873e-7];
Fy = [1.3762e-7,-6.5552e-11,2.9239e-7,4.4823e-8,0;
      5.3211e-7,3.8932e-7,6.624e-7,3.8698e-7,5.2512e-7;
      1.3873e-7,-6.3089e-13,-2.9259e-22,6.3089e-13,-1.3873e-7];
Fxy = [5.3211e-7,0,0,0,-4.6868e-7;
       1.3762e-7,-1.0952e-13,-2.272e-10,1.046e-13,-1.3619e-7;
       1.3873e-7,2.8922e-12,2.7747e-7,2.8922e-12,1.3873e-7];
    
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

%% 2D Truss Completion Script-- FEA/ANSYS Verification (9 UC)
%   This script is used to process raw results from the ANSYS verification 
%   simulations.  This script is specified for Case 1, 3x3 UC, 7x7 grid.

% Initialization
clc;    
close all; 
clear;

% Values from ANSYS
Fx = [[],[],[],[],[],[],[];
      [],[],[],[],[],[],[];
      [],[],[],[],[],[],[]];
Fy = [[],[],[],[],[],[],[];
      [],[],[],[],[],[],[];
      [],[],[],[],[],[],[]];
Fxy = [[],[],[],[],[],[],[];
       [],[],[],[],[],[],[];
       [],[],[],[],[],[],[]];
    
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



