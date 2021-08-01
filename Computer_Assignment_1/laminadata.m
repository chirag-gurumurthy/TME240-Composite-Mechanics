function [Q, alpha, E_T] = laminadata( E_f, E_m, nu_f, nu_m, V_f, V_m,...
    alpha_f, alpha_m, xsi_E, xsi_G )
%LAMINADATA CREATES THE NESSESSARY MATERIAL DATA FOR A LAMINA FROM FIBRE
%AND MATRIX MATERIAL DATA
%   Input:
%     E_f           -       Young's modulus for fibre
%     E_m           -       Young's modulus for matrix
%     nu_f          -       Poisson's ratio for fibre
%     nu_m          -       Poisson's ratio for matrix
%     V_f           -       Volumefraction for fibre
%     V_m           -       Volumefraction for matrix
%     alpha_f       -       Thermal expansion coeficient for fibre
%     alpha_m       -       Thermal expansion coeficient for matrix
%     xsi_E         -       xsi term in Halpin-Tsai for Young's modulus
%     xsi_G         -       xsi term in Halpin-Tsai for Shear modulus
% 
%   Output:
%     Q             -       3x3 Q matrix of a laminate
%     alpha         -       3x1 vector with the thermal expansion
%                           coeficients in local coordinates
% 

% ===== RESULTS: =====
% === E_T === (using ROM)
E_L     = E_f*V_f+E_m*V_m;

% === E_T === (using Halpin-Tsai)
eta_E   = (E_f/E_m - 1)/(E_f/E_m + xsi_E);
E_T     = E_m*(1+xsi_E*eta_E*V_f)/(1-eta_E*V_f);

% === nu ===
nu_LT   = nu_f*V_f+nu_m*V_m; %Major
nu_TL   = nu_LT*E_T/E_L; %Minor

% === G_TL === (Using Halpin-Tsai)
G_f     = E_f/(2+2*nu_f);
G_m     = E_m/(2+2*nu_m);

eta_G   = (G_f/G_m - 1)/(G_f/G_m + xsi_G);
G_LT    = G_m*(1+xsi_G*eta_G*V_f)/(1-eta_G*V_f);

% === Q ===
Q       = [E_L/(1-nu_LT*nu_TL) nu_TL*E_L/(1-nu_LT*nu_TL) 0;...
           nu_LT*E_T/(1-nu_LT*nu_TL) E_T/(1-nu_LT*nu_TL) 0;...
           0 0 G_LT]; 

% === alpha ===       
alpha_L = 1/E_L*(alpha_f*E_f*V_f + alpha_m*E_m*V_m);
alpha_T = (1 + nu_f)*alpha_f*V_f + (1 + nu_m)*alpha_m*V_m - alpha_L*nu_LT;

alpha=[alpha_L alpha_T 0]';

end

