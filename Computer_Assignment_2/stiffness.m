function [T1,T2,Q,T_1tilde,Q_tilde] = stiffness(angles, E_L, E_T, G_LT, G_TT, nu_LT, nu_TL)

%input parameters
%angles     = fibre orientation
%angles_r   = angles*pi/180; for angles to be in radians

%output parameters
T1=[cosd(angles)^2 sind(angles)^2 2*sind(angles)*cosd(angles);
    sind(angles)^2 cosd(angles)^2 -2*sind(angles)*cosd(angles);
    -sind(angles)*cosd(angles) ...
    sind(angles)*cosd(angles) cosd(angles)^2-sind(angles)^2];

T2=[cosd(angles)^2 sind(angles)^2 sind(angles)*cosd(angles);
    sind(angles)^2 cosd(angles)^2 -sind(angles)*cosd(angles);
    -2*sind(angles)*cosd(angles) ...
    2*sind(angles)*cosd(angles) cosd(angles)^2-sind(angles)^2];


 Q=[(E_L/(1-(nu_LT*nu_TL))) (nu_TL*E_L/(1-(nu_LT*nu_TL))) 0;
     (nu_LT*E_T/(1-(nu_LT*nu_TL))) (E_T/(1-(nu_LT*nu_TL))) 0;
     0 0 G_LT];
 
 T_1tilde= [cosd(angles) -sind(angles);
            sind(angles)  cosd(angles)];
 
 Q_tilde=  inv(T_1tilde)*[G_TT  0;
                           0    G_LT]*T_1tilde;   
               
end
