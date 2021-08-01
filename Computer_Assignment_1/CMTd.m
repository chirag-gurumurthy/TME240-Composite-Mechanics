function [ T1 , T2] = CMTd(theta)
%CMTd calculates the T1 and T2 transformation matrices using theta in
%degrees
%
%   INPUT
%   [theta]      -       angle between global x-direction and the laminar
%                        direction. Theta should be given in degrees.
%
%   NB! - If the angle theta is given in radians use the function: 
%  [ T1 , T2] = CMTr(theta) 

T1=[cosd(theta)^2 sind(theta)^2 2*sind(theta)*cosd(theta);
    sind(theta)^2 cosd(theta)^2 -2*sind(theta)*cosd(theta);
    -sind(theta)*cosd(theta) sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2];

T2=[cosd(theta)^2 sind(theta)^2 sind(theta)*cosd(theta);
    sind(theta)^2 cosd(theta)^2 -sind(theta)*cosd(theta);
    -2*sind(theta)*cosd(theta) 2*sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2];

end
