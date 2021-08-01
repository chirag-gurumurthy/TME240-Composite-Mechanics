function [sigL, sigT, tauLT, zcoord] = laminatestress(ex, ey, h, angles, Q, ed)
% -------------------------------------------------------------
%  PURPOSE: None.
%
%   INPUT:  ex  = [x1, x2, x3, x4] element coordinates
%           ey  = [y1, y2, y3, y4]
%           h   = [-hlam/2; .... hlam/2] z-coordinates of interfaces
%           ang = [ang1; ang2; ang3; ... angn] lamina directions
%           Q   = [E_L/(1+\nu_T*nu_TL) ...] material stiffness matrix 
%                 (in local coordinate system)
%           ed  = [u1; v1; u2; v2; u3; v3; u4; v4;
%                 w1; w2; w3; w4;
%                 phix1; phiy1; phix2; phiy2; phix3; phiy3; phix4; phiy4] element displacements
%
%   OUTPUT: In plane stresses in the *middle of the element* as
%          sigL = [sigL_h_min ... sigL_hmax]      Note! Two values per ply
%          sigT = [sigT_h_min ... sigT_hmax]      Note! Two values per ply
%         tauLT = [tauLT_h_min ... tauLT_hmax]   Note! Two values per ply
%        zcoord = [h0 h1 h1 h2 h2 h3 h3 ... hn]
%
% Fredrik Ekre, 2016, 2018

[xiet] = gaussquad(4);      % calculating positions and weights for 
                            % integration points for the shear contribution
p=0;
for j = 1:length(angles)
for o=0:1
    p=p+1;
for i = 1:4
    xi = xiet(1, i);       % position of integration point in xi-direction
    eta = xiet(2, i);      % position of integration point in eta-direction

    [Nvec, dNdxi] =shapequad(xi, eta, 4);
    % Nvec = [N1 N2 N3 N4]
    % dNdxi = [dN1/dxi   dN2/dxi   dN3/dxi   dN4/dxi;
    %          dN1/deta  dN2/deta  dN3/deta  dN4/deta]

    % to calculate jacobian
    Elem_co=[ex ey];
    JT = dNdxi*Elem_co; % calculate JT
    detJ = det(JT);
    if (detJ < 0)
        error('Jacobian not invertable')
    end
    
    % to calculate N
    Bbar=JT\dNdxi;
    N= [Nvec(1,1) 0 Nvec(1,2) 0 Nvec(1,3) 0 Nvec(1,4) 0;
        0 Nvec(1,1) 0 Nvec(1,2) 0 Nvec(1,3) 0 Nvec(1,4)]; 
    
    % to calculate B
    Bb= [Bbar(1,1) 0 Bbar(1,2) 0 Bbar(1,3) 0 Bbar(1,4) 0;
         0 Bbar(2,1) 0 Bbar(2,2) 0 Bbar(2,3) 0 Bbar(2,4);
         Bbar(2,1) Bbar(1,1) Bbar(2,2) Bbar(1,2) Bbar(2,3) ...
         Bbar(1,3) Bbar(2,4) Bbar(1,4)]; 
     
    % MIDPLANE
    % strain
    u      = ed(1:8,1);   % assign u1; v1; u2; v2; u3; v3; u4; v4 to a 
                          % matrix 'u' which is used to calculate strain 
                          % in the middle element of the plate
    eps_0  = Bb*u;        %strain in midplane
    
    % curvature
    phi    = ed(13:20,1); % assign phix1; phiy1; phix2; phiy2; phix3; phiy3;
                          % phix4; phiy4 to a matrix 'phi' which is used to
                          % calculate curvature in the middle element of 
                          % the plate
    kurv   = Bb*phi;    %curvature
     

     
    zcoord(p)=h(j+o);
    epsilon_global = eps_0 + zcoord(p)*kurv;
     
    theta          = angles(j);
    %%%%%%%%%%%%%%% transform T_2 %%%%%%%%%%%%%%%%%%
    T2 = [cosd(theta)^2 sind(theta)^2 sind(theta)*cosd(theta);
    sind(theta)^2 cosd(theta)^2 -sind(theta)*cosd(theta);
    -2*sind(theta)*cosd(theta) 2*sind(theta)*cosd(theta)...
    cosd(theta)^2-sind(theta)^2];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sigma_local     = Q*T2*epsilon_global;
     
    sigL_local(i)   = sigma_local(1,:);
    sigT_local(i)   = sigma_local(2,:);
    tauLT_local(i)  = sigma_local(3,:);
     
end 
    
    sigL(p)  = mean(sigL_local);
    sigT(p)  = mean(sigT_local);
    tauLT(p) = mean(tauLT_local);
end

end
 
end

function [N, B] = shapequad(xi, eta, nn)
    if nn == 4
        N = 1/4*[(xi-1)*(eta-1), -(xi+1)*(eta-1), (xi+1)*(eta+1), -(xi-1)*(eta+1)];
        B = 1/4*[eta-1, -(eta-1), eta+1, -(eta+1);
            xi-1, -(xi+1), xi+1, -(xi-1)];
    else
        error('Only 4-noded quadrilateral implemented')
    end
end

function [xiet, w] = gaussquad(ngauss)
    if ngauss == 1
        xiet = [0; 0];
        w = [4];
    elseif ngauss == 4
        xiet = [-1  1  1 -1;
                -1 -1  1  1] / sqrt(3);
        w = [1 1 1 1];
    else
        error('Only 1x1 and 2x2 integration implemented')
    end
end





