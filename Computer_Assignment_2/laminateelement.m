function [ke, fe] = laminateelement(ex, ey, ep, A, B, D, A_tilde, eq)
% -------------------------------------------------------------
%  PURPOSE
%   Calculate the stiffness matrix for a 4-node Mindlin plate element.
%
%   INPUT:  ex = [x1, x2, x3, x4]      element coordinates
%           ey = [y1, y2, y3, y4]
%
%           ep = [irb irs]          irb - integration rule bending
%                                   irs - integration rule shear
%
%           A, B, D                 constitutive matrices for CLT
%
%           A_tilde                 constitutive matrix for shear (2 x 2)
%                                   (including shear correction factor)
%
%           eq = [p]                vertical load/unit area
%
%    OUTPUT: ke :  element stiffness matrix (20 x 20)
%            fe :  equivalent nodal forces (20 x 1)
%  -------------------------------------------------------------
%
%   Element dofs: u1, v1, u2, v2, u3, v3, u4, v4, w1, w2, w3, w4,
%                 phix1, phiy1, phix2, phiy2, phix3, phiy3, phix4, phiy4
%  -------------------------------------------------------------

p = eq(1);
if length(ep) < 2
    irb = 2;
    irs = 1;
else
    irb = ep(1);
    irs = ep(2);
end

ngpb = irb^2;             % number of integration points for the bending
                          % and in-plane contributions
ngps = irs^2;             % number of gauss points for the shear
                          % contribution

% initialize workspace
Bb = zeros(3,2*4);
Bbar = zeros(2,4);
N = zeros(2,2*4);
Nbar = zeros(1,4);

ke = zeros(20,20);
fe = zeros(20,1);

udofs = 1:2*4;                      % indices for in-plane displacement
                                    % degree-of-freedoms
wdofs = udofs(end)+1:udofs(end)+4;  % indices for out-of-plane 
                                    % degree-of-freedoms
phidofs = wdofs(end)+1:20;          % indices for rotational 
                                    % degree-of-freedoms


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate the shear contributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xiet, gw] = gaussquad(ngps);     % calculating positions and weights
                                  % for integration points for the shear 
                                  % contribution
for i = 1:ngps
    xi = xiet(1, i);                  % position of integration point 
                                      % in xi-direction
    eta = xiet(2, i);                 % position of integration point in 
                                      % eta-direction

    [Nvec, dNdxi] =shapequad(xi, eta, 4);
    % Nvec = [N1 N2 N3 N4]
    % dNdxi = [dN1/dxi   dN2/dxi   dN3/dxi   dN4/dxi;
    %          dN1/deta  dN2/deta  dN3/deta  dN4/deta]

    % to calculate jacobian
    Elem_co=[ex ey];
    JT = dNdxi*Elem_co;                          % calculate JT
    detJ = det(JT);
    if (detJ < 0)
        error('Jacobian not invertable')
    end
    
    % to calculate N
    Bbar=JT\dNdxi;   % transform elemental B_tilde matrix to 
                     % global B_tilde matrix
    N= [Nvec(1,1) 0 Nvec(1,2) 0 Nvec(1,3) 0 Nvec(1,4) 0;
        0 Nvec(1,1) 0 Nvec(1,2) 0 Nvec(1,3) 0 Nvec(1,4)]; 
    
    % to calculate B
    Bb= [Bbar(1,1) 0 Bbar(1,2) 0 Bbar(1,3) 0 Bbar(1,4) 0;
         0 Bbar(2,1) 0 Bbar(2,2) 0 Bbar(2,3) 0 Bbar(2,4);
         Bbar(2,1) Bbar(1,1) Bbar(2,2) Bbar(1,2) Bbar(2,3)...
         Bbar(1,3) Bbar(2,4) Bbar(1,4)];
    
    Wt=detJ*gw(i);                 
    Kphiphi = 5/6*N'*A_tilde*N;    
    Kphiw = 5/6*N'*A_tilde*Bbar;   
    Kwphi = 5/6*Bbar'*A_tilde*N;                         
    Kww = 5/6*Bbar'*A_tilde*Bbar;                           
    
    ke(phidofs, phidofs) = ke(phidofs, phidofs) + Kphiphi * Wt;
    ke(phidofs, wdofs) = ke(phidofs, wdofs)     + Kphiw   * Wt;
    ke(wdofs, phidofs) = ke(wdofs, phidofs)     + Kwphi   * Wt;
    ke(wdofs, wdofs) = ke(wdofs, wdofs)         + Kww     * Wt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate the remaining parts of stiffness and the force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xiet, gw] = gaussquad(ngpb);
for i = 1:ngpb
    xi = xiet(1, i);
    eta = xiet(2, i);

    [Nvec, dNdxi] = shapequad(xi, eta, 4);
    % Nvec = [N1 N2 N3 N4]
    % dNdxi = [dN1/dxi   dN2/dxi   dN3/dxi   dN4/dxi;
    %          dN1/deta  dN2/deta  dN3/deta  dN4/deta]
    
    % to calculate jacobian
    Elem_co=[ex ey];          % assign co-ordinate to a matrix to extract
                              % displacement components in next section
    JT = dNdxi*Elem_co;       % calculate Jacobian
    detJ = det(JT);           % determinant of Jacobian
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
    
    
    Wt = detJ * gw(i);                      
    Kuu = Bb'*A*Bb;           
    Kuphi = Bb'*B*Bb;    
    Kphiu = Bb'*B*Bb;   
    Kphiphi = Bb'*D*Bb; 

    ke(udofs, udofs)    = ke(udofs, udofs)    + Kuu     * Wt;
    ke(udofs, phidofs)  = ke(udofs, phidofs)  + Kuphi   * Wt;
    ke(phidofs, udofs)  = ke(phidofs, udofs)  + Kphiu   * Wt;
    ke(phidofs,phidofs) = ke(phidofs,phidofs) + Kphiphi * Wt;

    % computing force due to vertical load p
    if nargout > 1
        fe(wdofs) = fe(wdofs) + (p * Wt) * Nvec';      
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of element routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper-functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N, B] = shapequad(xi, eta, nn)
    if nn == 4
        N = 1/4*[(xi-1)*(eta-1), -(xi+1)*(eta-1), (xi+1)*(eta+1),...
            -(xi-1)*(eta+1)];
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
