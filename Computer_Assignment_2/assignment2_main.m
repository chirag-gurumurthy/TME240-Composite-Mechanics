clear all;
close all;
clc;
% Geometry
Lx = 0.4;                  % length in x-direction
Ly = 0.6;                  % length in y-direction
nx = 21;                   % number of elements in x-direction
ny = 21;                   % number of elements in y-direction
thickness = 0.001;         % The total thickness of the plate
nl = 10;                   % no. of layers

angles = [0 0 45 90 -45 -45 90 45 0 0];     % Orientation of the fibres

E_L = 181e9;            % Youngs modulus in longitudinal direction[Pa]
E_T = 10.3e9;           % Youngs modulus in transverse direction[Pa]
nu_LT = 0.28;           % poissons ratio in LT direction[-]
G_LT = 7.17e9;          % shear modulus in LT direction[Pa]
G_TT = 3.5e9;           % shear modulus in TT or TL direction[Pa]
nu_TL = nu_LT*E_T/E_L;  % poissons ratio in TL direction

%Q matrix
 Q=[(E_L/(1-(nu_LT*nu_TL))) (nu_TL*E_L/(1-(nu_LT*nu_TL))) 0;
     (nu_LT*E_T/(1-(nu_LT*nu_TL))) (E_T/(1-(nu_LT*nu_TL))) 0;
     0 0 G_LT];
 
% Applied load [Pa]
p = -5e3;

% Type of boundary condition
bc = 'clamped';
%bc = 'simplysupported'

% generate the FE mesh (do not use too many elements because this is slow!)
[Enod, Coord, Edof] = rectangularmesh(0, Lx, 0, Ly, nx, ny);

% Enode: for each column i (corresponding to element i), the rows j=1-4 
%        give the
% nodes
%
%           elem_1   elem_2   ...     elem_nel
% Enode = [   n1       n1     ...        n1
%             n2       n2     ...        n2
%             n3       n3     ...        n3
%             n4       n4     ...        n4];

% Coord: Each row i gives the x and y coordinate of global node i
% Coord =   [x1 y1
%            x2 y2
%            .  .
%            .  .
%            xnno ynno];

% Edof: Each column corresponds to the degrees of freedom for an element
%       such that Edof(:,i) contains the dofs for the ith element in the 
%       following order:
%       u1, v1, u2, v2, u3, v3, u4, v4, w1, w2, w3, w4,
%       phix1, phiy1, phix2, phiy2, phix3, phiy3, phix4, phiy4

nel = length(Edof(1, :));   %number of elements
nno = length(Coord(:, 1));  %number of nodes

%close all
figure(1);
hold on;

% plot the undeformed mesh
for i = 1:nel
    plot(Coord(Enod([1 2 3 4 1], i), 1), Coord(Enod([1 2 3 4 1], i), 2));
end

% Generate A, B, D, A_tilde matrices
[A, B, D, A_tilde] = foo(E_L, E_T, nu_LT,nu_TL, G_LT, G_TT, angles,...
    thickness);

h=zeros(1,length(angles)+1);
m=1;
for n=thickness/2:-thickness/length(angles):-thickness/2
    h(m)=n;
    m=m+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes1 = find(Coord(:,1)==min(Coord(:,1)));
nodes2 = find(Coord(:,2)==min(Coord(:,2)));
nodes3 = find(Coord(:,1)==max(Coord(:,1)));
nodes4 = find(Coord(:,2)==max(Coord(:,2)));

% degrees-of-freedom on left most boundary x = 0
cdofx1 = nodes1*5-4;
cdofy1 = nodes1*5-3;
cdofz1 = nodes1*5-2;
cdofphix1 = nodes1*5-1;
cdofphiy1 = nodes1*5;

% degrees-of-freedom on lower boundary y = 0
cdofx2 = nodes2*5-4;
cdofy2 = nodes2*5-3;
cdofz2 = nodes2*5-2;
cdofphix2 = nodes2*5-1;
cdofphiy2 = nodes2*5;

% degrees-of-freedom on right most boundary x = Lx
cdofx3 = nodes3*5-4;
cdofy3 = nodes3*5-3;
cdofz3 = nodes3*5-2;
cdofphix3 = nodes3*5-1;
cdofphiy3 = nodes3*5;

% degrees-of-freedom on upper boundary y = Ly
cdofx4 = nodes4*5-4;
cdofy4 = nodes4*5-3;
cdofz4 = nodes4*5-2;
cdofphix4 = nodes4*5-1;
cdofphiy4 = nodes4*5;


switch (bc)
    case 'clamped'
        % Example: clamped plate
        cdof = [cdofx1;cdofx2;cdofx3;cdofx4;cdofy1;cdofy2;cdofy3;...
            cdofy4;cdofz1;cdofz2;cdofz3;cdofz4;...
            cdofphix1;cdofphix2;cdofphix3;cdofphix4;cdofphiy1;...
            cdofphiy2;cdofphiy3;cdofphiy4];
        cdof = intersect(cdof,cdof);
        cdof = [cdof zeros(length(cdof),1)];
    case 'simplysupported'
        % % Example: simply supported plate
        cdof = [cdofx1;cdofx2;cdofx3;cdofx4;cdofy1;cdofy2;cdofy3;...
            cdofy4;cdofz1;cdofz2;cdofz3;cdofz4];
        cdof = intersect(cdof,cdof);
        cdof = [cdof zeros(length(cdof),1)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END Define BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initialize global stiffness matrix and force vector
K = spalloc(5*nno, 5*nno, 20*5*nno);
f = zeros(5*nno,1);

% Assemble the element contributions
for i=1:nel
    edof = Edof(:, i);

    [ke, fe] = laminateelement(Coord(Enod(:, i), 1), Coord(Enod(:, i),...
        2), [2, 1], A, B, D, A_tilde, p);

    K(edof, edof) = K(edof, edof) + ke;
    f(edof) = f(edof) + fe;
end

% Invoke the boundary conditions and solve
free = setdiff(1:(5*nno), cdof(:,1)');
a = zeros(5*nno, 1);
a(free)=K(free, free) \ f(free);

% Plot the resulting deformations
figure(1);
hold on;
s = max(abs(a(3:5:end)));
s = 0.1/s;

for i=1:nel
    plot3(Coord(Enod([1 2 3 4 1],i),1),Coord(Enod([1 2 3 4 1],i),2),...
        s*a(Edof([9:12,9], i)),'r')
end

title(sprintf('displacement scale factor is %g, max deflection is %g, and error is %g [%%]',...
    s, max(abs(a(3:5:end))), (0.026943-max(abs(a(3:5:end))))/0.026943*100))
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtask 2: calculate stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to extract 
    [nie,n]=size(Edof);
    for i = 1:nie
        ed(i,1:n)=a(Edof(i,:))';
    end
 
    
% plot stress for middle element
middle_element = ceil(n/2);     % evaluating middle element of the plate 
dismid= ed(:,middle_element);   % assigning displaccement of middle 
                                % element of the plate
ex=Coord(Enod(:, nel/2+0.5), 1);
ey=Coord(Enod(:, nel/2+0.5), 2);
[sigL, sigT, tauLT, zCoord] = laminatestress(ex, ey, h, angles, Q, ed);
figure
subplot(1, 3, 1)
plot(sigL, zCoord)
legend('\sigma_L')
subplot(1, 3, 2)
plot(sigT, zCoord)
legend('\sigma_T')
subplot(1, 3, 3)
plot(tauLT, zCoord)
legend('\sigma_{LT}')

% plot stress on domain
sigL = zeros(length(sigL), nel);
sigT = zeros(length(sigT), nel);
tauLT = zeros(length(tauLT), nel);

for i = 1:nel
    ex=Coord(Enod(:, i), 1);
    ey=Coord(Enod(:, i), 2);
    [sigLe, sigTe, tauLTe, zcoord] = laminatestress(ex, ey, h,...
        angles, Q, ed);
    sigL(:, i) = sigLe';
    sigT(:, i) = sigTe';
    tauLT(:, i) = tauLTe';
end
stressplot(Coord, Enod, sigL, sigT, tauLT)
