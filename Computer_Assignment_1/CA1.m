clc
clear all
close all
%% CA1
% Created by Oscar Holke and Chirag Gurmurthy during parts of January 
% and February of 2019 as part of the course TME240 Composite Mechanics 
% at Chalmers University of Technology.

%% Input
% === Design ===
theta   = [45 90 0 0 90 45]; % lam. orientation (write out all)
sym     = 1; % sym=1 if symetric laminate, sym=2 if not symetric laminate
V_f     = 0.6; % Volume fraction of fibres
V_m     = 0.4; % Volume fraction of matrix

% = Geometry =
th      = 1.2; % define hight of laminate [unit=mm]

% = Fibres =
E_f     = 350e+03; % Young's modulus for fibre [unit=MPa]
nu_f    = 0.2; % Poissons ratio for fibre
rho_f   = 2e-03; % Density of fibre [unit=g/mm^3]
alpha_f = -1e-06;

% = Matrix =
E_m     = 3.5e+03; % Young's modulus for matrix [unit=MPa]
nu_m    = 0.35; % Poissons ratio for matrix
rho_m   = 1.2e-03; % Density of matrix [unit=g/mm^3]
alpha_m = 50e-06;

% = Halpin-Tsai =
xsi_E   = 2; % Choose the value for xsi in Halpin-Tsai method for 
             % Young's modulus
xsi_G   = 1; % Choose the value for xsi in Halpin-Tsai method for 
             % shear modulus

% === Load case ===
T             = [125 25]; % temperature [T(1)=starting temp 
                          % and T(2)=end temp]
eps_0_given   = [0.0005 0 0]; % initial added strain
N             = [0 0 0]'; %% Added mechanical force [unit=N/m]
M             = [-75 -50 0]'; % Added mechanical moment [unit=Nm/m]

% == set load cases ==
L_case  = 1; % 1=thermal only, 2=combined load cases



%% Setup
[Q, alpha, E_T] = laminadata( E_f, E_m, nu_f, nu_m, V_f, V_m,...
    alpha_f, alpha_m, xsi_E, xsi_G );

% = create h vector =
h=zeros(1,length(theta)+1);
m=1;
for n=th/2:-th/length(theta):-th/2
    h(m)=n;
    m=m+1;
end

% = setup zero matrices =
A          =   zeros(size(Q));
B          =   A;
D          =   A;
N          =   zeros(size(A,1),1);

%% Solve strain and stress problems
for i=1:length(theta)
    % = create Q_bar =
    [ T1 , T2]  =   CMTd(theta(i));
    Q_bar       =   T1\Q*T2;
    % = create A, B, D matrices =
    A           =   A + Q_bar*(h(i)-h(i+1));
    if sym==2
        B       =   B + 1/2*Q_bar*(h(i)^2-h(i+1)^2);
    end
    D           =   D + 1/3*Q_bar*(h(i)^3-h(i+1)^3);
end

% = compute N from given strain = (not general for all cases)
if L_case==1
    N           =  [0 0 0]';
    M           =  [0 0 0]';
elseif L_case==2
    ep0x       =   eps_0_given(1);
    ep_yxy     =   A([2 3],[2 3])\(N([2 3])-A([2 3],1)*ep0x);
    N(1)       =   A(1,1)*ep0x+A(1,[2 3])*ep_yxy;
end

% = compute strain =
ep0         =   (A - B*inv(D)*B)\(N - B*inv(D)*M);
k           =   D\(M - ep0);

p=0;    
for i=1:length(theta)
    [ T1 , T2]  =   CMTd(theta(i));
    Q_bar       =   T1\Q*T2;
    alpha_g     =   T2\alpha;
    for o=0:1
        p=p+1;
        z(p)=h(i+o);
        ep          =   ep0 + z(p)*k;
        ep_M(:,p)   =   ep-alpha_g*(T(2)-T(1));
    
    % = compute stress =
        sigma(:,p)  =   Q_bar*ep_M(:,p);
    end
end

p=0;
for i=1:length(theta)
    [T1] = CMTd(theta(i));
    for j=1:2
        p=p+1;
        sigma_LT(:,p)=T1*sigma(:,p);
    end
end

%% Plotting
% == plot sigma_global ==
labely = 'z';
labelx_g = {'sigma_x','sigma_y','tau_{xy}'};
legend_plot = {'Stressdistribution inn global coordinates',...
    'Stressdistribution inn lokal coordinates'};

figure(1)
for i=1:size(sigma,1)
    y=-z;
    x=sigma(i,:);
    ax(i)=subplot(3,1,i);
    plot(ax(i),x,y);
    if i==1
    title(legend_plot{1});
    end
    xlabel(labelx_g{i})
    ylabel(labely)
end

% axis([ax(1) ax(2) ax(3)],[-1700 1700 -0.6 0.6])

% == plot sigma_TL ==
labelx_lok = {'sigma_L','sigma_T','tau_{LT}'};

figure(2)
for i=1:size(sigma,1)
    y=-z;
    x=sigma_LT(i,:);
    ax(i)=subplot(3,1,i);
    plot(ax(i),x,y);
    if i==1
    title(legend_plot{2});
    end
    xlabel(labelx_lok{i})
    ylabel(labely)
end

% axis([ax(1) ax(2) ax(3)],[-1700 1700 -0.6 0.6])









