function [A, B, D, A_tilde] = foo(E_L, E_T, nu_LT,nu_TL, G_LT, G_TT,...
    angles, thickness)
% midplane
% = create h vector =
h=zeros(1,length(angles)+1);
m=1;
for n=(thickness/2):(-thickness/length(angles)):(-thickness/2)
    h(m)=n;
    m=m+1;
end   
p=0; 
for i=1:length(angles)
   for o=0:1
        p=p+1;
        z(p)=h(i+o);
   end
end
% extensional stiffness matrix
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
A_tilde=zeros(2,2);
for i=1:length(angles)
[T1,T2,Q,T_1tilde,Q_tilde] = stiffness(angles(i),E_L,E_T,G_LT,G_TT,...
    nu_LT,nu_TL);

Q_b = T1\Q*T2; %global stiffness matrix
A= A+Q_b*(h(i)-h(i+1)); %extensional stiffness matrix
B=B+0.5*(Q_b*(h(i)^2-h(i+1)^2)); %coupling stiffness matrix
D=D+0.3334*(Q_b*(h(i)^3-h(i+1)^3)); %bending stiffness matrix

A_tilde= A_tilde + ((1)*[Q_tilde(2,2) Q_tilde(2,1);
                             Q_tilde(1,2) Q_tilde(1,1)]*(h(i)-h(i+1)));           
end

end




