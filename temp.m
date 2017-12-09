

function u = temp(qd,vd,ad,vd_measured)
Kv = diag([1 1 1]);
Kp = diag([1 1 1]);
%       0 1 0;
%       0 0 1;];
% Kp =  [1 0 0;
%        0 1 0;
%        0 0 1;];

aq = ad - Kv*(vd) - Kp*(qd)
 [M,C,G]=dynamics(qd);
%  display(vd_measured)
 u = M*aq + C*vd_measured + G
 
function [M,C,G] = dynamics(qd)
I = [0.006757, 0.0006036, 0.0015514;
    0.001745, 0.0005596, 0.00006455;
    0.00706657, 0.0006254, 0.0015708];
m = [1.01992, 0.3519, 0.22772]; % masses of the link in kg
l = [0.25107, 0.191, 0.37843];% lengths of links in m
r = [0.10424, 0.14550, 0.203];

t1= qd(1);
t2= qd(2);
t3= qd(3);

M11 = I(5)*(sin(t2)^2) + I(6)*(cos(t2+t3)^2)+ m(2)*I(7)*(cos(t2)^2)*(r(2)^2)+ m(3)*(l(1)*cos(t2)+ r(2)*cos(t2+t3))^2;
M12 = 0;
M13 = 0;
M21 = 0;

M22 = I(2) + I(3) + m(2)*(l(1))^2 + m(2)*(r(1))^2 + m(3)*(r(3))^2 +2 * m(3)*l(1)*r(2)* cos(t3);


M23 = I(3) + m(3)*(r(2))^2 + m(3)*l(1)*r(1)*cos(t3);


M31 = 0;
M32 = I(3) + m(3)*(r(1))^2 + m(3)*l(1)*r(2)*cos(t3);
M33 = I(3) + m(3)* (r(2))^2;


M =  [ M11 M12 M13;
       M21 M22 M23;
       M31 M32 M33];
%% Gravity is a bitch
g = 9.81;

G1 = 0;


G2 = -(m(2)*r(1) + m(3)*l(1))*g* cos(t2)- m(3)*r(2)*cos(t2 + t3);

G3 = -m(3)*g*r(2)*cos(t2 + t3);

G = [G1;G2;G3];
%% coriolis matrix
t23 = t2 + t3;
C = zeros(3);

gamma_001 = 0.5*( I(8) - I(9) - m(2)*r(2)^2)*sin(2*t2) + 0.5*( I(5)-I(9) )*sin(2*t23)-m(3)*(l(2)*cos(t2) + r(3)*cos(t23))*(l(2)*sin(t2) + r(3)*sin(t23));

gamma_002 = 0.5*(I(6) - I(9))*cos(2*t23) -m(3)*r(3)*sin(t23)*(l(2)*cos(t2) + r(3)*cos(t23));

gamma_010 = 0.5*(I(5) - I(8) - m(2)*r(2)^2)*cos(2*t2) + 0.5*( I(6) - I(9))*cos(2*t2)-m(3)*(l(2)*cos(t2) + r(3)*cos(t23))*(l(2)*sin(t2) + r(3)*sin(t23));


gamma_020 = 0.5*(I(6) - I(9))*sin(2*t23) -m(3)*r(3)*sin(t23)*(l(2)*cos(t2) + r(3)*cos(t23));

gamma_100 = 0.5*(I(8) - I(5) + m(2)*r(2)^2)*sin(2*t2) + 0.5*(I(9) - I(6))*sin(2*t23) + m(3)*( l(2)*cos(t2) + r(3)*cos(t23))*(l(2)*sin(t2) + r(3)*sin(t23));

gamma_112 = -l(2)*m(3)*r(2)*sin(t3);

gamma_121 = gamma_112;
gamma_122 = gamma_112;


gamma_200 = 0.5*(I(9) - I(6)) + m(2)*r(3)*sin(t23)*(l(2)*cos(t2) + r(3)*cos(t23));
gamma_211 = l(2)*m(3)*r(3)*sin(t3);


C(1) = gamma_001 + gamma_002;
C(4) = gamma_010;
C(7) = gamma_020;
C(2) = gamma_100;
C(5) = gamma_112;
C(8) = gamma_121 + gamma_122;

C(3) = gamma_200;
C(6) = gamma_211;


% qdd =  inv(M)*(u - C - G);
end    
end