syms theta1 theta2 theta3 L1 L2 L3
DH =[[0,  pi/2, L1, theta1]
     [L2, 0,    0,  theta2] 
     [L3,-pi/2, 0,  theta3]] 
 %% Frame transformations
 T1 = FrameTransform(DH(1,:));
 T2 = FrameTransform(DH(2,:));
 T3 = FrameTransform(DH(3,:));
 
 T01 = T1;
 T02 = simplify(T01*T2);
 T03 = simplify(T02*T3);
 
J1= simplify(jacobian(T01(1:3,4),[theta1]));
J2= simplify(jacobian(T02(1:3,4),[theta1 theta2]));
J3= simplify(jacobian(T03(1:3,4),[theta1 theta2 theta3]));

syms dq1 dq2 dq3 
v_m1 = simplify(J1 * dq1);
v_m2 = simplify(J2 * [dq1;dq2]);
v_m3 = simplify(J3 * [dq1;dq2;dq3]);

syms m1 m2 m3 I1 I2 I3
k1 =simplify(( 0.5 * m1 * (v_m1.' * v_m1)) + (0.5 * I1 * (dq1)^2));
k2 =simplify(( 0.5 * m2 * (v_m2.' * v_m2)) + (0.5 * I2 * (dq1+dq2)^2));
k3 =simplify(( 0.5 * m3 * (v_m3.' * v_m3)) + (0.5 * I3 * (dq1+dq2+dq3)^2));

K = k1 +k2+ k3 ;

g = 9.8;

P1 = simplify(m1 * g * T01(3,4));
P2 = simplify(m2 * g * T02(3,4));
P3 = simplify(m3 * g * T03(3,4));

P = P1 + P2 + P3 ;

L = simplify ( K - P );
syms th1(t) th2(t) th3(t)
syms ddq1 ddq2 ddq3

dL_dq1 = diff(L,dq1);
dL_dq1_sub= subs(dL_dq1,[theta1 theta2 theta3 dq1 dq2 dq3],[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t)]);
dL1 = diff(dL_dq1_sub,t);
dL1_a = subs(dL1,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t)...
         diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t),t,t)],[theta1 theta2 theta3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
     
dL_dq2 = diff(L,dq2);
dL_dq2_sub= subs(dL_dq2,[theta1 theta2 theta3 dq1 dq2 dq3],[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t)]);
dL2 = diff(dL_dq2_sub,t);
dL2_a = subs(dL2,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t)...
         diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t),t,t)],[theta1 theta2 theta3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
     
dL_dq3 = diff(L,dq3);
dL_dq3_sub= subs(dL_dq3,[theta1 theta2 theta3 dq1 dq2 dq3],[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t)]);
dL3 = diff(dL_dq3_sub,t);
dL3_a = subs(dL1,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t)...
         diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t),t,t)],[theta1 theta2 theta3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

dL1_b = diff(L,theta1);
dL2_b = diff(L,theta2);
dL3_b = diff(L,theta3);

Tau1 = dL1_a - dL1_b;
Tau2 = dL2_a - dL2_b;
Tau3 = dL3_a - dL3_b;

M11 = simplify(Tau1 - subs(Tau1,ddq1,0))/ddq1;
M12 = simplify(Tau1 - subs(Tau1,ddq2,0))/ddq2;
M13 = simplify(Tau1 - subs(Tau1,ddq3,0))/ddq3;
M21 = simplify(Tau2 - subs(Tau2,ddq1,0))/ddq1;
M22 = simplify(Tau2 - subs(Tau2,ddq2,0))/ddq2;
M23 = simplify(Tau2 - subs(Tau2,ddq3,0))/ddq3;
M31 = simplify(Tau3 - subs(Tau3,ddq1,0))/ddq1;
M32 = simplify(Tau3 - subs(Tau3,ddq2,0))/ddq2;
M33 = simplify(Tau3 - subs(Tau3,ddq3,0))/ddq3;

M = simplify([M11 M12 M13;
     M21 M22 M23;
     M31 M32 M33])
 
G1 = subs(Tau1, [dq1 dq2 dq3 ddq1 ddq2 ddq3],[0 0 0 0 0 0]);
G2 = subs(Tau2, [dq1 dq2 dq3 ddq1 ddq2 ddq3],[0 0 0 0 0 0]);
G3 = subs(Tau3, [dq1 dq2 dq3 ddq1 ddq2 ddq3],[0 0 0 0 0 0]);

G = simplify([G1;G2;G3])

C1 = Tau1- (M(1,:)* [ddq1 ddq2 ddq3].' + G1);
C2 = Tau2- (M(2,:)* [ddq1 ddq2 ddq3].' + G2);
C3 = Tau3- (M(3,:)* [ddq1 ddq2 ddq3].' + G3);

C = simplify([C1;C2;C3])
function [ x ] = FrameTransform( y )
a = y(1);
alpha = y(2) * (pi/180);
d = y(3);
theta = y(4) * (pi/180);
RotzTheta = [cos(theta), -sin(theta), 0, 0;
             sin(theta), cos(theta), 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1];
TranszD = [1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, d;
            0, 0, 0, 1];
TransxA = [1,0,0,a;
           0,1,0,0;
           0,0,1,0;
           0,0,0,1];
RotxAlpha = [1,0,0,0;
             0, cos(alpha), -sin(alpha), 0;
             0 sin(alpha), cos(alpha), 0;
             0,0, 0, 1];
x = RotzTheta * TranszD * TransxA * RotxAlpha;
end
