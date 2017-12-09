

% function [qd,vd,ad] = fcn()
qi1= 0;
qf1= 10;
qi2= 0;
qf2= 10;
qi3= 0;
qf3= 10;
vi1 = 0;
vi2 = 0;
vi3 = 0;
vf1 = 0;
vf2 = 0;
vf3 = 0;
d1 = [qi1,qf1,vi1,vf1,0,5];
[qd1,vd1,ad1] = cubic(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6));
d2 = [qi2,qf2,vi2,vf2,0,5];
[qd2,vd2,ad2] = cubic(d2(1),d2(2),d2(3),d2(4),d2(5),d2(6));
d3 = [qi3,qf3,vi3,vf3,0,5];
[qd3,vd3,ad3] = cubic(d3(1),d3(2),d3(3),d3(4),d3(5),d3(6));
qd = [qd1;qd2;qd3]'
vd = [vd1;vd2;vd3]' 
ad = [ad1;ad2;ad3]'
function [qd,vd,ad] =cubic(qi, qf, vi,vf,ti,tf)

t = tf-ti;
c = ones(size(t)); 
A=[1,ti, ti^2,ti^3;
    0,1,2*ti,3*ti^2;
    1,tf, tf^2,tf^3;
    0,1,2*tf,3*tf^2];

B=[qi;vi;qf;vf];
a=A\B;
x= 1;
qd = a(1).*c + a(2).*x +a(3).*x.^2 + a(4).*x.^3 ;
vd = a(2).*c +2*a(3).*x +3*a(4).*x.^2 ; 
ad = 2*a(3).*c + 6*a(4).*x;
% qd = a(1).*c + a(2).*t +a(3).*t.^2 + a(4).*t.^3 ;
% vd = a(2).*c +2*a(3).*t +3*a(4).*t.^2 ; 
% ad = 2*a(3).*c + 6*a(4).*t;
end


