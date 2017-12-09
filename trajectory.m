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
t= linspace(0,5,100*5);
subplot(3,1,1);
plot(t,qd1);
xlabel('Time');
ylabel('Angle(deg)for Link 1');
title('Position Trajectory of Link1 ');
subplot(3,1,2)
plot(t,vd1);
xlabel('Time');
ylabel('Velocity for Link 1');
title('Velocity Trajectory of Link1');
subplot(3,1,3)
plot(t,ad1);
xlabel('Time');
ylabel('Acceleration for Link 1');
title(' Acceleration Trajectory of Link1');
figure
subplot(3,1,1);
plot(t,qd2);
xlabel('Time');
ylabel('Angle for Link 2');
title('Position Trajectory of Link2 ');
subplot(3,1,2);
plot(t,vd2);
xlabel('Time');
ylabel('Velocity for Link 2');
title('Velocity Trajectory of Link2');
subplot(3,1,3);
plot(t,ad2);
xlabel('Time');
ylabel('Acceleration for Link 2');
title(' Acceleration Trajectory of Link2');
figure;
subplot(3,1,1);
plot(t,qd3);
xlabel('Time');
ylabel('Angle for Link 3');
title('Position Trajectory of Link3 ');
subplot(3,1,2);
plot(t,vd3);
xlabel('Time');
ylabel('Velocity for Link 3');
title('Velocity Trajectory of Link3');
subplot(3,1,3);
plot(t,ad3);
xlabel('Time');
ylabel('Acceleration for Link 3');
title(' Acceleration Trajectory of Link3');

[x1,x2,y1,y2,x3,y3] =Forward_Kinematics(qd1,qd2,qd3);
Link1=[0 x1 y1];
Link2=[0 x2 y2];
Link3=[0 x3 y3];
figure
for i = 1:10:500
plot([x1(i),x2(i),x3(i)],[y1(i),y2(i),y3(i)],'-o','LineWidth',1,'MarkerSize',5,'MarkerFaceColor',[1 1 1]);
xlabel('x');
ylabel('y');
axis([-1000 1000 -1000 1000])
pause (0.01);
end

function [x1,y1,x2,y2,x3,y3] = Forward_Kinematics(t1,t2,t3)
l1 = 0.25107; 
l2 = 0.191;
l3 = 0.37843;
x1 = l1*cosd(t1);
y1 = l1*sind(t1);
x2 = l1*cosd(t1) + l2*cosd(t1+t2);
y2 = l1*sind(t1) + l2*sind(t1+t2);
x3 = l1*cosd(t1) + l2*cosd(t1+t2)+l3*cosd(t1+t2+t3);
y3 = l1*sind(t1) + l2*sind(t1+t2)+l3*cosd(t1+t2+t3);
end



function [qd,vd,ad] =cubic(qi, qf, vi,vf,ti,tf)

t = linspace(ti,tf,100*(tf-ti));
c = ones(size(t)); 
A=[1,ti, ti^2,ti^3;
    0,1,2*ti,3*ti^2;
    1,tf, tf^2,tf^3;
    0,1,2*tf,3*tf^2];

B=[qi;vi;qf;vf];
a=A\B;

qd = a(1).*c + a(2).*t +a(3).*t.^2 + a(4).*t.^3 ;
vd = a(2).*c +2*a(3).*t +3*a(4).*t.^2 ; 
ad = 2*a(3).*c + 6*a(4).*t;
end

