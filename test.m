%  A prorgam to find the torques of three joints for movement of
%  end-effector from source location to goal location
%  PID controller is used to control the torqe calculation and movement
%  cubic polynomial velocity principle is applied for trajectory Planning
%  By : Pranaya Ranjan Behera, Electronics and Instrumentation Engg.
%  
%--------------------------------------------------------------------------

clear all;
clc
g=9.8 ;          % in mitre/sec^2
m2= 12;             % m2 and m3 are in Kg; mass m1=0;
m3= 8;
l= .68;           % length in mitre  length1=0; length2=length3=l;
kp= 1;
kv= 5;

ls=4*l*l; 
count=0;
while(count==0)
disp('Insert the source location of end-effector(in centi-miter)')
xs1= input('Xs=');
ys1= input('Ys=');
zs1= input('Zs=');  
xs=xs1/100; ys=ys1/100;  zs=zs1/100;  % converting into miter
% check singularity
ds2=(xs^2)+(ys^2)+(zs^2);
if (ls<ds2)
    disp('*** *** *** That is a Singular Point *** *** ***'); 
else
    count=1;
end   
end

while(count==1)
disp('Insert the destination(goal) location of end-effector(in centi-miter)')
xg1= input('Xg=');
yg1= input('Yg=');
zg1= input('Zg=');    
xg=xg1/100; yg=yg1/100;  zg=zg1/100; % converting into miter
% check singularity
dg2=(xg^2)+(yg^2)+(zg^2);
if (ls<dg2)
    disp('*** *** *** That is a Singular Point *** *** ***');
else
    count=0;
end
end
     
while(count==0)     
t= input('Insert time of movement(in seconds)=');
if (t<.999)    % time of travel should be minimun 1sec.
    disp('*** *** *** Time should be at least 1 second *** *** ***');
else
    count=1;
end
end


%  1 unit = 50ms, i.e.  each 50ms torques will be caculated

qs1=atan2(ys,xs)

qs3=acos(((((xs* cos(qs1))+ (ys* sin(qs1)))^2) + (zs^2)-(2*l*l) )/ (2*l*l))

qs2=asin(((zs*(1+cos(qs3))^2)-(((zs*(1+cos(qs3))^2)^2)- 2*(1+cos(qs3))*((zs^2)-(l^2)*(sin(qs3))^2))^0.5)/(2*l*(1+cos(qs3))))

%----------------------------------------------------------------

qg1=atan2(yg,xg)

qg3=acos(((((xg* cos(qg1))+ (yg* sin(qg1)))^2) + (zg^2)-(2*l*l) )/ (2*l*l))

qg2=asin(((zg*(1+cos(qg3))^2)-(((zg*(1+cos(qg3))^2)^2)- 2*(1+cos(qg3))*((zg^2)-(l^2)*(sin(qg3))^2))^0.5)/(2*l*(1+cos(qg3))))


%-----------------------------------------------------------------
q1=q_cubic_fn(qs1, qg1, t);
q2=q_cubic_fn(qs2, qg2, t);
q3=q_cubic_fn(qs3, qg3, t);

qd11=qd_cubic_fn(qs1, qg1, t);
qd22=qd_cubic_fn(qs2, qg2, t);
qd33=qd_cubic_fn(qs3, qg3, t);

qdd1=qdd_cubic_fn(qs1, qg1, t);
qdd2=qdd_cubic_fn(qs2, qg2, t);
qdd3=qdd_cubic_fn(qs3, qg3, t);



 for i=1:20*t
     qd1=qd11(i);
     qd2=qd22(i);
     qd3=qd33(i);
     qs1= qd1^2;
     qs2= qd2^2;
     qs3= qd3^2;           %  square of qd ;
     
     c1= cos(q1(i));
     s1= sin(q1(i));
     c2= cos(q2(i));
     s2= sin(q2(i));
     c3= cos(q3(i));
     s3= sin(q3(i));
     s223=sin(q3(i)+ q2(i)+q2(i));
     c12= cos(q1(i)+ q2(i));
     s12= sin(q1(i)+ q2(i));
     c23= cos(q3(i)+ q2(i));
     s23= sin(q3(i)+ q2(i));
     
 % Define 3 Error function for 3 joints
     
     
    
     gq1(i)= .5*m3*l*s1*s3*g + .5*m2*l*s1*s2*c2*g ;
     
     hq1(i)= .5*m2*l*l*((.5*s1*c2*s2*s2*qs1) + ((2/3)*c2*s2*qd1*qd2) + (.5*s1*c2*(c2*c2*qs1 + qs2)) + (c2*s2*qd1*qd2*(c1*c2-(5/3)))) + .5*m2*l*l*((s1*s3*(s23+c2- .5*c23)*qs1) + ( 2*c1*c3*s2 - (2/3)*s223 + 4*s2*c2 )*qd1*qd2 + (2/3)*s223*qd1*qd3 );
     
%--------------------------------------------------------------      
    
     
     gq2(i)= .5*m3*l*g*(c1- 2*c23 - 2*c2) + .5*m2*l*g*(c1- 2*c2) ;
     
     hq2(i)= .5*m3*l*l*((2*s3+c1*s2)*qs2 + 3*(s3- .5*c1*s23)*(qd1+qd2)^2  - ((5/3)*c23 - 2*s23 -2*c2)*s23*qs1 - 2*s2*(c23 - s23 -c2)*qs1 ) + .5*m2*l*l*((1/3)*c2*s2*qs1 + .5*c2*s2*qs2);
     
%----------------------------------------------------------------

     
     
     gq3(i)= m3*l*g*(c1 - 2*c23)/2;
     
     hq3(i)=m3*l*l*((2*s3 - c1*s2)*qs2 - 1.5*c1*s23*(qd2+qd3)^2 - (s23*((4/3)*c23 - 2*s23 - 2*c2)));
     
%----------------------------------------------------------------  
 
     er1(i)=0;                    % Define error Functions: (desired values - values from sensors)
     er2(i)=0;
     er3(i)=0;  
     erd1(i)=0;
     erd2(i)=0;
     erd3(i)=0;
     
     
     k1=(qdd1(i)+(kp*er1(i))+(kv*erd1(i)));
     k2=(qdd2(i)+(kp*er2(i))+(kv*erd2(i)));
     k3=(qdd3(i)+(kp*er3(i))+(kv*erd3(i)));
%----------------------------------------------------------------
    
     mq1(i)= (((m2*l*l*(((s2/3)-(c1*c2*s2)/4)*s2 + (.5 - .25*c1*c2)*(c2^2))) + ( .5*m3*l*l*( (c23 - 2*c2 - c1*c3)*c2 + .5*c1*c3*c23 )))*k1) + ((.25*m2*l*l*s1*s2 - .5*m3*l*l*((2/3)*s2*s23 + s1*(s2-c2)*(c3+ .5*s1 )))*k2) + ((.25*m3*l*l*s1*(s2-2*c2))*k3);
          
     
     mq2(i)= ((((m3*l*l*(1- c3 - c1*c2 - .5*c1*c23 ))+( m2*l*l*(.5 - .5*c1*c2)))/2)*k2) + (( m3*l*l*(1-.5*c1*c23 -c3))*k3);
          
     
     mq3(i)= (( (4*c3 - 2*c1*c23 - c1*c23)*m3*l*l/4)*k2) + ((-c1*c23*m3*l*l/4)*k3);
        
%----------------------------------------------------------------
    
     tau1(i)= ( mq1(i)) + gq1(i) + hq1(i) ;
     
     tau2(i)= ( mq2(i)) + gq2(i) + hq2(i) ;
     
     tau3(i)= ( mq3(i)) + gq3(i) + hq3(i) ;
  
     
 end
 tau1;
 tau2;
 tau3;
 
 nn=size(tau1);
 n=nn(2);
 j=1:1:n;
 stem(tau1, j/20);  % As 1 unit= 50ms i.e. 1sec/20 = 50ms
 
 subplot(4,1,1)
 stem(j/20,tau2 );
 subplot(4,1,2)
 stem(j/20,tau2 );
 subplot(4,1,3)
 stem(j/20, tau3);
 %  A function to find the angular velocity from angular displacement by
%  cubic polynimial principle
% A trajectory Planning

% q=cumsum(qd/20)
% 
% i=size(qd);
% j=1:1:i(2);
% stem(j/20,qd);
%  A function to find the angular velocity from angular displacement by
%  cubic polynimial principle
% A trajectory Planning


function qdd = qdd_cubic_fn(qs, qg, t)

a=qs;
%b=0;
c=3*(qg-qs)/(t*t);
d=-2*(qg-qs)/(t*t*t);

% 1 unit = 50ms
% -----------------------------------------------------------------------------------------------------------------

tt= ceil(20*t);
for i=1:tt
        qdd(i)=((2*c)+(6*d*i/20));
end
qdd;

%for i=1:tt
 %       qd(i)=((2*c*i/20)+(3*d*i*i/400));
%end
%qd
%for i=1:tt
 %       q(i)=(a+(c*i*i/400)+(d*i*i*i/8000));
%end
%q
%  A function to find the angular velocity from angular displacement by
%  cubic polynimial principle
% A trajectory Planning
end

function qd = qd_cubic_fn(qs, qg, t)

a=qs;
%b=0;
c=3*(qg-qs)/(t*t);
d=-2*(qg-qs)/(t*t*t);

% 1 unit = 50ms
% -----------------------------------------------------------------------------------------------------------------

tt= ceil(20*t);
for i=1:tt
        qd(i)=((2*c*i/20)+(3*d*i*i/400));
end
qd;

%for i=1:tt
%        qdd(i)=((2*c)+(6*d*i/20));
%end
%qdd


%qd
%for i=1:tt
 %       q(i)=(a+(c*i*i/400)+(d*i*i*i/8000));
%end
%q
%  A function to find the angular velocity from angular displacement by
%  trapezoidal angular velocity principle
% A trajectory Planning with trapezoidal velovity for rotion more than pi/4
% and tringular velocity profile less than equal to pi/4 

end
function qd = qd_trap_fn(q)


   
 % 1 unit = 50ms
 % -----------------------------------------------------------------------------------------------------------------
 if(q>0)
if (q>pi/4)
    ti= (4*q)/pi;
 ti20=floor(20*ti);
 ti21=floor(20*(ti+1));
 disp('case 1')
    
    for i=1:20
        qd(i)=(pi/4)*i/20;
    end
    
    for i=21:ti20
        qd(i)=pi/4;
    end
    for i=(ti20)+1:(ti21)
    qd(i)=(pi/4)*(ti21-i)/20;
    end
    
else (q<=pi/4)
    ti=abs(sqrt(4*q/pi));
    ti20=floor(20*ti);
    disp('case 2')
     for i=1:ti20
        qd(i)=(pi/4)*(i/20);
     end
     for i=(ti20)+1:2*ti20
        qd(i)=(pi/4)*(2*ti20-i)/20;
     end
end
else
   if (q>=-pi/4 )
      q=abs(q);
      disp('case 3')
         ti20=floor(20*abs(sqrt(4*q/pi)));
     for i=1:ti20
        qd(i)=(-pi/4)*(i/20);
     end
     for i=(ti20)+1:2*ti20
       qd(i)=(-pi/4)*(2*ti20-i)/20;
     end
     
   else 
    q=abs(q);
    ti= (4*q)/pi;
    ti20=floor(20*ti);
    ti21=floor(20*(ti+1));
 
 disp('case 4')
    for i=1:20
        qd(i)=(-pi/4)*(i/20);
    end

    for i=21:ti20
        qd(i)=-pi/4;
    end
    for i=(ti20+1):(ti21)
    qd(i)=(-pi/4)*(ti21-i)/20;
    end
   end
    
 end
end


function q = q_cubic_fn(qs, qg, t)

a=qs;
%b=0;
c=3*(qg-qs)/(t*t);
d=-2*(qg-qs)/(t^3);

% 1 unit = 50ms
% -----------------------------------------------------------------------------------------------------------------

tt= ceil(20*t);

for i=1:tt
        q(i)=(a+(c*i*i/400)+(d*i*i*i/8000));
end

%for i=1:tt
%        qdd(i)=((2*c)+(6*d*i/20));
%end
%qdd

%for i=1:tt
 %       qd(i)=((2*c*i/20)+(3*d*i*i/400));
%end
%qd
end