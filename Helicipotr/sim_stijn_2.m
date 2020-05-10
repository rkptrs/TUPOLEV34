clear 
%INITIAL DATA HELICOPTER 
g=9.81;	 
cla=0.115 * 180 / pi; 
volh=0.092;	%blade solidity	 
lok=6; 
cds=3.0; 
mass=17650 * 0.45359237; 
rho=1.225; 
diam=2*24 * 0.3048; 
iy=86693.05; 
mast=1.2; 
omega_r = 221.2848; 
omega = omega_r/diam*2; 
vtip=omega*(diam/2); 
area=pi/4*diam^2; 
tau=.1;		%time constant in dynamiCs inflow!!! 
collect(1)=4.415*pi/180; 
longit(1)=0*pi/180; 
 
%initial values; 
u1 = 46.3; 
u2 = 36; 
u3 = 46.3; 
u4 = 56.9; 
u_list = [u2, u3, u4];
theta_wish_1 = -2.8722; 
theta_wish_2 = -1.7384; 
theta_wish_3 = -2.8722; 
theta_wish_4 = -4.2861; 
theta_list = [theta_wish_2, theta_wish_3, theta_wish_4];
collect_list = [5.139*pi/180, 5.799*pi/180, 6.769*pi/180];
c_des = 0;
h_des = 0;


t0=0; 
u0=u1; 
w0=0; 
q0=0; 
pitch0=theta_wish_1*pi/180; 
x0=0; 
labi0=sqrt(mass*g/(area*2*rho))/vtip; 
c0 = 0;

t(1)=t0; 
u(1)=u0; 
w(1)=w0; 
q(1)=q0; 
pitch(1)=pitch0; 
x(1)=x0; 
labi(1)=labi0; 
z(1)=0; 
delta_th(1) = 0; 
c(1) = c0;
corrc(1) = 0;
 
%INTEGRATION  
aantal=12000; 
teind=120; 
stap=(teind-t0)/aantal; 

% gains
K1 = 2;
K2 = 0.5;
K3 = 0.03;
K4 = 0.01;

%altitude hold
K1c = 0.4;
K2c = 0.0;
K3c = 0.03;

for j = 0:2

if j > 0
    t(j*aantal + 1) = t(j*aantal);
    u(j*aantal + 1)=u(j*aantal);
    w(j*aantal + 1)=w(j*aantal);
    q(j*aantal + 1)=q(j*aantal);
    pitch(j*aantal + 1)=pitch(j*aantal);
    x(j*aantal + 1)=x(j*aantal);
    labi(j*aantal + 1)=labi(j*aantal);
    z(j*aantal + 1)=z(j*aantal);
    %delta_th(j*aantal + 1) = delta_th(j*aantal);
end
    

    
for i=aantal*j+1:aantal*(j+1) 
    
theta_wish = -5*pi/180;
theta_wish = (theta_wish_1*pi/180 - pitch(i)) + 0.02*u(i);
longit(i)=K1*(pitch(i) - theta_list(j+1)*pi/180) + K2*q(i) + K3*(u_list(j+1) -u(i));% + K4*w(i);
% collect(i)= collect_list(j+1);
collect(i) = collect(1);
    



    
%NO LAW FOR COLLECTIVE 
 
c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i)); 
c_des(i) = 0;
h(i)=-z(i); 
% collect(i)=collect(1); 


c_des(i) = K3c * (h_des - h(i));
c_des(i) = 0;
longit(i)=longit(i) + K1c*(c_des(i) - c(i));% + K5*corrc(i);
collect(i) = collect_list(j+1) + K1c*(c_des(i) - c(i));


 
%Defining the differential equations 
 
%defining the nondimensional notations 
qdiml(i)=q(i)/omega; 
vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip; 
if u(i)==0 	if w(i)>0 	phi(i)=pi/2; 
        else phi(i)=-pi/2;end 
else 
phi(i)=atan(w(i)/u(i)); 
end 
if u(i)<0 
phi(i)=phi(i)+pi; 
end 
alfc(i)=longit(i)-phi(i); 
 
mu(i)=vdiml(i)*cos(alfc(i)); 
labc(i)=vdiml(i)*sin(alfc(i)); 
 
%a1 Flapping calculi 
teller(i)=-16/lok*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i)); 
a1(i)=teller(i)/(1-.5*mu(i)^2); 
 
%the thrust coefficient 
ctelem(i)=cla*volh/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i))); 
%Thrust coefficient from Glauert 
alfd(i)=alfc(i)-a1(i); 
ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*... 
sin(alfd(i))+labi(i))^2); 
 
%Equations of motion 
labidot(i)=ctelem(i);  
thrust(i)=labidot(i)*rho*vtip^2*area; 
helling(i)=longit(i)-a1(i); 
vv(i)=vdiml(i)*vtip; 		%it is 1/sqrt(u^2+w^2) 
 
udot(i)=-g*sin(pitch(i))-cds/mass*.5*rho*u(i)*vv(i)+... 
thrust(i)/mass*sin(helling(i))-q(i)*w(i); 
 
wdot(i)=g*cos(pitch(i))-cds/mass*.5*rho*w(i)*vv(i)-... 
thrust(i)/mass*cos(helling(i))+q(i)*u(i); 
 
qdot(i)=-thrust(i)*mast/iy*sin(helling(i)); 
 
pitchdot(i)=q(i); 
 
xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i)); 
 
zdot(i)=-c(i); 
 
labidot(i)=(ctelem(i)-ctglau(i))/tau; 
corrdot(i)=u1-u(i); 
corrcdot(i)=c_des(i)-c(i); 
%delta_th(i+1) = delta_th(i) + stap*deltathdot(i);

u(i+1)=u(i)+stap*udot(i); 
w(i+1)=w(i)+stap*wdot(i); 
q(i+1)=q(i)+stap*qdot(i); 
pitch(i+1)=pitch(i)+stap*pitchdot(i); 
x(i+1)=x(i)+stap*xdot(i); 
labi(i+1)=labi(i)+stap*labidot(i); 
z(i+1)=z(i)+stap*zdot(i); 
t(i+1)=t(i)+stap; 
corr(i+1) = corr(i) + stap*corrdot(i); 
corrc(i+1) = corrc(i) + stap*corrcdot(i);
end 
end

% plot(t, pitch * 180 / pi); hold on; plot(t(1:8000), h)
plot(t, u); hold on; plot(t, pitch *180 /pi); grid on; hold on; plot(t(1:36000), c);plot(t,z); 

