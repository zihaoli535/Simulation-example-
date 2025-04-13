function dx=Example(t,x,xd)
%% System parameters (Inverted Pendulum System)
m1 = 2;     % kg (pendulum end mass)
m2 = 2.5;   % kg
J1 = 5;     % kg*m^2 (moment of inertia)
J2 = 6.25;  % kg*m^2
g = 9.81;   % m/s^2 (gravitational acceleration)
s = 100;    % N/m (spring constant)
b = 0.4;    % m (distance between pendulum hinges)
l = 0.5;    % m (spring natural length)
r = 0.5;    % m (height of pendulum)

%% Control parameters
Tf = 2;         % Settling time
varpi = 0.05; % Steady-state accuracy 
iota = 2;      % Positive integers (2ι > n+1)
c = 50;         % Transformation parameters
a = 0.5;       % Safety distance

%% Command filter parameters
wn = 80;  % Natural frequency
zeta = 0.85; % Damping ratio

%% Error compensation parameters
k11 = 5; k12 = 35;         % Subsystem 1 gains
k21 = 6; k22 = 35;         % Subsystem 2 gains

%% Adaptive law parameters and Control Coefficient
sigma1 = 0.05; sigma2 = 0.05; 
b1 = 1; b2 = 1; 
ze1=10;ze2=10;o1=0.005;o2=0.005;
g11=1;g12=1/5;g21=1;g22=1/6.25; %g12=1/J1 g22=1/J2

%% Reference signals
y1r=0.5*sin(t)+0.5*sin(0.5*t);
dy1r=0.5*cos(t)+0.25*cos(0.5*t);
y2r=sin(0.5*t)*sin(t);
dy2r=0.5*cos(0.5*t)*sin(t)+sin(0.5*t)*cos(t);

%% RBF NN parameters
num_nodes = 6;             % Number of NN nodes
centers = linspace(-2, 2, num_nodes); % Centers evenly spaced in [-2,2]
sigma = 2;                 % Width of Gaussian functions
X1=[x(1);x(2);x(5);x(7)];
X2=[x(3);x(4);x(6);x(8)];
phi1 = zeros(num_nodes,1);
for j = 1:num_nodes
    phi1(j) = exp(-norm(X1 - centers(j))^2 / sigma^2);
end
phi2 = zeros(num_nodes,1);
for j = 1:num_nodes
    phi2(j) = exp(-norm(X2 - centers(j))^2 / sigma^2);
end

%% Time-Varying Constraint Function
if t<=Tf
     B=(1/t-1/Tf)^(2*iota)+varpi;
     dB=2*iota*(1/t-1/Tf)^(2*iota-1)*(-1/(t^2));
else
    B=varpi;
    dB=0;
end

%% Subsystem 1 Error Dynamics
z11=x(1)-y1r;
r1=c*(B^2-z11^2);
if r1>0&&r1<=a
    h1=1-(r1/a-1)^(2*iota);
    p1=1/h1-4*c*iota*(r1/a-1)^(2*iota-1)*z11^2/(a*h1^2);
    q1=4*c*iota*(r1/a-1)^(2*iota-1)*B*dB*z11/(a*h1^2);
else
    h1=1;
    p1=1;
    q1=0;
end

%% Subsystem 2 Error Dynamics
z21=x(3)-y2r;
r2=c*(B^2-z21^2);
if (r2>0)&&(r2<=a)
    h2=1-(r2/a-1)^(2*iota);
    p2=1/h2-4*c*iota*(r2/a-1)^(2*iota-1)*z21^2/(a*h2^2);
    q2=4*c*iota*(r2/a-1)^(2*iota-1)*B*dB*z21/(a*h2^2);
else
    h2=1;
    p2=1;
    q2=0;
end

%% System Nonlinearities
f11=0; f12=((m1*g*r/J1-s*r^2/(4*J1))*sin(x(1))+s*r*(l-b)/(2*J1)); 
f21=0; f22=((m2*g*r/J2-s*r^2/(4*J2))*sin(x(3))+s*r*(l-b)/(2*J2));

%% Subsystem 1 Control Design
z11=x(1)-y1r;  s11=z11/h1;  XI11=s11-x(9);
alp11=-(k11*s11+1.5*p1*XI11+p1*x(5)*f11-p1*dy1r+q1+p1^2*x(9))/(p1*g11);

% Filter dynamics
dx13=wn*x(14);
dx14=-2*zeta*wn*x(14)-wn*(x(13)-alp11);

% Sliding surface and control law
s12=x(2)-x(13);  XI12=s12-x(10);
u1=-(k12*s12+1.5*XI12+x(5)*f12+g11*s11*p1-dx13+(XI12*x(7)*norm(phi1))/(XI12^2+sigma1))/g12;

%% Subsystem 2 Control Design
z21=x(3)-y2r;  s21=z21/h2;  XI21=s21-x(11);
alp21=-(k21*s21+1.5*p2*XI21+p2*x(6)*f21-p2*dy2r+q2+p2^2*x(11))/(p2*g21);

% Filter dynamics
dx15=wn*x(16);
dx16=-2*zeta*wn*x(16)-wn*(x(15)-alp21);

% Sliding surface and control law
s22=x(4)-x(15);  XI22=s22-x(12);
u2=-(k22*s22+1.5*XI22+x(6)*f22+g21*s21*p2-dx15+(XI22*x(8)*norm(phi2))/(XI22^2+sigma2))/g22;

%% System Dynamics Equations
dx1=x(2);% Subsystem 1 position
dx2=u1/J1+(m1*g*r/J1-s*r^2/(4*J1))*sin(x(1))+s*r*(l-b)/(2*J1)+s*r^2*sin(x(3))/(4*J1)+0.1*xd(1,1)*sin(xd(2,2))/J1+0.5*sin(t);

dx3=x(4); % Subsystem 2 position
dx4=u2/J2+(m2*g*r/J2-s*r^2/(4*J2))*sin(x(3))+s*r*(l-b)/(2*J2)+s*r^2*sin(x(1))/(4*J2)+0.1*xd(3,3)*cos(xd(4,4))/J1+0.2*cos(0.5*t);

%% Parameter Adaptation Laws
dx5=p1*f11*XI11+f12*XI12-o1*x(5);
dx6=p2*f21*XI21+f22*XI22-o2*x(6);
dx7=ze1*XI12^2*norm(phi1)/(XI12^2+sigma1)-b1*ze1*x(7);
dx8=ze2*XI22^2*norm(phi2)/(XI22^2+sigma2)-b2*ze2*x(8);

%% Auxiliary System Dynamics
dv11=-k11*x(9)+g11*p1*x(10)+g11*p1*(x(13)-alp11)-p1^2*x(9);
dv12=-k12*x(10)-g11*x(9);
dv21=-k21*x(11)+g21*p2*x(12)+g21*p2*(x(15)-alp21)-p2^2*x(11);
dv22=-k22*x(12)-g21*x(11);

%% State Derivatives Vector
dx=[dx1;dx2;dx3;dx4;dx5;dx6;dx7;dx8;dv11;dv12;dv21;dv22;dx13;dx14;dx15;dx16];