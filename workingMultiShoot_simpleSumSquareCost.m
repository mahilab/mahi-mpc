clear all
close all
clc


addpath('C:\MatlabHelpers\Casadi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
T = 10;
N = 200;
L = 1;
m = 1;
g = 9.81;
qA = MX.sym('qA');
qB = MX.sym('qB');
qA_dot = MX.sym('qA_dot');
qB_dot = MX.sym('qB_dot');

x = [qA;qB;qA_dot;qB_dot];
TA = MX.sym('TA');
TB = MX.sym('TB');
%TB = 0;
u = [TA;TB];
%u = [TA];
qA_ddot = -(TA - TB - TB*cos(qB) + L^2*m*qA_dot^2*sin(qB) + ...
    L^2*m*qB_dot^2*sin(qB) - 2*L*g*m*cos(qA) + ...
    L^2*m*qA_dot^2*cos(qB)*sin(qB) + 2*L^2*m*qA_dot*qB_dot*sin(qB) +  ...
    L*g*m*cos(qA + qB)*cos(qB))/(L^2*m*(cos(qB)^2 - 2));

qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + ...
    2*L*g*m*cos(qA + qB) + 3*L^2*m*qA_dot^2*sin(qB) + ...
    L^2*m*qB_dot^2*sin(qB) - 2*L*g*m*cos(qA) + ...
    2*L^2*m*qA_dot^2*cos(qB)*sin(qB) + L^2*m*qB_dot^2*cos(qB)*sin(qB) - ...
    2*L*g*m*cos(qA)*cos(qB) + 2*L^2*m*qA_dot*qB_dot*sin(qB) + ...
    L*g*m*cos(qA + qB)*cos(qB) + ...
    2*L^2*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(L^2*m*(cos(qB)^2 - 2));




ode = [qA_dot;qB_dot;qA_ddot;qB_ddot]; 

f = Function('f',{x,u},{ode},{'x','u'},{'ode'});

intg_options = struct;
intg_options.tf = T/N;
intg_options.simplify = true;
intg_options.number_of_finite_elements = 4;

dae = struct;
dae.x = x;
dae.p = u;
dae.ode = f(x,u);
intg = integrator('intg','rk',dae,intg_options);
res = intg('x0',x,'p',u);
x_next = res.xf;

F = Function('F',{x,u},{x_next},{'x','u'},{'x_next'});

sim = F.mapaccum(N);

x0 = [0;0;0;0];
%uAct = [ones(1,N)*3*9.81;ones(1,N)*9.81];
%res = sim(x0,uAct);
% 
% figure
tgrid = linspace(0,T,N+1);
% plot(tgrid,full([x0 res]));
% legend('qA','qB','qA_dot','qB_dot')
% xlabel('t (s)')



opti = casadi.Opti();
x = opti.variable(4,N+1);
u = opti.variable(2,N);
Q = opti.parameter(4,4);
opti.set_value(Q,[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1])
R = opti.parameter(2,2);
opti.set_value(R,[1,0;0,1])
x0param = opti.parameter(4,1);
opti.set_value(x0param,[0;0;0;0]);
cost = sumsqr(x) + sumsqr(u);
%x(:,1) = x0param;
% for k = 1:N
%     cost = cost + 1000*((x(:,k)-[pi/2*sin(tgrid(k));pi/2*sin(tgrid(k));pi/2*cos(tgrid(k));pi/2*cos(tgrid(k))])'...
%         *Q*(x(:,k)-[pi/2*sin(tgrid(k));pi/2*sin(tgrid(k));pi/2*cos(tgrid(k));pi/2*cos(tgrid(k))])) + u(:,k)'*R*u(:,k);
% end
% cost = cost + ...
%     1000*(x(:,N+1)-[pi/2*sin(tgrid(N+1));pi/2*sin(tgrid(N+1));pi/2*cos(tgrid(N+1));pi/2*cos(tgrid(N+1))])'...
%     *Q*(x(:,N+1)-[pi/2*sin(tgrid(N+1));pi/2*sin(tgrid(N+1));pi/2*cos(tgrid(N+1));pi/2*cos(tgrid(N+1))]);
J = Function('J',{x,u},{cost},{'x','u'},{'cost'});
resCost = J('x',x,'u',u);
%myCost = J(xTest,uTest);
%opti.minimize( 1000*sumsqr(xTest-[pi/2,0,0,0]') + sumsqr(uTest) );
opti.minimize(J(x,u));
for k = 1:N
    opti.subject_to(x(:,k+1)==F(x(:,k),u(:,k)));
end
opti.subject_to(x(:,1)==x0param);

opti.solver('ipopt');
%opti.solver('sqpmethod',struct('qpsol','qrqp'));

y = 7;
sol = opti.solve();

% M = opti.to_function('M',{x0param},{u(:,1)},{'x'},{'u_opt'});
% 
% x_log = [];
% u_log = [];
% 
% x0Test = [0,0,0,0]';
% x = x0Test;
% for i = 1:4*N
%     u = full(M(x));
%     
%     u_log(:,i) = u;
%     x_log(:,i) = x;
%     
%     x = full(F(x,u))+ [0;rand*0.02;0;rand*0.02];
% end
figure
hold on
plot(tgrid,sol.value(x(1:2,:)))
%stairs(tgrid,[sol.value(uTest(1,:)) nan],'-.');
%stairs(tgrid,[sol.value(uTest(2,:)) nan],'-.');
xlabel('t (s)')
legend('qA','qB')%,'qA_dot','qB_dot','TA','TB')
visualize(sol.value(x))
