clear all
close all
clc


% Add casadi to path - dependent on personal machine
addpath('C:\MatlabHelpers\Casadi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

% T is the time horizon (This is MPC code so time horizon of 1 sec but
% total time is 10 secs. As shown in the loop later)
T = 1;
% N is the number of shooting segments per horizon. So 10 segments per 1
% sec.
N = 10;
%% Model Properties
L = 1; % Length of the pendulum
m = 1; % Mass of the pendulum
g = 9.81; % Gravity
%% Casadi Stuff
% Create the Casadi variables (MX vs SX means something but not entirely
% sure what. I think MX is more general, so used here)
qA = MX.sym('qA');
qB = MX.sym('qB');
qA_dot = MX.sym('qA_dot');
qB_dot = MX.sym('qB_dot');
% Combine the Casadi variables into a state vector
x = [qA;qB;qA_dot;qB_dot];
% Casadi variables for the input
TA = MX.sym('TA');
TB = MX.sym('TB');
% The TB=0 line is for an underactuated case
%TB = 0;
% Combine into an input vector
u = [TA;TB];
% in the underactuated case u=TA;
%u = [TA];

% Define the derivative of the state (in this case both second derivatives
% of position)
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

% Create a vector called 'ode' equivalent to xdot
ode = [qA_dot;qB_dot;qA_ddot;qB_ddot]; 
% makes a casadi function that takes x and u as inputs and outputs ode
% (xdot)
% Functions are of the form:
%   Function('name',{inputs},{outputs},{input names},{output names}) 
f = Function('f',{x,u},{ode},{'x','u'},{'ode'});

% Create a variable for intg_options
% in C++ this is a dict variable
intg_options = struct;
% Set the final time of the integrator to be T/N --- time for 1 shooting
% segment
intg_options.tf = T/N; 
intg_options.simplify = true;
intg_options.number_of_finite_elements = 4;

% Another variable for a differential algebraic equation
% Another area for potential changes, lots of ways to define this I think.
% See below where we make an integrator with dae and intg_options as inputs
dae = struct;
dae.x = x;
dae.p = u; % parameters of the dae are our control inputs
dae.ode = f(x,u); % The ode of our dae is just our derivative function from above
% creates an integrator (think function) for our dae
intg = integrator('intg','rk',dae,intg_options);

% res is the "result" of the integrator as a structure
res = intg('x0',x,'p',u);
x_next = res.xf; %This gets us our next state (symbolic still I think)

% Creates a function that takes x and u as input and x_next as output
% x_next is at end of 1 shooting segment
F = Function('F',{x,u},{x_next},{'x','u'},{'x_next'});

% This commented out line will create a simulation
% Lines uAct and res will add numbers to the simulation
%sim = F.mapaccum(N);

x0 = [0;0;0;0];
%uAct = [ones(1,N)*3*9.81;ones(1,N)*9.81];
%res = sim(x0,uAct);
% 
% figure
% Creates a time vector
tgrid = linspace(0,T,N+1);
dt = tgrid(2)-tgrid(1);
%%% These commented out lines would plot the simulation from above
% plot(tgrid,full([x0 res]));
% legend('qA','qB','qA_dot','qB_dot')
% xlabel('t (s)')


% I think the opti class respresents an optimization problem
opti = casadi.Opti();
% Add a 4xN+1 matrix (think N+1 state vectors concatenated
x = opti.variable(4,N+1);
% Same with control, but 1 less control input than state
u = opti.variable(2,N);
% Assume cost function is J=x'Qx + u'Ru
% Add the Q matrix as a parameter (this might not be necessary)
Q = opti.parameter(4,4);
% Set the value of the Q matrix
opti.set_value(Q,[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1])
% Same with R
R = opti.parameter(2,2);
opti.set_value(R,[1,0;0,1])
% The initial state as a parameter
x0param = opti.parameter(4,1);
% initial time as a parameter
t0 = opti.parameter(1,1);
%initialGuess = MX(2,N);
%opti.set_initial(u, initialGuess);
%opti.set_value(x0param,[0;0;0;0]);

% initialize the cost as 0
cost = 0;
%x(:,1) = x0param;
% This loops through each of the shooting segments
for k = 1:N
    % for MPC, t0 is different in every optimization problem.
    % This matter for a tracking problem where tracking a function of time
    % Get the current time for the current shooting segment independent of
    % where we are in the MPC
    tCurrent = tgrid(k) + t0;
    % Get the cost --- (x-x_des)'*Q*(x-x_des) + u'Ru
    cost = cost + 1000*((x(:,k)-[pi/2*sin(tCurrent);pi/2*sin(tCurrent);pi/2*cos(tCurrent);pi/2*cos(tCurrent)])'...
        *Q*(x(:,k)-[pi/2*sin(tCurrent);pi/2*sin(tCurrent);pi/2*cos(tCurrent);pi/2*cos(tCurrent)]));% + u(:,k)'*R*u(:,k);
end
% add a final cost
cost = cost + ...
    1000*(x(:,N+1)-[pi/2*sin(tgrid(N+1)+t0);pi/2*sin(tgrid(N+1)+t0);pi/2*cos(tgrid(N+1)+t0);pi/2*cos(tgrid(N+1)+t0)])'...
    *Q*(x(:,N+1)-[pi/2*sin(tgrid(N+1)+t0);pi/2*sin(tgrid(N+1)+t0);pi/2*cos(tgrid(N+1)+t0);pi/2*cos(tgrid(N+1)+t0)]);
% Create a function that takes x and u as inputs and J as the output
J = Function('J',{x,u},{cost},{'x','u'},{'cost'});

%myCost = J(xTest,uTest);
%opti.minimize( 1000*sumsqr(xTest-[pi/2,0,0,0]') + sumsqr(uTest) );
% Give the optimization problem our cost function
opti.minimize(J(x,u));

% This loop adds the constraints on each shooting segment
% The next point has to be equal to the point you get when you integrate
% the dynamics from where you started. I think this is easier to draw out
% on paper
for k = 1:N
    opti.subject_to(x(:,k+1)==F(x(:,k),u(:,k)));
end
% The first state has to be equal to the x0 param.
opti.subject_to(x(:,1)==x0param);

%opti.solver('ipopt');
% Creates an option structure
opts = struct;
% Give it a solver (they used qrqp in the video for mpc), could be ipopt
opts.qpsol = 'qrqp';
% This changes the print settings
opts.print_header = false;
opts.print_iteration = false;
opts.print_time = false;
opts.qpsol_options.print_iter = false;
opts.qpsol_options.print_header = false;
opts.qpsol_options.print_info = false;
%opti.solver('sqpmethod',struct('qpsol','qrqp'));
%sol = opti.solve();
% Give the optimization problem the solver we want
opti.solver('sqpmethod',opts);
%%%% AT THIS POINT THE CODE IS JUST A SINGLE MULTIPLE SHOOTING CODE
% COULD CHANGE ABOVE CODE TO BE COLLOCATION OR ANY OTHER OPTIMIZATION
% APRROACH THEORETICALLY

%% MPC
% this creates a function M that takes x0 and t0 and gives the next control
% input? Not entirely clear to me right now
M = opti.to_function('M',{x0param,t0},{u(:,1)},{'x','t0'},{'u_opt'});
% initialize matrices as logs
x_log = [];
u_log = [];
% This sets the initial position I think
x0Test = [0,0,0,0]';
x = x0Test;
t0 = 0;
guess = [0;0];
tic
for i = 1:10*N
    tic
    % Gets the next control input (not sure why full)
     u = full(M(x,t0));
     toc
     % Add to the log
     u_log(:,i) = u;
     % Update the next initial guess for the optimal control
     guess = u;
     % Add to the log
     x_log(:,i) = x;
     % Next 'initial time'
     t0 = t0 + dt;
     % simulate to get the next state using the ideal control found (can
     % add noise here)
     x = full(F(x,u)) + [0;0;rand*.1;rand*.1];%+ [0;rand*0.1;0;rand*0.1];
end
toc

%% Plotting
figure
hold on
%plot(tgrid,sol.value(x(1:2,:)))
tVec = 0:dt:T*10;
plot(tVec(1:end-1),x_log(1:2,:))
%stairs(tgrid,[sol.value(uTest(1,:)) nan],'-.');
%stairs(tgrid,[sol.value(uTest(2,:)) nan],'-.');
xlabel('t (s)')
ylabel('Angle (rad)')
legend('qA','qB')%,'qA_dot','qB_dot','TA','TB')
%visualize(sol.value(x))
plot(tVec(1:end-1),pi/2*sin(tVec(1:end-1)))
legend('qA','qB','Desired')