clc; clear all; close all;
syms qA qB qA_dot qB_dot TA TB

L = 1; m = 1; g = 9.81;

x = [qA; qB; qA_dot; qB_dot];
u = [TA; TB];

qA_ddot = -(TA - TB - TB*cos(qB) + L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));
qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + 2*L*g*m*cos(qA + qB) + 3*L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + 2*L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + L*L*m*qB_dot*qB_dot*cos(qB)*sin(qB) - 2*L*g*m*cos(qA)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));

x_dot = [qA_dot;qB_dot;qA_ddot;qB_ddot];

A = jacobian(x_dot,x);
B = jacobian(x_dot,u);

%% general setup variables
h = 0.002; % s
sim_time = 0:h:0.25;
x0 = [0;0;0;0];

%%
x_test = [ 0;0; 0; 0];
u_test = [0; 0];

a_test = vpa(subs(A,[x;u],[x_test;u_test]),6)
b_test = vpa(subs(B,[x;u],[x_test;u_test]),6)
x_dot_test = vpa(subs(x_dot,[x;u],[x_test;u_test]),6)

%%
x_test = [ 0.012566;-0.012566; 6.286450; -6.286450];
% x_test = [0.00628315;-0.00628309;6.28305;-6.28278];
x_init_test = [0.00628315;-0.00628309;6.28305;-6.28278];

u_test = [9451.340333; 3150.037249];
% u_test = [9458.794556; 3152.724932];
u_init_test = [9458.794556; 3152.724932];

a_test = (subs(A,[x;u],[x_init_test;u_init_test]));
b_test = (subs(B,[x;u],[x_init_test;u_init_test]));
x_dot_init_test = vpa(subs(x_dot,[x;u],[x_init_test;u_init_test]),6)

x_dot_now = vpa(F_xy_lin(a_test,b_test,x_test,u_test,x_dot_init_test,x_init_test,u_init_test),6);
scaled_x_dot_now = vpa(x_dot_now*h,6);
x_lin_res = x_test - x_init_test + scaled_x_dot_now;  % main equation
display(vpa(x_dot_now',6))
display(vpa(scaled_x_dot_now',6))
display(vpa(x_lin_res',6))
% standard is qA 0.025139, qB -0.025139, qA_dot   12.568, qB_dot -12.5677
% looking for qA 0.018856, qB -0.018855, qA_dot 6.284949, qB_dot -6.284949

% IT IS USING x_test - x_init_test!!!!!
%% nonlinear version
clear x_nonlin;
x_nonlin(:,1) = x0;

for i = 1:length(sim_time)
    U_nonlin(:,i) = get_U(sim_time(i),x_nonlin(:,i));
    k_1 = F_xy(x_nonlin(:,i),U_nonlin(:,i));
    k_2 = F_xy(x_nonlin(:,i)+0.5*h*k_1,U_nonlin(:,i));
    k_3 = F_xy(x_nonlin(:,i)+0.5*h*k_2,U_nonlin(:,i));
    k_4 = F_xy(x_nonlin(:,i)+k_3*h,U_nonlin(:,i));
    x_nonlin(:,i+1) = x_nonlin(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end
close all;
figure(1)
plot(sim_time,x_nonlin(1:2,(1:end-1)),'--'); hold on;
plot(sim_time,[sin(2*pi*sim_time);-sin(2*pi*sim_time)],'k')

% visualize(x_nonlin(1:2,:));

%% linearized version
clear x_lin;
x_lin(:,1) = x0;

A_num = vpa(subs(A,[x;u],[x_lin(:,1);[0;0]]),5);
B_num = vpa(subs(B,[x;u],[x_lin(:,1);[0;0]]),5);
x_dot_init = F_xy(x_lin(:,1),get_U(0,x_lin(:,1)));

integrator = "rk";

for i = 1:length(sim_time)
    U_lin(:,i) = U_nonlin(:,i);%get_U(sim_time(i),x_lin(:,i));
    if integrator == "rk"
        k_1 = F_xy_lin(A_num,B_num,x_lin(:,i),U_lin(:,i),x_dot_init,x_lin(:,1),U_lin(:,1));
        k_2 = F_xy_lin(A_num,B_num,x_lin(:,i)+0.5*h*k_1,U_lin(:,i),x_dot_init,x_lin(:,1),U_lin(:,1));
        k_3 = F_xy_lin(A_num,B_num,x_lin(:,i)+0.5*h*k_2,U_lin(:,i),x_dot_init,x_lin(:,1),U_lin(:,1));
        k_4 = F_xy_lin(A_num,B_num,x_lin(:,i)+k_3*h,U_lin(:,i),x_dot_init,x_lin(:,1),U_lin(:,1));
        x_lin(:,i+1) = x_lin(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    else
        x_lin(:,i+1) = x_lin(:,i) + F_xy_lin(A_num,B_num,x_lin(:,i),U_lin(:,i),x_dot_init,x_lin(:,1),U_lin(:,1))*h;  % main equation
    end
end
% figure(2)
plot(sim_time,x_lin(1:2,(1:end-1)));
xlabel('Time (s)');
ylabel('position (rad)');

%%
function x_dot = F_xy(X,U)
    L = 1; m = 1; g = 9.81;
    qA = X(1);
    qB = X(2);
    qA_dot = X(3);
    qB_dot = X(4);
    TA = U(1);
    TB = U(2);

    qA_ddot = -(TA - TB - TB*cos(qB) + L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));
    qB_ddot = (TA - 3*TB + TA*cos(qB) - 2*TB*cos(qB) + 2*L*g*m*cos(qA + qB) + 3*L*L*m*qA_dot*qA_dot*sin(qB) + L*L*m*qB_dot*qB_dot*sin(qB) - 2*L*g*m*cos(qA) + 2*L*L*m*qA_dot*qA_dot*cos(qB)*sin(qB) + L*L*m*qB_dot*qB_dot*cos(qB)*sin(qB) - 2*L*g*m*cos(qA)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*sin(qB) + L*g*m*cos(qA + qB)*cos(qB) + 2*L*L*m*qA_dot*qB_dot*cos(qB)*sin(qB))/(L*L*m*(cos(qB)*cos(qB) - 2));

    x_dot = [qA_dot;qB_dot;qA_ddot;qB_ddot];
end

function x_dot_lin = F_xy_lin(A,B,X,U,x_dot_init,x_init,u_init)
    x_dot_lin = A*(X-x_init)+B*(U-u_init) + x_dot_init;
end

function U = get_U(t,x)
    Kp = [1000, 0; 0, 1000];
    Kd = [ 10, 0; 0, 10];
    x_des = [sin(2*pi*t); -sin(2*pi*t)]; 
    x_dot_des = [2*pi*cos(2*pi*t); -2*pi*(cos(2*pi*t))];
    U = Kp*(x_des-x(1:2))+Kd*(x_dot_des-x(3:4));

%     U(1,1) = cos(t);
%     U(2,1) = cos(t);
end