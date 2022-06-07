#pragma once

#include <Mahi/Mpc/LinearModel.hpp>
#include <iostream>

class DoublePendulumModel : public LinearModel
{
private:
    double L = 1;
    double m = 1;
    double g = 9.81;
public:
    DoublePendulumModel(){}

    Eigen::MatrixXd get_A(const std::vector<double>& x,const std::vector<double>& u) override {
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4,4);

        const double qA = x[0];
        const double qB = x[1];
        const double qA_dot = x[2];
        const double qB_dot = x[3];

        const double TA = u[0];
        const double TB = u[1];

        const double t2 = cos(qB);
        const double t3 = sin(qA);
        const double t4 = sin(qB);
        const double t5 = qA+qB;
        const double t6 = L*L;
        const double t7 = qA_dot*2.0;
        const double t8 = qB*2.0;
        const double t9 = qB_dot*2.0;
        const double t10 = qA_dot*qA_dot;
        const double t11 = qB_dot*qB_dot;
        const double t17 = 1.0/L;
        const double t20 = 1.0/m;
        const double t12 = cos(t8);
        const double t13 = t2*t2;
        const double t14 = t2*t2*t2;
        const double t15 = sin(t8);
        const double t16 = t4*t4*t4;
        const double t18 = 1.0/t6;
        const double t19 = sin(t5);
        const double t22 = t7+t9;
        const double t23 = L*g*m*t3*2.0;
        const double t21 = TA*t15;
        const double t24 = t12-3.0;
        const double t25 = t13-2.0;
        const double t30 = L*g*m*t3*t13*3.0;
        const double t31 = m*t6*t11*t14;
        const double t32 = m*qB_dot*t6*t7*t14;
        const double t26 = -t21;
        const double t27 = 1.0/t24;
        const double t28 = 1.0/t25;
        const double t33 = -t30;
        const double t29 = t28*t28;
        A(0,2) = 1.0;
        A(1,3) = 1.0;
        A(2,0) = -g*t17*t27*(t3*3.0-sin(qB+t5));
        A(2,1) = t18*t20*t29*(t23+t26+t31+t32+t33+TB*t4*3.0+TB*t15-TB*t16-m*t6*t10*2.0+m*t6*t10*t13*3.0+m*t6*t10*t14);
        A(2,2) = -t28*(t4*t7+t4*t9+t2*t4*t7);
        A(2,3) = (t4*t22)/(t4*t4+1.0);
        A(3,0) = g*t17*t28*(t3*2.0-t19*2.0+t2*t3*2.0-t2*t19);
        A(3,1) = -t18*t20*t29*(t23+t26+t31+t32+t33-TA*t4*3.0+TA*t16+TB*t4*6.0+TB*t15*3.0-TB*t16*2.0-m*t6*t10*4.0-m*t6*t11*2.0-m*qA_dot*qB_dot*t6*4.0+m*t6*t10*t13*6.0+m*t6*t10*t14*3.0+m*t6*t11*t13*3.0-L*g*m*t3*t14*2.0+m*qA_dot*qB_dot*t6*t13*6.0);
        A(3,2) = t27*(qA_dot*t4*1.2E+1+qA_dot*t15*4.0+qB_dot*t4*4.0+t9*t15);
        A(3,3) = t22*t28*(t4+t2*t4);

        return A;
    }

    Eigen::MatrixXd get_B(const std::vector<double>& x, const std::vector<double>& u) override {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(4,2);

        const double qA = x[0];
        const double qB = x[1];
        const double qA_dot = x[2];
        const double qB_dot = x[3];

        const double t2 = cos(qB);
        const double t5 = 1.0/(L*L);
        const double t6 = 1.0/m;
        const double t3 = t2+1.0;
        const double t4 = t2*t2;
        const double t7 = t4-2.0;
        const double t8 = 1.0/t7;
        const double t9 = t3*t5*t6*t8;

        B(2,0) = -t5*t6*t8;
        B(2,1) = t9;
        B(3,0) = t9;
        B(3,1) = -t5*t6*t8*(t3*2.0+1.0);
        
        return B;
    }

    Eigen::MatrixXd get_Xdot(const std::vector<double>& x, const std::vector<double>& u) override {
        Eigen::MatrixXd Xdot = Eigen::MatrixXd::Zero(4,1);

        const double qA = x[0];
        const double qB = x[1];
        const double qA_dot = x[2];
        const double qB_dot = x[3];

        const double TA = u[0];
        const double TB = u[1];

        const double t2 = cos(qA);
        const double t3 = cos(qB);
        const double t4 = sin(qB);
        const double t5 = qA+qB;
        const double t6 = L*L;
        const double t7 = qA_dot*qA_dot;
        const double t8 = qB_dot*qB_dot;
        const double t12 = 1.0/m;
        const double t9 = t3*t3;
        const double t10 = cos(t5);
        const double t11 = 1.0/t6;
        const double t13 = L*g*m*t2*2.0;
        const double t16 = m*t4*t6*t8;
        const double t17 = m*qA_dot*qB_dot*t4*t6*2.0;
        const double t14 = t9-2.0;
        const double t15 = -t13;
        const double t19 = L*g*m*t3*t10;
        const double t18 = 1.0/t14;
        Xdot(0,0) = qA_dot;
        Xdot(1,0) = qB_dot;
        Xdot(2,0) = -t11*t12*t18*(TA-TB+t15+t16+t17+t19-TB*t3+m*t4*t6*t7+m*t3*t4*t6*t7);
        Xdot(3,0) = t11*t12*t18*(TA-TB*3.0+t15+t16+t17+t19+TA*t3-TB*t3*2.0+t3*t16+t3*t17+m*t4*t6*t7*3.0+L*g*m*t10*2.0-L*g*m*t2*t3*2.0+m*t3*t4*t6*t7*2.0);
        
        return Xdot;
    }

};