#include <iostream>
// #include <casadi/casadi.hpp>
#include <C:/MatlabHelpers/Casadi/casadi-windows-matlabR2016a-v3.5.5/include/casadi/casadi.hpp>
#include <Mahi/Util/Print.hpp>


int main() {
    casadi::SX q1 = casadi::SX::sym("q1");
    casadi::SX q2 = casadi::SX::sym("q2");
    casadi::SX q3 = casadi::SX::sym("q3");
    casadi::SX q4 = casadi::SX::sym("q4");

    casadi::SX Icxx0 = casadi::SX::sym("Icxx0");
    casadi::SX Icyy0 = casadi::SX::sym("Icyy0");
    casadi::SX Iczz0 = casadi::SX::sym("Iczz0");
    casadi::SX Icxy0 = casadi::SX::sym("Icxy0");
    casadi::SX Icxz0 = casadi::SX::sym("Icxz0");
    casadi::SX Icyz0 = casadi::SX::sym("Icyz0");
    casadi::SX Pcx0 = casadi::SX::sym("Pcx0");
    casadi::SX Pcy0 = casadi::SX::sym("Pcy0");
    casadi::SX Pcz0 = casadi::SX::sym("Pcz0");

    casadi::SX Icxx1 = casadi::SX::sym("Icxx1");
    casadi::SX Icyy1 = casadi::SX::sym("Icyy1");
    casadi::SX Iczz1 = casadi::SX::sym("Iczz1");
    casadi::SX Icxy1 = casadi::SX::sym("Icxy1");
    casadi::SX Icxz1 = casadi::SX::sym("Icxz1");
    casadi::SX Icyz1 = casadi::SX::sym("Icyz1");
    casadi::SX Pcx1 = casadi::SX::sym("Pcx1");
    casadi::SX Pcy1 = casadi::SX::sym("Pcy1");
    casadi::SX Pcz1 = casadi::SX::sym("Pcz1");

    casadi::SX Icxx2 = casadi::SX::sym("Icxx2");
    casadi::SX Icyy2 = casadi::SX::sym("Icyy2");
    casadi::SX Iczz2 = casadi::SX::sym("Iczz2");
    casadi::SX Icxy2 = casadi::SX::sym("Icxy2");
    casadi::SX Icxz2 = casadi::SX::sym("Icxz2");
    casadi::SX Icyz2 = casadi::SX::sym("Icyz2");
    casadi::SX Pcx2 = casadi::SX::sym("Pcx2");
    casadi::SX Pcy2 = casadi::SX::sym("Pcy2");
    casadi::SX Pcz2 = casadi::SX::sym("Pcz2");

    casadi::SX Icxx3 = casadi::SX::sym("Icxx3");
    casadi::SX Icyy3 = casadi::SX::sym("Icyy3");
    casadi::SX Iczz3 = casadi::SX::sym("Iczz3");
    casadi::SX Icxy3 = casadi::SX::sym("Icxy3");
    casadi::SX Icxz3 = casadi::SX::sym("Icxz3");
    casadi::SX Icyz3 = casadi::SX::sym("Icyz3");
    casadi::SX Pcx3 = casadi::SX::sym("Pcx3");
    casadi::SX Pcy3 = casadi::SX::sym("Pcy3");
    casadi::SX Pcz3 = casadi::SX::sym("Pcz3");

    casadi::SX m0 = casadi::SX::sym("m0");
    casadi::SX m1 = casadi::SX::sym("m1");
    casadi::SX m2 = casadi::SX::sym("m2");
    casadi::SX m3 = casadi::SX::sym("m3");


    casadi::SX M00 = Icxx3+Icyy1+Iczz0+Iczz2+m1*(9.0/4.0E+2)+m2*(9.0/4.0E+2)+m3*(9.0/4.0E+2)+(Pcx0*Pcx0)*m0+(Pcx1*Pcx1)*m1+(Pcx2*Pcx2)*m2+(Pcy0*Pcy0)*m0+(Pcy2*Pcy2)*m2+(Pcy3*Pcy3)*m3+(Pcz1*Pcz1)*m1+(Pcz3*Pcz3)*m3+Icxx1*pow(cos(q1),2.0)-Icxx3*pow(cos(q1),2.0)-Icxx3*pow(cos(q3),2.0)-Icyy1*pow(cos(q1),2.0)+Icyy2*pow(cos(q1),2.0)+Icyy3*pow(cos(q3),2.0)-Iczz2*pow(cos(q1),2.0)+Iczz3*pow(cos(q1),2.0)+Icxy1*sin(q1*2.0)-Icxy3*sin(q3*2.0)-Pcx1*m1*sin(q1)*(3.0/1.0E+1)+Icxx2*pow(cos(q1),2.0)*pow(cos(q2),2.0)+Icxx3*pow(cos(q1),2.0)*pow(cos(q3),2.0)-Icyy2*pow(cos(q1),2.0)*pow(cos(q2),2.0)+Icyy3*pow(cos(q1),2.0)*pow(cos(q2),2.0)-Icyy3*pow(cos(q1),2.0)*pow(cos(q3),2.0)-Iczz3*pow(cos(q1),2.0)*pow(cos(q2),2.0)-(Pcx1*Pcx1)*m1*pow(cos(q1),2.0)+(Pcx3*Pcx3)*m3*pow(cos(q1),2.0)+(Pcx3*Pcx3)*m3*pow(cos(q3),2.0)+(Pcy1*Pcy1)*m1*pow(cos(q1),2.0)-(Pcy2*Pcy2)*m2*pow(cos(q1),2.0)-(Pcy3*Pcy3)*m3*pow(cos(q3),2.0)+(Pcz2*Pcz2)*m2*pow(cos(q1),2.0)-(Pcz3*Pcz3)*m3*pow(cos(q1),2.0)-Pcy1*m1*cos(q1)*(3.0/1.0E+1)+Pcz2*m2*cos(q1)*(3.0/1.0E+1)+Icxy3*cos(q1)*cos(q2)*sin(q1)*2.0-Icxz2*cos(q1)*cos(q2)*sin(q1)*2.0+Icxx3*pow(cos(q1),2.0)*pow(cos(q2),2.0)*pow(cos(q3),2.0)-Icyy3*pow(cos(q1),2.0)*pow(cos(q2),2.0)*pow(cos(q3),2.0)+Icyz2*cos(q1)*sin(q1)*sin(q2)*2.0+Icxy2*pow(cos(q1),2.0)*cos(q2)*sin(q2)*2.0+Icxy3*pow(cos(q1),2.0)*cos(q3)*sin(q3)*2.0+Pcx1*Pcy1*m1*sin(q1*2.0)-Pcx3*Pcy3*m3*sin(q3*2.0)+Pcy3*m3*cos(q1)*cos(q3)*(3.0/1.0E+1)-Pcx2*m2*cos(q2)*sin(q1)*(3.0/1.0E+1)+Pcx3*m3*cos(q1)*sin(q3)*(3.0/1.0E+1)-(Pcx2*Pcx2)*m2*pow(cos(q1),2.0)*pow(cos(q2),2.0)-(Pcx3*Pcx3)*m3*pow(cos(q1),2.0)*pow(cos(q3),2.0)+(Pcy2*Pcy2)*m2*pow(cos(q1),2.0)*pow(cos(q2),2.0)-(Pcy3*Pcy3)*m3*pow(cos(q1),2.0)*pow(cos(q2),2.0)+(Pcy3*Pcy3)*m3*pow(cos(q1),2.0)*pow(cos(q3),2.0)+(Pcz3*Pcz3)*m3*pow(cos(q1),2.0)*pow(cos(q2),2.0)+Pcy2*m2*sin(q1)*sin(q2)*(3.0/1.0E+1)-Pcz3*m3*sin(q1)*sin(q2)*(3.0/1.0E+1)-(Pcx3*Pcx3)*m3*pow(cos(q1),2.0)*pow(cos(q2),2.0)*pow(cos(q3),2.0)+(Pcy3*Pcy3)*m3*pow(cos(q1),2.0)*pow(cos(q2),2.0)*pow(cos(q3),2.0)-Icyz3*cos(q1)*cos(q3)*sin(q1)*sin(q2)*2.0-Icxz3*cos(q1)*sin(q1)*sin(q2)*sin(q3)*2.0-Icxy3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q1)*4.0-Icxz3*pow(cos(q1),2.0)*cos(q2)*cos(q3)*sin(q2)*2.0+Icyz3*pow(cos(q1),2.0)*cos(q2)*sin(q2)*sin(q3)*2.0-Pcx3*m3*cos(q2)*cos(q3)*sin(q1)*(3.0/1.0E+1)+Pcy3*m3*cos(q2)*sin(q1)*sin(q3)*(3.0/1.0E+1)+Icxy3*pow(cos(q1),2.0)*pow(cos(q2),2.0)*cos(q3)*sin(q3)*2.0+Pcy2*Pcz2*m2*cos(q1)*sin(q1)*sin(q2)*2.0+Pcx2*Pcy2*m2*pow(cos(q1),2.0)*cos(q2)*sin(q2)*2.0+Pcx3*Pcy3*m3*pow(cos(q1),2.0)*cos(q3)*sin(q3)*2.0+Icxx3*cos(q1)*cos(q2)*cos(q3)*sin(q1)*sin(q3)*2.0-Icyy3*cos(q1)*cos(q2)*cos(q3)*sin(q1)*sin(q3)*2.0+Pcx3*Pcy3*m3*cos(q1)*cos(q2)*sin(q1)*2.0-Pcx2*Pcz2*m2*cos(q1)*cos(q2)*sin(q1)*2.0-Pcx3*Pcy3*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q1)*4.0-Pcx3*Pcz3*m3*pow(cos(q1),2.0)*cos(q2)*cos(q3)*sin(q2)*2.0+Pcy3*Pcz3*m3*pow(cos(q1),2.0)*cos(q2)*sin(q2)*sin(q3)*2.0+Pcx3*Pcy3*m3*pow(cos(q1),2.0)*pow(cos(q2),2.0)*cos(q3)*sin(q3)*2.0-(Pcx3*Pcx3)*m3*cos(q1)*cos(q2)*cos(q3)*sin(q1)*sin(q3)*2.0+(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*cos(q3)*sin(q1)*sin(q3)*2.0-Pcy3*Pcz3*m3*cos(q1)*cos(q3)*sin(q1)*sin(q2)*2.0-Pcx3*Pcz3*m3*cos(q1)*sin(q1)*sin(q2)*sin(q3)*2.0;
    casadi::SX M01 = Icxy2*cos(q1)-Icxz1*cos(q1)+Icyz1*sin(q1)-Icxz3*cos(q1)*cos(q3)-Icyz2*cos(q2)*sin(q1)+Icyz3*cos(q1)*sin(q3)+Icxy3*sin(q1)*sin(q2)-Icxz2*sin(q1)*sin(q2)-Icxy2*cos(q1)*pow(cos(q2),2.0)*2.0+Icxx2*cos(q1)*cos(q2)*sin(q2)-Icyy2*cos(q1)*cos(q2)*sin(q2)+Icyy3*cos(q1)*cos(q2)*sin(q2)+Icyz3*cos(q2)*cos(q3)*sin(q1)-Iczz3*cos(q1)*cos(q2)*sin(q2)+Pcx2*Pcy2*m2*cos(q1)-Pcx1*Pcz1*m1*cos(q1)+Icxz3*cos(q2)*sin(q1)*sin(q3)+Pcy1*Pcz1*m1*sin(q1)+Icxz3*cos(q1)*pow(cos(q2),2.0)*cos(q3)*2.0-Icyz3*cos(q1)*pow(cos(q2),2.0)*sin(q3)*2.0-Icxy3*pow(cos(q3),2.0)*sin(q1)*sin(q2)*2.0-(Pcx2*Pcx2)*m2*cos(q1)*cos(q2)*sin(q2)+(Pcy2*Pcy2)*m2*cos(q1)*cos(q2)*sin(q2)-(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*sin(q2)+(Pcz3*Pcz3)*m3*cos(q1)*cos(q2)*sin(q2)-Pcx3*Pcz3*m3*cos(q1)*cos(q3)-Pcy2*Pcz2*m2*cos(q2)*sin(q1)+Pcy3*Pcz3*m3*cos(q1)*sin(q3)+Icxx3*cos(q3)*sin(q1)*sin(q2)*sin(q3)-Icyy3*cos(q3)*sin(q1)*sin(q2)*sin(q3)+Pcx3*Pcy3*m3*sin(q1)*sin(q2)-Pcx2*Pcz2*m2*sin(q1)*sin(q2)+Icxx3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)-Icyy3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)-Pcx2*Pcy2*m2*cos(q1)*pow(cos(q2),2.0)*2.0+Pcx3*Pcz3*m3*cos(q2)*sin(q1)*sin(q3)-(Pcx3*Pcx3)*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)+(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)+Pcx3*Pcz3*m3*cos(q1)*pow(cos(q2),2.0)*cos(q3)*2.0-Pcy3*Pcz3*m3*cos(q1)*pow(cos(q2),2.0)*sin(q3)*2.0-Pcx3*Pcy3*m3*pow(cos(q3),2.0)*sin(q1)*sin(q2)*2.0+Icxy3*cos(q1)*cos(q2)*cos(q3)*sin(q2)*sin(q3)*2.0+Pcy3*Pcz3*m3*cos(q2)*cos(q3)*sin(q1)-(Pcx3*Pcx3)*m3*cos(q3)*sin(q1)*sin(q2)*sin(q3)+(Pcy3*Pcy3)*m3*cos(q3)*sin(q1)*sin(q2)*sin(q3)+Pcx3*Pcy3*m3*cos(q1)*cos(q2)*cos(q3)*sin(q2)*sin(q3)*2.0;
    casadi::SX M02 = Icxx3*sin(q1)+Iczz2*sin(q1)+Pcy2*m2*sin(q2)*(3.0/2.0E+1)-Pcz3*m3*sin(q2)*(3.0/2.0E+1)+Icxy3*cos(q1)*cos(q2)-Icxz2*cos(q1)*cos(q2)+(Pcx2*Pcx2)*m2*sin(q1)+(Pcy2*Pcy2)*m2*sin(q1)+(Pcy3*Pcy3)*m3*sin(q1)+(Pcz3*Pcz3)*m3*sin(q1)+Icyz2*cos(q1)*sin(q2)-Icxx3*pow(cos(q3),2.0)*sin(q1)+Icyy3*pow(cos(q3),2.0)*sin(q1)-Pcx2*m2*cos(q2)*(3.0/2.0E+1)-Icyz3*cos(q1)*cos(q3)*sin(q2)-Icxy3*cos(q3)*sin(q1)*sin(q3)*2.0-Icxz3*cos(q1)*sin(q2)*sin(q3)-Icxy3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*2.0+(Pcx3*Pcx3)*m3*pow(cos(q3),2.0)*sin(q1)-(Pcy3*Pcy3)*m3*pow(cos(q3),2.0)*sin(q1)-Pcx3*m3*cos(q2)*cos(q3)*(3.0/2.0E+1)+Pcy3*m3*cos(q2)*sin(q3)*(3.0/2.0E+1)+Icxx3*cos(q1)*cos(q2)*cos(q3)*sin(q3)-Icyy3*cos(q1)*cos(q2)*cos(q3)*sin(q3)+Pcx3*Pcy3*m3*cos(q1)*cos(q2)-Pcx2*Pcz2*m2*cos(q1)*cos(q2)+Pcy2*Pcz2*m2*cos(q1)*sin(q2)-Pcx3*Pcy3*m3*cos(q3)*sin(q1)*sin(q3)*2.0-Pcx3*Pcz3*m3*cos(q1)*sin(q2)*sin(q3)-Pcx3*Pcy3*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*2.0-(Pcx3*Pcx3)*m3*cos(q1)*cos(q2)*cos(q3)*sin(q3)+(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*cos(q3)*sin(q3)-Pcy3*Pcz3*m3*cos(q1)*cos(q3)*sin(q2);
    casadi::SX M03 = -Icyz3*cos(q3)*sin(q1)+Iczz3*cos(q1)*sin(q2)-Icxz3*sin(q1)*sin(q3)-Icxz3*cos(q1)*cos(q2)*cos(q3)+(Pcx3*Pcx3)*m3*cos(q1)*sin(q2)+(Pcy3*Pcy3)*m3*cos(q1)*sin(q2)+Icyz3*cos(q1)*cos(q2)*sin(q3)+Pcy3*m3*cos(q3)*sin(q2)*(3.0/2.0E+1)+Pcx3*m3*sin(q2)*sin(q3)*(3.0/2.0E+1)-Pcy3*Pcz3*m3*cos(q3)*sin(q1)-Pcx3*Pcz3*m3*sin(q1)*sin(q3)-Pcx3*Pcz3*m3*cos(q1)*cos(q2)*cos(q3)+Pcy3*Pcz3*m3*cos(q1)*cos(q2)*sin(q3);
    casadi::SX M10 = Icxy2*cos(q1)-Icxz1*cos(q1)+Icyz1*sin(q1)-Icxz3*cos(q1)*cos(q3)-Icyz2*cos(q2)*sin(q1)+Icyz3*cos(q1)*sin(q3)+Icxy3*sin(q1)*sin(q2)-Icxz2*sin(q1)*sin(q2)-Icxy2*cos(q1)*pow(cos(q2),2.0)*2.0+Icxx2*cos(q1)*cos(q2)*sin(q2)-Icyy2*cos(q1)*cos(q2)*sin(q2)+Icyy3*cos(q1)*cos(q2)*sin(q2)+Icyz3*cos(q2)*cos(q3)*sin(q1)-Iczz3*cos(q1)*cos(q2)*sin(q2)+Pcx2*Pcy2*m2*cos(q1)-Pcx1*Pcz1*m1*cos(q1)+Icxz3*cos(q2)*sin(q1)*sin(q3)+Pcy1*Pcz1*m1*sin(q1)+Icxz3*cos(q1)*pow(cos(q2),2.0)*cos(q3)*2.0-Icyz3*cos(q1)*pow(cos(q2),2.0)*sin(q3)*2.0-Icxy3*pow(cos(q3),2.0)*sin(q1)*sin(q2)*2.0-(Pcx2*Pcx2)*m2*cos(q1)*cos(q2)*sin(q2)+(Pcy2*Pcy2)*m2*cos(q1)*cos(q2)*sin(q2)-(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*sin(q2)+(Pcz3*Pcz3)*m3*cos(q1)*cos(q2)*sin(q2)-Pcx3*Pcz3*m3*cos(q1)*cos(q3)-Pcy2*Pcz2*m2*cos(q2)*sin(q1)+Pcy3*Pcz3*m3*cos(q1)*sin(q3)+Icxx3*cos(q3)*sin(q1)*sin(q2)*sin(q3)-Icyy3*cos(q3)*sin(q1)*sin(q2)*sin(q3)+Pcx3*Pcy3*m3*sin(q1)*sin(q2)-Pcx2*Pcz2*m2*sin(q1)*sin(q2)+Icxx3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)-Icyy3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)-Pcx2*Pcy2*m2*cos(q1)*pow(cos(q2),2.0)*2.0+Pcx3*Pcz3*m3*cos(q2)*sin(q1)*sin(q3)-(Pcx3*Pcx3)*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)+(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*sin(q2)+Pcx3*Pcz3*m3*cos(q1)*pow(cos(q2),2.0)*cos(q3)*2.0-Pcy3*Pcz3*m3*cos(q1)*pow(cos(q2),2.0)*sin(q3)*2.0-Pcx3*Pcy3*m3*pow(cos(q3),2.0)*sin(q1)*sin(q2)*2.0+Icxy3*cos(q1)*cos(q2)*cos(q3)*sin(q2)*sin(q3)*2.0+Pcy3*Pcz3*m3*cos(q2)*cos(q3)*sin(q1)-(Pcx3*Pcx3)*m3*cos(q3)*sin(q1)*sin(q2)*sin(q3)+(Pcy3*Pcy3)*m3*cos(q3)*sin(q1)*sin(q2)*sin(q3)+Pcx3*Pcy3*m3*cos(q1)*cos(q2)*cos(q3)*sin(q2)*sin(q3)*2.0;
    casadi::SX M11 = Icxx2+Icyy3+Iczz1+(Pcx1*Pcx1)*m1+(Pcx3*Pcx3)*m3+(Pcy1*Pcy1)*m1+(Pcy2*Pcy2)*m2+(Pcz2*Pcz2)*m2+(Pcz3*Pcz3)*m3-Icxx2*pow(cos(q2),2.0)+Icxx3*pow(cos(q3),2.0)+Icyy2*pow(cos(q2),2.0)-Icyy3*pow(cos(q2),2.0)-Icyy3*pow(cos(q3),2.0)+Iczz3*pow(cos(q2),2.0)-Icxy2*sin(q2*2.0)+Icxy3*sin(q3*2.0)-Icxx3*pow(cos(q2),2.0)*pow(cos(q3),2.0)+Icyy3*pow(cos(q2),2.0)*pow(cos(q3),2.0)+(Pcx2*Pcx2)*m2*pow(cos(q2),2.0)-(Pcx3*Pcx3)*m3*pow(cos(q3),2.0)-(Pcy2*Pcy2)*m2*pow(cos(q2),2.0)+(Pcy3*Pcy3)*m3*pow(cos(q2),2.0)+(Pcy3*Pcy3)*m3*pow(cos(q3),2.0)-(Pcz3*Pcz3)*m3*pow(cos(q2),2.0)+Icxz3*cos(q2)*cos(q3)*sin(q2)*2.0-Icyz3*cos(q2)*sin(q2)*sin(q3)*2.0-Icxy3*pow(cos(q2),2.0)*cos(q3)*sin(q3)*2.0-Pcx2*Pcy2*m2*sin(q2*2.0)+Pcx3*Pcy3*m3*sin(q3*2.0)+(Pcx3*Pcx3)*m3*pow(cos(q2),2.0)*pow(cos(q3),2.0)-(Pcy3*Pcy3)*m3*pow(cos(q2),2.0)*pow(cos(q3),2.0)-Pcy3*Pcz3*m3*cos(q2)*sin(q2)*sin(q3)*2.0-Pcx3*Pcy3*m3*pow(cos(q2),2.0)*cos(q3)*sin(q3)*2.0+Pcx3*Pcz3*m3*cos(q2)*cos(q3)*sin(q2)*2.0;
    casadi::SX M12 = -Icyz2*cos(q2)+Icxy3*sin(q2)-Icxz2*sin(q2)+Icyz3*cos(q2)*cos(q3)+Icxz3*cos(q2)*sin(q3)-Icxy3*pow(cos(q3),2.0)*sin(q2)*2.0-Pcy2*Pcz2*m2*cos(q2)+Icxx3*cos(q3)*sin(q2)*sin(q3)-Icyy3*cos(q3)*sin(q2)*sin(q3)+Pcx3*Pcy3*m3*sin(q2)-Pcx2*Pcz2*m2*sin(q2)+Pcy3*Pcz3*m3*cos(q2)*cos(q3)-(Pcx3*Pcx3)*m3*cos(q3)*sin(q2)*sin(q3)+(Pcy3*Pcy3)*m3*cos(q3)*sin(q2)*sin(q3)+Pcx3*Pcz3*m3*cos(q2)*sin(q3)-Pcx3*Pcy3*m3*pow(cos(q3),2.0)*sin(q2)*2.0;
    casadi::SX M13 = -Iczz3*cos(q2)-(Pcx3*Pcx3)*m3*cos(q2)-(Pcy3*Pcy3)*m3*cos(q2)-Icxz3*cos(q3)*sin(q2)+Icyz3*sin(q2)*sin(q3)-Pcx3*Pcz3*m3*cos(q3)*sin(q2)+Pcy3*Pcz3*m3*sin(q2)*sin(q3);
    casadi::SX M20 = Icxx3*sin(q1)+Iczz2*sin(q1)+Pcy2*m2*sin(q2)*(3.0/2.0E+1)-Pcz3*m3*sin(q2)*(3.0/2.0E+1)+Icxy3*cos(q1)*cos(q2)-Icxz2*cos(q1)*cos(q2)+(Pcx2*Pcx2)*m2*sin(q1)+(Pcy2*Pcy2)*m2*sin(q1)+(Pcy3*Pcy3)*m3*sin(q1)+(Pcz3*Pcz3)*m3*sin(q1)+Icyz2*cos(q1)*sin(q2)-Icxx3*pow(cos(q3),2.0)*sin(q1)+Icyy3*pow(cos(q3),2.0)*sin(q1)-Pcx2*m2*cos(q2)*(3.0/2.0E+1)-Icyz3*cos(q1)*cos(q3)*sin(q2)-Icxy3*cos(q3)*sin(q1)*sin(q3)*2.0-Icxz3*cos(q1)*sin(q2)*sin(q3)-Icxy3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*2.0+(Pcx3*Pcx3)*m3*pow(cos(q3),2.0)*sin(q1)-(Pcy3*Pcy3)*m3*pow(cos(q3),2.0)*sin(q1)-Pcx3*m3*cos(q2)*cos(q3)*(3.0/2.0E+1)+Pcy3*m3*cos(q2)*sin(q3)*(3.0/2.0E+1)+Icxx3*cos(q1)*cos(q2)*cos(q3)*sin(q3)-Icyy3*cos(q1)*cos(q2)*cos(q3)*sin(q3)+Pcx3*Pcy3*m3*cos(q1)*cos(q2)-Pcx2*Pcz2*m2*cos(q1)*cos(q2)+Pcy2*Pcz2*m2*cos(q1)*sin(q2)-Pcx3*Pcy3*m3*cos(q3)*sin(q1)*sin(q3)*2.0-Pcx3*Pcz3*m3*cos(q1)*sin(q2)*sin(q3)-Pcx3*Pcy3*m3*cos(q1)*cos(q2)*pow(cos(q3),2.0)*2.0-(Pcx3*Pcx3)*m3*cos(q1)*cos(q2)*cos(q3)*sin(q3)+(Pcy3*Pcy3)*m3*cos(q1)*cos(q2)*cos(q3)*sin(q3)-Pcy3*Pcz3*m3*cos(q1)*cos(q3)*sin(q2);
    casadi::SX M21 = -Icyz2*cos(q2)+Icxy3*sin(q2)-Icxz2*sin(q2)+Icyz3*cos(q2)*cos(q3)+Icxz3*cos(q2)*sin(q3)-Icxy3*pow(cos(q3),2.0)*sin(q2)*2.0-Pcy2*Pcz2*m2*cos(q2)+Icxx3*cos(q3)*sin(q2)*sin(q3)-Icyy3*cos(q3)*sin(q2)*sin(q3)+Pcx3*Pcy3*m3*sin(q2)-Pcx2*Pcz2*m2*sin(q2)+Pcy3*Pcz3*m3*cos(q2)*cos(q3)-(Pcx3*Pcx3)*m3*cos(q3)*sin(q2)*sin(q3)+(Pcy3*Pcy3)*m3*cos(q3)*sin(q2)*sin(q3)+Pcx3*Pcz3*m3*cos(q2)*sin(q3)-Pcx3*Pcy3*m3*pow(cos(q3),2.0)*sin(q2)*2.0;
    casadi::SX M22 = Icxx3+Iczz2+(Pcx2*Pcx2)*m2+(Pcy2*Pcy2)*m2+(Pcy3*Pcy3)*m3+(Pcz3*Pcz3)*m3-Icxx3*pow(cos(q3),2.0)+Icyy3*pow(cos(q3),2.0)-Icxy3*sin(q3*2.0)+(Pcx3*Pcx3)*m3*pow(cos(q3),2.0)-(Pcy3*Pcy3)*m3*pow(cos(q3),2.0)-Pcx3*Pcy3*m3*sin(q3*2.0);
    casadi::SX M23 = -Icyz3*cos(q3)-Icxz3*sin(q3)-Pcy3*Pcz3*m3*cos(q3)-Pcx3*Pcz3*m3*sin(q3);
    casadi::SX M30 = -Icyz3*cos(q3)*sin(q1)+Iczz3*cos(q1)*sin(q2)-Icxz3*sin(q1)*sin(q3)-Icxz3*cos(q1)*cos(q2)*cos(q3)+(Pcx3*Pcx3)*m3*cos(q1)*sin(q2)+(Pcy3*Pcy3)*m3*cos(q1)*sin(q2)+Icyz3*cos(q1)*cos(q2)*sin(q3)+Pcy3*m3*cos(q3)*sin(q2)*(3.0/2.0E+1)+Pcx3*m3*sin(q2)*sin(q3)*(3.0/2.0E+1)-Pcy3*Pcz3*m3*cos(q3)*sin(q1)-Pcx3*Pcz3*m3*sin(q1)*sin(q3)-Pcx3*Pcz3*m3*cos(q1)*cos(q2)*cos(q3)+Pcy3*Pcz3*m3*cos(q1)*cos(q2)*sin(q3);
    casadi::SX M31 = -Iczz3*cos(q2)-(Pcx3*Pcx3)*m3*cos(q2)-(Pcy3*Pcy3)*m3*cos(q2)-Icxz3*cos(q3)*sin(q2)+Icyz3*sin(q2)*sin(q3)-Pcx3*Pcz3*m3*cos(q3)*sin(q2)+Pcy3*Pcz3*m3*sin(q2)*sin(q3);
    casadi::SX M32 = -Icyz3*cos(q3)-Icxz3*sin(q3)-Pcy3*Pcz3*m3*cos(q3)-Pcx3*Pcz3*m3*sin(q3);
    casadi::SX M33 = Iczz3+(Pcx3*Pcx3)*m3+(Pcy3*Pcy3)*m3;


    casadi::SX M_row0 = casadi::SX::horzcat({M00,M01,M02,M03});
    casadi::SX M_row1 = casadi::SX::horzcat({M10,M11,M12,M13});
    casadi::SX M_row2 = casadi::SX::horzcat({M20,M21,M22,M23});
    casadi::SX M_row3 = casadi::SX::horzcat({M30,M31,M32,M33});


    casadi::SX M = casadi::SX::vertcat({M_row0,M_row1,M_row2,M_row3});

    //auto M_inv = casadi::SX::inv(M);

    casadi::SX tempVar = casadi::SX::sym("tempVar",4,4);
    casadi::Function F = casadi::Function("F",{tempVar},{casadi::SX::inv(tempVar)},{"M"},{"M_inv"});


    auto invRes =\
    F({M});
    // auto invResRes = invRes["M_inv"];
    std::cout << invRes << std::endl;



    return 0;
}