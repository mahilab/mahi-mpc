#pragma once
#include <vector>
#include <casadi/casadi.hpp>
#include <Mahi/Casadi/mass_properties.hpp>
using namespace casadi;

SX get_G(const SX& qs){
	SX G = SX::zeros(4,1); 

	auto q0 = qs(0,0);
	auto q1 = qs(1,0);
	auto q2 = qs(2,0);
	auto q3 = qs(3,0);
	auto qd0 = qs(4,0);
	auto qd1 = qs(5,0);
	auto qd2 = qs(6,0);
	auto qd3 = qs(7,0);

	auto t2 = cos(q0);
	auto t3 = cos(q1);
	auto t4 = cos(q2);
	auto t5 = cos(q3);
	auto t6 = sin(q0);
	auto t7 = sin(q1);
	auto t8 = sin(q2);
	auto t9 = sin(q3);
	G(0,0) = g*(m1*t6*(-3.0/2.0E+1)-m2*t6*(3.0/2.0E+1)-m3*t6*(3.0/2.0E+1)+Pcx0*m0*t6+Pcy0*m0*t2+Pcz1*m1*t2+Pcx2*m2*t2*t8+Pcx1*m1*t6*t7+Pcy2*m2*t2*t4+Pcy1*m1*t3*t6-Pcz3*m3*t2*t4-Pcz2*m2*t3*t6+Pcx2*m2*t4*t6*t7+Pcx3*m3*t2*t5*t8-Pcx3*m3*t3*t6*t9-Pcy3*m3*t3*t5*t6-Pcy2*m2*t6*t7*t8-Pcy3*m3*t2*t8*t9+Pcz3*m3*t6*t7*t8+Pcx3*m3*t4*t5*t6*t7-Pcy3*m3*t4*t6*t7*t9);
	G(1,0) = -g*t2*(Pcx1*m1*t3-Pcy1*m1*t7+Pcz2*m2*t7+Pcx2*m2*t3*t4+Pcx3*m3*t7*t9-Pcy2*m2*t3*t8+Pcy3*m3*t5*t7+Pcz3*m3*t3*t8+Pcx3*m3*t3*t4*t5-Pcy3*m3*t3*t4*t9);
	G(2,0) = g*(Pcx2*m2*t4*t6-Pcy2*m2*t6*t8+Pcz3*m3*t6*t8+Pcx2*m2*t2*t7*t8+Pcx3*m3*t4*t5*t6+Pcy2*m2*t2*t4*t7-Pcy3*m3*t4*t6*t9-Pcz3*m3*t2*t4*t7+Pcx3*m3*t2*t5*t7*t8-Pcy3*m3*t2*t7*t8*t9);
	G(3,0) = g*m3*(Pcx3*t2*t3*t5-Pcx3*t6*t8*t9-Pcy3*t2*t3*t9-Pcy3*t5*t6*t8+Pcx3*t2*t4*t7*t9+Pcy3*t2*t4*t5*t7);

	return G;
}