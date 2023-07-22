#ifndef STENC_FUNCS_HH
#define STENC_FUNCS_HH

inline double Uxp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t11*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (lamb*(-(t23*t23) + t22*t33))/(D*(dx*dx)) + 
   (mu*(-(t23*t23) + t22*t33))/(D*(dx*dx)) + 
   (mu*t11*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uc_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*t11*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) - 
   (2*mu*t11*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t11*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t11*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t11*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t11*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t11*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t11*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*lamb*(-(t23*t23) + t22*t33))/(D*(dx*dx)) - 
   (2*mu*(-(t23*t23) + t22*t33))/(D*(dx*dx)) - 
   (2*mu*t11*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uxm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t11*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (lamb*(-(t23*t23) + t22*t33))/(D*(dx*dx)) + 
   (mu*(-(t23*t23) + t22*t33))/(D*(dx*dx)) + 
   (mu*t11*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t11*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t11*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Uym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t11*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t11*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Uzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t11*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t11*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Uzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t11*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t11*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Uxpyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*t11*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t11*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxmym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*t11*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t11*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxpym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t11*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*t11*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t11*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxmyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t11*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*t11*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t11*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxpzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*t11*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t11*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t11*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxmzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*t11*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t11*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t11*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxpzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*t11*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t11*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t11*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxmzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*t11*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t11*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t11*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uypzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t11*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t11*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uymzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t11*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t11*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t11*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uypzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t11*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t11*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t11*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uymzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t11*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t11*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t11*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vxp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t12*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vc_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*t12*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) - 
   (2*mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t12*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*lamb*(t13*t23 - t12*t33))/(D*(dy*dy)) - 
   (2*mu*(t13*t23 - t12*t33))/(D*(dy*dy)) - 
   (2*mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t12*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vxm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t12*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t12*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (lamb*(t13*t23 - t12*t33))/(D*(dy*dy)) + 
   (mu*(t13*t23 - t12*t33))/(D*(dy*dy)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Vym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t12*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (lamb*(t13*t23 - t12*t33))/(D*(dy*dy)) + 
   (mu*(t13*t23 - t12*t33))/(D*(dy*dy)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Vzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Vzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Vxpyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) + 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) + 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxmym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) + 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) + 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxpym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) - 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) - 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxmyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) - 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dy) - 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxpzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxmzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxpzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxmzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vypzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (lamb*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vymzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (lamb*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vypzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (lamb*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vymzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (lamb*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dy*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wxp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t13*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wc_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*t13*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) - 
   (2*mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*lamb*(-(t13*t22) + t12*t23))/(D*(dz*dz)) - 
   (2*mu*(-(t13*t22) + t12*t23))/(D*(dz*dz)) - 
   (2*mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t13*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t13*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wxm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t13*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t13*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Wym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t13*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Wzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (lamb*(-(t13*t22) + t12*t23))/(D*(dz*dz)) + 
   (mu*(-(t13*t22) + t12*t23))/(D*(dz*dz)) + 
   (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Wzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (lamb*(-(t13*t22) + t12*t23))/(D*(dz*dz)) + 
   (mu*(-(t13*t22) + t12*t23))/(D*(dz*dz)) + 
   (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Wxpyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxmym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxpym_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxmyp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxpzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) + 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxmzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) + 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxpzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) - 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxmzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (lamb*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) - 
   (mu*(-(t23*t23) + t22*t33))/(4.*D*dx*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wypzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dy*dz) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dy*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wymzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dy*dz) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dy*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wypzm_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dy*dz) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dy*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wymzp_1(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dy*dz) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dy*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uxp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (lamb*(t13*t23 - t12*t33))/(D*(dx*dx)) + 
   (mu*(t13*t23 - t12*t33))/(D*(dx*dx)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t12*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uc_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*t12*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) - 
   (2*mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t12*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*lamb*(t13*t23 - t12*t33))/(D*(dx*dx)) - 
   (2*mu*(t13*t23 - t12*t33))/(D*(dx*dx)) - 
   (2*mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t12*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uxm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (lamb*(t13*t23 - t12*t33))/(D*(dx*dx)) + 
   (mu*(t13*t23 - t12*t33))/(D*(dx*dx)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t12*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t12*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Uym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t12*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t12*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Uzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Uzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t12*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t12*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Uxpyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) + 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) + 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxmym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) + 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) + 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxpym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) - 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) - 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxmyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) - 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dx*dy) - 
   (mu*t12*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t12*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxpzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dx*dz) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dz) + 
   (mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxmzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dx*dz) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dz) + 
   (mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxpzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dx*dz) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dz) - 
   (mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxmzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dx*dz) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dz) - 
   (mu*t12*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uypzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uymzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uypzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uymzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t12*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t12*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t12*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vxp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t22*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t22*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vc_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*t22*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) - 
   (2*mu*t22*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t22*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t22*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t22*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*lamb*(-(t13*t13) + t11*t33))/(D*(dy*dy)) - 
   (2*mu*(-(t13*t13) + t11*t33))/(D*(dy*dy)) - 
   (2*mu*t22*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t22*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t22*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t22*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vxm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t22*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t22*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (lamb*(-(t13*t13) + t11*t33))/(D*(dy*dy)) + 
   (mu*(-(t13*t13) + t11*t33))/(D*(dy*dy)) + 
   (mu*t22*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t22*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Vym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (lamb*(-(t13*t13) + t11*t33))/(D*(dy*dy)) + 
   (mu*(-(t13*t13) + t11*t33))/(D*(dy*dy)) + 
   (mu*t22*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t22*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Vzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t22*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t22*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Vzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t22*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t22*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Vxpyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*t22*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t22*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxmym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) + 
   (mu*t22*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t22*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxpym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t22*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*t22*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t22*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxmyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t22*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dy) - 
   (mu*t22*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t22*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxpzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t22*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t22*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxmzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t22*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t22*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t22*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxpzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t22*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t22*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t22*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxmzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t22*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t22*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t22*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vypzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*t22*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t22*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t22*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vymzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*t22*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t22*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t22*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vypzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*t22*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t22*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t22*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vymzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*t22*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t22*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t22*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wxp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t23*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wc_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t23)/(D*D*(dz*dz)) - 
   (2*lamb*(t12*t13 - t11*t23))/(D*(dz*dz)) - 
   (2*mu*(t12*t13 - t11*t23))/(D*(dz*dz)) - 
   (2*mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t23*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t23*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wxm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t23*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t23*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Wym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t23*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Wzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t23)/(D*D*(dz*dz)) + 
   (lamb*(t12*t13 - t11*t23))/(D*(dz*dz)) + 
   (mu*(t12*t13 - t11*t23))/(D*(dz*dz)) + 
   (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Wzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t23)/(D*D*(dz*dz)) + 
   (lamb*(t12*t13 - t11*t23))/(D*(dz*dz)) + 
   (mu*(t12*t13 - t11*t23))/(D*(dz*dz)) + 
   (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Wxpyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxmym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxpym_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxmyp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxpzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dz) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxmzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dz) + 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxpzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dz) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxmzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (lamb*(t13*t23 - t12*t33))/(4.*D*dx*dz) - 
   (mu*(t13*t23 - t12*t33))/(4.*D*dx*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wypzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) + 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wymzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) + 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wypzm_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) - 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wymzp_2(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (lamb*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) - 
   (mu*(-(t13*t13) + t11*t33))/(4.*D*dy*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uxp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(D*(dx*dx)) + 
   (mu*(-(t13*t22) + t12*t23))/(D*(dx*dx)) + 
   (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t13*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uc_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*t13*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) - 
   (2*mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*lamb*(-(t13*t22) + t12*t23))/(D*(dx*dx)) - 
   (2*mu*(-(t13*t22) + t12*t23))/(D*(dx*dx)) - 
   (2*mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t13*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t13*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uxm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(D*(dx*dx)) + 
   (mu*(-(t13*t22) + t12*t23))/(D*(dx*dx)) + 
   (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t13*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Uyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t13*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Uym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t13*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t13*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Uzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Uzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22)))/(D*D*(dz*dz)) + 
   (mu*t13*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t13*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Uxpyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dx*dy) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dy) + 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxmym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dx*dy) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dy) + 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxpym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dx*dy) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dy) - 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxmyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dx*dy) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dx*dy) - 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t13*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Uxpzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) + 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) + 
   (mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxmzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) + 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) + 
   (mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxpzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) - 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) - 
   (mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uxmzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) - 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dx*dz) - 
   (mu*t13*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Uypzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uymzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uypzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Uymzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*t13*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t13*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t13*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vxp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t23*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vc_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t23)/(D*D*(dz*dz)) - 
   (2*lamb*(t12*t13 - t11*t23))/(D*(dy*dy)) - 
   (2*mu*(t12*t13 - t11*t23))/(D*(dy*dy)) - 
   (2*mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) - 
   (2*mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) - 
   (2*mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz)) - 
   (2*mu*t23*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t23*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vxm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dx*dx)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t23*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Vyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(D*(dy*dy)) + 
   (mu*(t12*t13 - t11*t23))/(D*(dy*dy)) + 
   (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t23*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Vym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(D*(dy*dy)) + 
   (mu*(t12*t13 - t11*t23))/(D*(dy*dy)) + 
   (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dy*dy)) + 
   (mu*t23*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t23*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Vzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t23)/(D*D*(dz*dz)) + 
   (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Vzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t23)/(D*D*(dz*dz)) + 
   (mu*t23*((t12*t13 - t11*t23)*(t12*t13 - t11*t23)))/(D*D*(dz*dz)) + 
   (mu*t23*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23)))/(D*D*(dz*dz));

}

inline double Vxpyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) + 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxmym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) + 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxpym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) - 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxmyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dy) - 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t23*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Vxpzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxmzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxpzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vxmzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(-(t12*t12) + t11*t22)*t23*(-(t13*t22) + t12*t23))/(2.*(D*D)*dx*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Vypzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) + 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) + 
   (mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vymzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) + 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) + 
   (mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) + 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vypzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) - 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) - 
   (mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Vymzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) - 
   (mu*(-(t12*t12) + t11*t22))/(4.*D*dy*dz) - 
   (mu*(-(t12*t12) + t11*t22)*t23*(t12*t13 - t11*t23))/(2.*(D*D)*dy*dz) - 
   (mu*t23*(t12*t13 - t11*t23)*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*t23*(-(t13*t22) + t12*t23)*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wxp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23))*t33)/(D*D*(dx*dx)) + 
   (mu*t33*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t33*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wc_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (-2*lamb*(-(t12*t12) + t11*t22))/(D*(dz*dz)) - 
   (2*mu*(-(t12*t12) + t11*t22))/(D*(dz*dz)) - 
   (2*mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t33)/(D*D*(dz*dz)) - 
   (2*mu*((t12*t13 - t11*t23)*(t12*t13 - t11*t23))*t33)/(D*D*(dy*dy)) - 
   (2*mu*((t12*t13 - t11*t23)*(t12*t13 - t11*t23))*t33)/(D*D*(dz*dz)) - 
   (2*mu*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23))*t33)/(D*D*(dx*dx)) - 
   (2*mu*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23))*t33)/(D*D*(dz*dz)) - 
   (2*mu*t33*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t33*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) - 
   (2*mu*t33*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy)) - 
   (2*mu*t33*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wxm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23))*t33)/(D*D*(dx*dx)) + 
   (mu*t33*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dx*dx)) + 
   (mu*t33*((-(t23*t23) + t22*t33)*(-(t23*t23) + t22*t33)))/(D*D*(dx*dx));

}

inline double Wyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((t12*t13 - t11*t23)*(t12*t13 - t11*t23))*t33)/(D*D*(dy*dy)) + 
   (mu*t33*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t33*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Wym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*((t12*t13 - t11*t23)*(t12*t13 - t11*t23))*t33)/(D*D*(dy*dy)) + 
   (mu*t33*((-(t13*t13) + t11*t33)*(-(t13*t13) + t11*t33)))/(D*D*(dy*dy)) + 
   (mu*t33*((t13*t23 - t12*t33)*(t13*t23 - t12*t33)))/(D*D*(dy*dy));

}

inline double Wzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t12*t12) + t11*t22))/(D*(dz*dz)) + 
   (mu*(-(t12*t12) + t11*t22))/(D*(dz*dz)) + 
   (mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t33)/(D*D*(dz*dz)) + 
   (mu*((t12*t13 - t11*t23)*(t12*t13 - t11*t23))*t33)/(D*D*(dz*dz)) + 
   (mu*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23))*t33)/(D*D*(dz*dz));

}

inline double Wzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t12*t12) + t11*t22))/(D*(dz*dz)) + 
   (mu*(-(t12*t12) + t11*t22))/(D*(dz*dz)) + 
   (mu*((-(t12*t12) + t11*t22)*(-(t12*t12) + t11*t22))*t33)/(D*D*(dz*dz)) + 
   (mu*((t12*t13 - t11*t23)*(t12*t13 - t11*t23))*t33)/(D*D*(dz*dz)) + 
   (mu*((-(t13*t22) + t12*t23)*(-(t13*t22) + t12*t23))*t33)/(D*D*(dz*dz));

}

inline double Wxpyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dy) + 
   (mu*t33*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t33*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxmym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (mu*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dy) + 
   (mu*t33*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) + 
   (mu*t33*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxpym_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dy) - 
   (mu*t33*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t33*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxmyp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(mu*(t12*t13 - t11*t23)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dy) - 
   (mu*t33*(-(t13*t13) + t11*t33)*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dy) - 
   (mu*t33*(t13*t23 - t12*t33)*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dy);

}

inline double Wxpzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dz) + 
   (mu*(t12*t13 - t11*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*(-(t13*t22) + t12*t23)*t33*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxmzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) + 
   (mu*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dz) + 
   (mu*(t12*t13 - t11*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) + 
   (mu*(-(t13*t22) + t12*t23)*t33*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxpzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dz) - 
   (mu*(t12*t13 - t11*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*(-(t13*t22) + t12*t23)*t33*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wxmzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*(-(t13*t22) + t12*t23))/(4.*D*dx*dz) - 
   (mu*(-(t12*t12) + t11*t22)*(-(t13*t22) + t12*t23)*t33)/(2.*(D*D)*dx*dz) - 
   (mu*(t12*t13 - t11*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dx*dz) - 
   (mu*(-(t13*t22) + t12*t23)*t33*(-(t23*t23) + t22*t33))/(2.*(D*D)*dx*dz);

}

inline double Wypzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23)*t33)/(2.*(D*D)*dy*dz) + 
   (mu*(t12*t13 - t11*t23)*t33*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*(-(t13*t22) + t12*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wymzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return (lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) + 
   (mu*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23)*t33)/(2.*(D*D)*dy*dz) + 
   (mu*(t12*t13 - t11*t23)*t33*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) + 
   (mu*(-(t13*t22) + t12*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wypzm_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23)*t33)/(2.*(D*D)*dy*dz) - 
   (mu*(t12*t13 - t11*t23)*t33*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*(-(t13*t22) + t12*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

inline double Wymzp_3(double t11, double t12, double t13,
			double t21, double t22, double t23,
			double t31, double t32, double t33,
			double mu, double lamb, double D,
			double dx, double dy, double dz) {
	
 return -(lamb*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*(t12*t13 - t11*t23))/(4.*D*dy*dz) - 
   (mu*(-(t12*t12) + t11*t22)*(t12*t13 - t11*t23)*t33)/(2.*(D*D)*dy*dz) - 
   (mu*(t12*t13 - t11*t23)*t33*(-(t13*t13) + t11*t33))/(2.*(D*D)*dy*dz) - 
   (mu*(-(t13*t22) + t12*t23)*t33*(t13*t23 - t12*t33))/(2.*(D*D)*dy*dz);

}

#endif
