#ifndef FIELD_HH
#define FIELD_HH

#include <cmath>
#include <cstdio>

using std::isnan;

/* Data structure for storing the field values at a grid point. */
struct Field {
    double u, v, w;                                 // x, y, z velocity
    double s11, s12, s13, s22, s23, s33;            // Components of the stress tensor
    double chi;                                     // Effective temperature
    double cu, cv, cw;                              // Change in velocity components
    double cs11, cs12, cs13, cs22, cs23, cs33;      // Change in stress components
    double cchi;                                    // Change in effective temperature
    double xi_x, xi_y, xi_z;                        // Components of the reference map
    double cxi_x, cxi_y, cxi_z;                     // Change in the reference map
    double sbar;                                    // sbar at the point.
    double ad_Dpl;                                  // Adaptively-computed Dpl value.

    Field() {};
    Field(double a) : 
        u(a), v(a), w(a), s11(a), s12(a), s13(a), s22(a), s23(a), s33(a),
        chi(a), cu(a), cv(a), cw(a), cs11(a), cs12(a), cs13(a), cs22(a), cs23(a), cs33(a),
        cchi(a), xi_x(a), xi_y(a), xi_z(a), cxi_x(a), cxi_y(a), cxi_z(a),
        sbar(a), ad_Dpl(a) {};

    /* Updates all components at the current grid point according to the
     * (assumed to be) previously calculated changes. */
    inline void update() {
        u += cu;
        v += cv;
        w += cw;
        s11 += cs11;
        s12 += cs12;
        s13 += cs13;
        s22 += cs22;
        s23 += cs23;
        s33 += cs33;
        chi += cchi;
        xi_x += cxi_x;
        xi_y += cxi_y;
        xi_z += cxi_z;

        cu=cv=cw=cs11=cs12=cs13=cs22=cs23=cs33=cchi=cxi_x=cxi_y=cxi_z=0;
    }

    /* Subtract another field from this one. */
    inline Field operator- (Field f) {
        Field nf(0);
        nf.u   = u - f.u;
        nf.v   = v - f.v;
        nf.w   = w - f.w;
        nf.s11 = s11 - f.s11;
        nf.s12 = s12 - f.s12;
        nf.s13 = s13 - f.s13;
        nf.s22 = s22 - f.s22;
        nf.s23 = s23 - f.s23;
        nf.s33 = s33 - f.s33;
        nf.chi = chi - f.chi;
        return nf;
    }

    inline bool operator== (Field f) {
        return ((f.u == u) && (f.v == v) && (f.w == w) && (f.s11 == s11) && (f.s12 == s12) \
                && (f.s13 == s13) && (f.s22 == s22) && (f.s23 == s23) && (f.s33 == s33) && (f.chi == chi));
    }

    inline bool operator!= (Field f) {
        return !(f == *this);
    }

    /* Add another Field to this one. */
    inline Field operator+ (Field f){
        Field nf(0);
        nf.u   = u + f.u;
        nf.v   = v + f.v;
        nf.w   = w + f.w;
        nf.s11 = s11 + f.s11;
        nf.s12 = s12 + f.s12;
        nf.s13 = s13 + f.s13;
        nf.s22 = s22 + f.s22;
        nf.s23 = s23 + f.s23;
        nf.s33 = s33 + f.s33;
        nf.chi = chi + f.chi;
        return nf;
    }

    /* Scale multiply corresponding field values. */
    inline Field operator* (Field f){
        Field nf(0);
        nf.u   = u*f.u;
        nf.v   = v*f.v;
        nf.w   = w*f.w;
        nf.s11 = s11*f.s11;
        nf.s12 = s12*f.s12;
        nf.s13 = s13*f.s13;
        nf.s22 = s22*f.s22;
        nf.s23 = s23*f.s23;
        nf.s33 = s33*f.s33;
        nf.chi = chi*f.chi;
        return nf;
    }

    /* Scale all values by a double. */
    inline Field operator* (double val){
        Field nf(0);
        nf.u   = val*u;
        nf.v   = val*v;
        nf.w   = val*w;
        nf.s11 = val*s11;
        nf.s12 = val*s12;
        nf.s13 = val*s13;
        nf.s22 = val*s22;
        nf.s23 = val*s23;
        nf.s33 = val*s33;
        nf.chi = val*chi;
        return nf;
    }

    /** Diagnostic function for finding weird values. */
    inline bool weird() {
        if(isnan(u)||fabs(u)>1e10) return true;
        if(isnan(v)||fabs(v)>1e10) return true;
        if(isnan(w)||fabs(w)>1e10) return true;
        if(isnan(s11)||fabs(s11)>1e10) return true;
        if(isnan(s12)||fabs(s12)>1e10) return true;
        if(isnan(s13)||fabs(s13)>1e10) return true;
        if(isnan(s22)||fabs(s22)>1e10) return true;
        if(isnan(s23)||fabs(s23)>1e10) return true;
        if(isnan(s33)||fabs(s33)>1e10) return true;
        if(isnan(chi)||fabs(chi)>1e10) return true;
        return false;
    }

    /** Diagnostic routine for printing out the velocity and stress. */
    inline void print() {
        printf("(%g,%g,%g)   (%g,%g,%g,%g,%g,%g)  [%.15g]\n",u,v,w,s11,s12,s13,s22,s23,s33,dev());
    }

    /** Calculates the pressure, given by minus a third of the trace of the stress
     * components. */
    inline double pressure() {
        return -(s11+s22+s33)*(1/3.);
    }

    /** Calculates the square of the deviatoric stress. */
    inline double devsq() {
        double p   = pressure(),
               t11 = s11 + p,
               t22 = s22 + p,
               t33 = s33 + p;
        return 0.5*(t11*t11+t22*t22 + t33*t33) + s12*s12 + s13*s13 + s23*s23;
    }

    /** Returns a particular field member indexed by a number. */
    inline double fval(int i) {
        switch(i) {
            case 0: return u;
            case 1: return v;
            case 2: return w;
            case 3: return s11;
            case 4: return s12;
            case 5: return s13;
            case 6: return s22;
            case 7: return s23;
            case 8: return s33;
            case 9: return pressure();
            case 10: return dev();
            case 11: return chi;
            case 12: return chi*21000;
            case 13: return xi_x;
            case 14: return xi_y;
            case 15: return xi_z;
        }
        return 0;
    }

    /** Calculates the deviatoric stress. */
    inline double dev() {
        double val = sqrt(devsq());
        sbar = val;
        return val;}
};

#endif
