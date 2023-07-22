#ifndef MAT_HH
#define MAT_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "vec.hh"

/* Forward declaration for a symmetric 3x3 matrix */
class sym_mat;

/* Basic class for a 3x3 matrix */
class mat {

	public:

        /* 3x3 matrix elements */
		double a11, a12, a13, a21, a22, a23, a31, a32, a33;

        /* Basic constructor */
		mat() {};

        /* Elementwise constructor */
		mat(double vals[9]) :  a11(vals[0]), a12(vals[1]), a13(vals[2]), 
                                a21(vals[3]), a22(vals[4]), a23(vals[5]),
                                a31(vals[6]), a32(vals[7]), a33(vals[8]) {};

        /* Set every element */
        inline mat operator= (double val){
            a11 = a12 = a13 = a21 = a22 = a23 = a31 = a32 = a33 = val;
            return *this;
        }

        /* Matrix Addition */
		inline mat operator+ (mat p) {
            // Array to hold the matrix entries
            double new_vals[9] = {a11 + p.a11, a12 + p.a12, a13 + p.a13, 
                                a21 + p.a21, a22 + p.a22, a23 + p.a23, 
                                a31 + p.a31, a32 + p.a32, a33 + p.a33};

            return mat(new_vals);
        }

        /* Matrix Subtraction */
		inline mat operator- (mat p) {
            double new_vals[] = {a11 - p.a11, a12 - p.a12, a13 - p.a13, 
                               a21 - p.a21, a22 - p.a22, a23 - p.a23, 
                               a31 - p.a31, a32 - p.a32, a33 - p.a33};

            return mat(new_vals); } 
        /* Scalar-Matrix Multiplication */
		inline mat operator* (double e) {
            double new_vals[] = {e*a11, e*a12, e*a13, 
                                e*a21, e*a22, e*a23, 
                                e*a31, e*a32, e*a33};
            return mat(new_vals);
        }
        /* Matrix-Matrix Multiplication*/
		inline mat operator*(mat e) {
            double prod11 = a11*e.a11 + a12*e.a21 + a13*e.a31;
            double prod12 = a11*e.a12 + a12*e.a22 + a13*e.a32;
            double prod13 = a11*e.a13 + a12*e.a23 + a13*e.a33;

            double prod21 = a21*e.a11 + a22*e.a21 + a23*e.a31;
            double prod22 = a21*e.a12 + a22*e.a22 + a23*e.a32;
            double prod23 = a21*e.a13 + a22*e.a23 + a23*e.a33;

            double prod31 = a31*e.a11 + a32*e.a21 + a33*e.a31;
            double prod32 = a31*e.a12 + a32*e.a22 + a33*e.a32;
            double prod33 = a31*e.a13 + a32*e.a23 + a33*e.a33;

            double new_vals[] = {prod11, prod12, prod13,
                                 prod21, prod22, prod23,
                                 prod31, prod32, prod33};

            return mat(new_vals);
        }

        /* Element Lookup */
        const double &operator[](int n) const{
            switch (n){
            case 0: return a11; break;
            case 1: return a12; break;
            case 2: return a13; break;
            case 3: return a21; break;
            case 4: return a22; break;
            case 5: return a23; break;
            case 6: return a31; break;
            case 7: return a32; break;
            case 8: return a33; break;
            default: return 0;
            }
        }

        /* Matrix-Vector Multiplication */
		inline vec operator* (vec e) {
            return vec(a11*e.x + a12*e.y + a13*e.z, a21*e.x + a22*e.y + a23*e.z, a31*e.x + a32*e.y + a33*e.z);
        }

        /* Matrix-Scalar Division */
		inline mat operator/ (double e) {
			double ei=1./e;
            double new_vals[] = {a11*ei, a12*ei, a13*ei, 
                                 a21*ei, a22*ei, a23*ei, 
                                 a31*ei, a32*ei, a33*ei};

			return mat(new_vals);
		}

        /* In-Place Matrix-Matrix Addition */
		inline void operator+= (mat p) {
			a11 += p.a11;
            a12 += p.a12;
            a13 += p.a13;
            a21 += p.a21;
            a22 += p.a22;
            a23 += p.a23;
            a31 += p.a31;
            a32 += p.a32;
            a33 += p.a33;
        }

        /* In-Place Matrix-Matrix Subtraction */
		inline void operator-= (mat p) {
			a11 -= p.a11;
            a12 -= p.a12;
            a13 -= p.a13;
            a21 -= p.a21;
            a22 -= p.a22;
            a23 -= p.a23;
            a31 -= p.a31;
            a32 -= p.a32;
            a33 -= p.a33;
		}

        /* In-Place Scalar-Matrix Multiplication */
		inline void operator*= (double e) {
			a11 *= e;
            a12 *= e;
            a13 *= e;
            a21 *= e;
            a22 *= e;
            a23 *= e;
            a31 *= e;
            a32 *= e;
            a33 *= e;
		}

        /* In-Place Scalar-Matrix Division */
		inline void operator/= (double e) {
			a11 /= e;
            a12 /= e;
            a13 /= e;
            a21 /= e;
            a22 /= e;
            a23 /= e;
            a31 /= e;
            a32 /= e;
            a33 /= e;
		}

        /* Change the matrix elements */
		inline void set(double vals[9]){
            a11 = vals[0]; a12 = vals[1]; a13 = vals[2];
            a21 = vals[3]; a22 = vals[4]; a23 = vals[5];
            a31 = vals[6]; a32 = vals[7]; a33 = vals[8];
		}

        /* Compute the trace */
		inline double trace() {return a11 + a22 + a33;}

        /* Compute the determinant */
		inline double det() {
            return -a13*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33;
        }

        /* Frobenius Norm */
		inline double mod_sq() {
            return a11*a11 + a12*a12 + a13*a13 
                + a21*a21 + a22*a22 + a23*a23 
                + a31*a31 + a32*a32 + a33*a33;
        }
        /* Matrix Transpose */
		inline mat transpose() {
            double new_vals[] = {a11, a21, a31,
                                 a12, a22, a32,
                                 a13, a23, a33};
            return mat(new_vals);
        }

        /* Matrix Inverse */        
		inline mat inverse() {
            // Inverse determinant
			double idet = 1.0/det();

            double new_vals[] = {(a22*a33 - a23*a32)*idet, (a13*a32 - a12*a33)*idet, (a12*a23 - a13*a22)*idet,
                                 (a23*a31 - a21*a33)*idet, (a11*a33 - a13*a31)*idet, (a13*a21 - a11*a23)*idet,
                                 (a21*a32 - a22*a31)*idet, (a12*a31 - a11*a32)*idet, (a11*a22 - a12*a21)*idet};

            // Analytical form for 3x3 matrices
			return mat(new_vals);
		}

        /* Matrix Inverse Transpose */
		inline mat inv_transpose() {
            // Inverse determinant 
			double idet = 1.0/det();

            double new_vals[] = {(a22*a33 - a23*a32)*idet, (a23*a31 - a21*a33)*idet, (a21*a32 - a22*a31)*idet,
                                 (a13*a32 - a12*a33)*idet, (a11*a33 - a13*a31)*idet, (a12*a31 - a11*a32)*idet,
                                 (a12*a23 - a13*a22)*idet, (a13*a21 - a11*a23)*idet, (a11*a22 - a12*a21)*idet};

            // Just the transpose of the above
			return mat(new_vals);
		}
};

/* Scalar-First Scalar-Matrix Multiplication */
inline mat operator*(const double e, mat f) {
    double new_vals[] = {e*f.a11, e*f.a12, e*f.a13, 
                         e*f.a21, e*f.a22, e*f.a23, 
                         e*f.a31, e*f.a32, e*f.a33};
    return mat(new_vals);
}

/* Scalar-MatInv Multiplication */
inline mat operator/(double e,mat f) {
	return e*f.inverse();
}

/* Matrix Negation */
inline mat operator-(mat f) {
    double new_vals[] = {-f.a11, -f.a12, -f.a13,
                         -f.a21, -f.a22, -f.a23,
                         -f.a31, -f.a32, -f.a33};
	return mat(new_vals);
}

/* Symmetric matrices - complete later */
/*
 *struct sym_mat {
 *    public:
 *        double a,b,d;
 *        sym_mat() {};
 *        sym_mat(double a_) : a(a_),b(0),d(a_) {};
 *        sym_mat(double a_,double b_,double d_) : a(a_),b(b_),d(d_) {};
 *        inline sym_mat operator+ (sym_mat p) {return sym_mat(a+p.a,b+p.b,d+p.d);}
 *        inline sym_mat operator- (sym_mat p) {return sym_mat(a-p.a,b-p.b,d-p.d);}
 *        inline sym_mat operator* (double e) {return sym_mat(a*e,b*e,d*e);}
 *        inline sym_mat operator*(sym_mat e) {return sym_mat(a*e.a+b*e.b,a*e.b+b*e.d,b*e.b+d*e.d);}
 *        inline vec operator* (vec e) {return vec(a*e.x+b*e.y,b*e.x+d*e.y);}
 *        inline sym_mat operator/ (double e) {
 *            double ei=1/e;
 *            return sym_mat(a*ei,b*ei,d*ei);
 *        }
 *        inline void operator+= (sym_mat p) {a+=p.a;b+=p.b;d+=p.d;}
 *        inline void operator-= (sym_mat p) {a-=p.a;b-=p.b;d-=p.d;}
 *        inline void operator*= (double e) {a*=e;b*=e;d*=e;}
 *        inline void operator/= (double e) {a/=e;b/=e;d/=e;}
 *        inline void set(double a_,double b_,double d_) {a=a_;b=b_;d=d_;}
 *        inline double devmod() {
 *            return sqrt(0.5*(a-d)*(a-d)+2*b*b);
 *        }
 *        inline double trace() {return a+d;}
 *        inline double det() {return a*d-b*b;}
 *        inline double mod_sq() {return a*a+2*b*b+d*d;}
 *        inline sym_mat transpose() {return sym_mat(a,b,d);}
 *        inline sym_mat inverse() {
 *            double idet=1.0/det();
 *            return sym_mat(idet*d,-idet*b,idet*a);
 *        }
 *        inline sym_mat inv_transpose() {
 *            return inverse();
 *        }
 *        inline void print_sym_mat() {printf(" [%g %g %g %g]",a,b,b,d);}
 *        void eigenvectors(double &l1,double &l2,mat &Lam);
 *};
 *
 *inline sym_mat mat::AAT() {
 *    return sym_mat(a*a+b*b,a*c+b*d,c*c+d*d);
 *}
 *
 *inline sym_mat mat::ATA() {
 *    return sym_mat(a*a+c*c,a*b+c*d,b*b+d*d);
 *}
 *
 *inline sym_mat mat::ATDA(double l1,double l2) {
 *    return sym_mat(a*a*l1+c*c*l2,a*b*l1+c*d*l2,b*b*l1+d*d*l2);
 *}
 *
 *inline sym_mat mat::ATSA(sym_mat s) {
 *    double e=a*s.a+c*s.b,f=a*s.b+c*s.d,
 *           g=b*s.a+d*s.b,h=b*s.b+d*s.d;
 *    return sym_mat(a*e+c*f,b*e+d*f,b*g+d*h);
 *}
 *
 *inline sym_mat operator*(const double e,sym_mat f) {
 *    return sym_mat(e*f.a,e*f.b,e*f.d);
 *}
 *
 *inline sym_mat operator/(double e,sym_mat f) {
 *    double idet=e/(f.a*f.d-f.b*f.b);
 *    return sym_mat(f.d*idet,-f.b*idet,f.a*idet);
 *}
 *
 *inline sym_mat operator-(sym_mat f) {
 *    return sym_mat(-f.a,-f.b,-f.d);
 *}
 */
#endif
