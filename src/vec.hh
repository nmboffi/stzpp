#ifndef VEC_HH
#define VEC_HH

#include <cmath>

/* Basic Three-Component Vector Class */
struct vec {

	public:

        /* Three vector components */
		double x, y, z;
        
        /* Empty vector */
		vec() {};

        /* Vector with all equal entries */
		vec(double x_) : x(x_), y(x_), z(x_) {};

        /* Vector with distinct entries */
		vec(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {};

        /* Set every element */
        inline vec operator= (double val){
            x = y = z = val;
            return *this;
        }

        /* Vector addition */
		inline vec operator+ (vec p) {
			return vec(x + p.x, y + p.y, z + p.z);
		}

        /* Vector subtraction */
		inline vec operator- (vec p) {
			return vec(x - p.x, y - p.y, z - p.z);
		}

        /* Scalar-Vector Multiplication */
		inline vec operator* (double a) {
			return vec(x*a, y*a, z*a);
		}

        /* Scalar-Vector Division */
		inline vec operator/ (double a) {
			return vec(x/a, y/a, z/a);
		}

        /* In-Place Vector Addition */
		inline void operator+= (vec p) {
			x += p.x;
            y += p.y;
            z += p.z;
		}

        /* In-Place Vector Subtraction */
		inline void operator-= (vec p) {
			x -= p.x;
            y -= p.y;
            z -= p.z;
		}

        /* In-Place Scalar-Vector Multiplication */
		inline void operator*= (double a) {
			x *= a;
            y *= a;
            z *= a;
		}

        /* In-Place Scalar-Vector Division */
		inline void operator/= (double a) {
			x /= a;
            y /= a;
            z /= a;
		}

        /* Cast a vector to a double by calculating the square magnitude */
		inline operator double() {
			return double(x*x + y*y + z*z);
		}

        /* Calculate the l2-norm */
		inline double magnitude() {
			return sqrt(x*x + y*y + z*z);
		}

        /* Normalize the vector */
		inline void normalize() {
            // Compute the magnitude
            double mag = sqrt(x*x + y*y + z*z);

            // Scale the components
			if(mag > 0) {
				x /= mag;
				y /= mag;
                z /= mag;
			}
		}
};

/* Scalar-first scalar-vector multiplication */
inline vec operator*(const double e, vec f) {
	return vec(e*f.x, e*f.y, e*f.z);
}

/* Vector negation */
inline vec operator-(vec f) {
	return vec(-f.x, -f.y, -f.z);
}

/* Non-member modulus function */
inline double mod_sq(vec a) {
    return a.x*a.x + a.y*a.y + a.z*a.z;
}

#endif
