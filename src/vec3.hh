#ifndef VEC_HH
#define VEC_HH

#include <cmath>

/* Basic Three-Component Vector Class */
struct vec3 {

	public:

        /* Three vector components */
		double x, y, z;
        
        /* Empty vector */
		vec3() {};

        /* Vector with all equal entries */
		vec3(double x_) : x(x_), y(x_), z(x_) {};

        /* Vector with distinct entries */
		vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {};

        void reset(double _x, double _y, double _z) {
            x = _x;
            y = _y;
            z = _z;
        }

        /* Set every element */
        inline vec3 operator= (double val){
            x = y = z = val;
            return *this;
        }

        /* Set every element */
        inline vec3 operator= (int val){
            x = y = z = val;
            return *this;
        }

        /* Vector addition */
		inline vec3 operator+ (vec3 p) {
			return vec3(x + p.x, y + p.y, z + p.z);
		}

        /* Vector subtraction */
		inline vec3 operator- (vec3 p) {
			return vec3(x - p.x, y - p.y, z - p.z);
		}

        /* Scalar-Vector Multiplication */
		inline vec3 operator* (double a) {
			return vec3(x*a, y*a, z*a);
		}

        /* Scalar-Vector Division */
		inline vec3 operator/ (double a) {
			return vec3(x/a, y/a, z/a);
		}

        /* In-Place Vector Addition */
		inline void operator+= (vec3 p) {
			x += p.x;
            y += p.y;
            z += p.z;
		}

        /* In-Place Vector Subtraction */
		inline void operator-= (vec3 p) {
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

        /* Element Lookup */
        const double &operator[](int n) const{
            switch (n) {
                case 0:  return x; break;
                case 1:  return y; break;
                case 2:  return z; break;
                default: return x;
            }
        }

        /* Element Lookup */
        const double &operator() (int n) const{
            switch (n) {
                case 0:  return x; break;
                case 1:  return y; break;
                case 2:  return z; break;
                default: return x;
            }
        }

        /* Cast a vector to a double by calculating the square magnitude */
		inline operator double() {
			return double(x*x + y*y + z*z);
		}

        inline operator float(){
            return float(x);
        }

        /* Calculate the l2-norm */
		inline double magnitude() {
			return sqrt(x*x + y*y + z*z);
		}

        /* Calculate the squared l2-norm */
		inline double modsq() {
			return x*x + y*y + z*z;
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
inline vec3 operator*(const double e, vec3 f) {
	return vec3(e*f.x, e*f.y, e*f.z);
}

/* Vector negation */
inline vec3 operator-(vec3 f) {
	return vec3(-f.x, -f.y, -f.z);
}

/* Vector dot product */
inline double operator*(vec3 v1, vec3 v2){
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}


/* Non-member modulus function */
inline double mod_sq(vec3 a) {
    return a.x*a.x + a.y*a.y + a.z*a.z;
}

#endif
