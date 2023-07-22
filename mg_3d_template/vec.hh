#ifndef VEC_HH
#define VEC_HH

#include <cmath>

#define pi4 0.785398163397448309615660845820
#define pi8 0.392699081698724154807830422910

struct vec {
	public:
		double x, y;
		vec() {};
		vec(int x_) : x(x_), y(x_) {};
		vec(double x_) : x(x_), y(x_) {};
		vec(double x_, double y_) : x(x_), y(y_) {};
		inline vec operator+ (vec p) {
			return vec(x + p.x, y + p.y);
		}
		inline vec operator- (vec p) {
			return vec(x - p.x, y - p.y);
		}
        inline vec operator= (double val){
            x = y = val;
            return *this;
        }
        inline vec operator= (int val){
            x = y = val;
            return *this;
        }
		inline vec operator* (double a) {
			return vec(x*a, y*a);
		}
        /*
         *inline double operator* (vec v){
         *    return x*v.x + y*v.y;
         *}
         */
		inline vec operator/ (double a) {
			return vec(x/a, y/a);
		}
		inline void operator+= (vec p) {
			x += p.x;
            y += p.y;
		}
		inline void operator-= (vec p) {
			x -= p.x;
            y -= p.y;
		}
		inline void operator*= (double a) {
			x *= a;
            y *= a;
		}
		inline void operator/= (double a) {
			x /= a;
            y /= a;
		}
		inline operator double() {
			return double(x*x + y*y);
		}
        inline operator float(){
            return float(x);
        }
		inline double magnitude() {
			return sqrt(x*x + y*y);
		}
		inline void normalize() {
			if(this -> magnitude() > 0) {
				x /= this->magnitude();
				y /= this->magnitude();
			}
		}
		inline int octant() {
			double angle = atan2(x,y);
			const double pi = 3.1415926535897932384626433832795;
			if(angle < 0)
				angle = 2 * pi + angle;
			return int(angle / pi4);
		}
};

inline vec operator*(const double e, vec f) {
	return vec(e*f.x,e*f.y);
}

inline vec operator-(vec f) {
	return vec(-f.x,-f.y);
}

inline double operator*(vec v1, vec v2){
    return v1.x*v2.x + v1.y*v2.y;
}

inline double mod_sq(vec a) {return a.x*a.x+a.y*a.y;}

#endif
