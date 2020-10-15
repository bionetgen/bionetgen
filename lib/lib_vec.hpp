// from: https://joncraton.org/blog/67/simple-vector-library-c/

#ifndef LIB_VEC3_H_
#define LIB_VEC3_H_

#include <cstdio>
#include <vector>
#include <cmath>

class vec3 {
    public:
        vec3();
        vec3(double, double, double);
        vec3(std::vector<double>);
        bool operator==(vec3 rhs);
        vec3 operator+(vec3 rhs);
        vec3 operator-(vec3 rhs);
        vec3 operator*(double scalar);
     
        vec3 cross(vec3 rhs);
        double dot(vec3 rhs);
        double length();

        double x;
        double y;
        double z;
};

#endif /* LIB_VEC3_H_ */