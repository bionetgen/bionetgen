// from: https://joncraton.org/blog/67/simple-vector-library-c/
#include "lib_vec.hpp"

#ifndef LIB_VEC_H_
#define LIB_VEC_H_

vec3::vec3(){ 
    x = 0;
    y = 0;
    z = 0;
}

vec3::vec3(std::vector<double> vec){ 
    if(vec.size() != 3)
        throw "Error in vec3 constructor: input vector size not 3";
    x = vec[0];
    y = vec[1];
    z = vec[2];
}

vec3::vec3(double _x, double _y, double _z){ 
    x = _x;
    y = _y;
    z = _z;
}

bool vec3::operator==(vec3 rhs) {
    return(x == rhs.x && y == rhs.y && z == rhs.z);
}

vec3 vec3::operator+(vec3 rhs) {
    return vec3( x + rhs.x, 
                 y + rhs.y, 
                 z + rhs.z);
}

vec3 vec3::operator-(vec3 rhs) {
    return vec3( x - rhs.x, 
                 y - rhs.y, 
                 z - rhs.z);
}

vec3 vec3::operator*(double scalar) {
    return vec3( x * scalar, 
                 y * scalar,
                 z * scalar);
}

double vec3::dot(vec3 rhs) {
    return (x * rhs.x + 
            y * rhs.y + 
            z * rhs.z);
}

vec3 vec3::cross(vec3 rhs) {
    return vec3( y * rhs.z - z * rhs.y,
                 z * rhs.x - x * rhs.z,
                 x * rhs.y - y * rhs.x);
}

double vec3::length() {
    return double(std::sqrt( x*x + y*y + z*z ));
}

#endif /* LIB_VEC_H_ */