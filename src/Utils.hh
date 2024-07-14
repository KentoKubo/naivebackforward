#ifndef UTILS_HH
#define UTILS_HH

#include <cmath>
#include <limits>

inline float xlog(const float& x){
    if (x == 0.0){
        return -std::numeric_limits<float>::infinity();
    }
    return log(x);
}

inline float xlogsumexp(const float& x, const float& y){
    if (x == -std::numeric_limits<float>::infinity()){
        return y;
    }
    if (y == -std::numeric_limits<float>::infinity()){
        return x;
    }
    if (x > y){
        return x + log(1 + exp(y - x));
    }
    return y + log(1 + exp(x - y));
}

inline float xlogsumexp(const float& x, const float& y, const float& z){
    return xlogsumexp(x, xlogsumexp(y, z));
}

#endif