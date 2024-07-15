#ifndef UTILS_HH
#define UTILS_HH

#include <cmath>
#include <limits>

inline double xlog(const double& x){
    if (x == 0.0){
        return -std::numeric_limits<double>::infinity();
    }
    return log(x);
}

inline double xlogsumexp(const double& x, const double& y){
    if (x == -std::numeric_limits<double>::infinity()){
        return y;
    }
    if (y == -std::numeric_limits<double>::infinity()){
        return x;
    }
    if (x > y){
        return x + log(1 + exp(y - x));
    }
    return y + log(1 + exp(x - y));
}

inline double xlogsumexp(const double& x, const double& y, const double& z){
    return xlogsumexp(x, xlogsumexp(y, z));
}

#endif