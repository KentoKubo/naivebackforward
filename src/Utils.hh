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

#endif