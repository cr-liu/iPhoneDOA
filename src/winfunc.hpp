#ifndef WINFUNC_HPP
#define WINFUNC_HPP

#include <eigen3/Eigen/Dense>

class WinFunc
{
public:
    WinFunc(const int length, const float alpha=0.54);
    Eigen::ArrayXf win;
    float& operator[](const int n);
};

#endif // WINFUNC_HPP
