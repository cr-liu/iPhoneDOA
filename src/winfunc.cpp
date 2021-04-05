#include "winfunc.hpp"

#include <cmath>

WinFunc::WinFunc(const int length, const float alpha)
{
    win.resize(length);
    for(int i=0; i<length; ++i)
        win(i) = alpha - (1-alpha)*cos(2*M_PI*i/(length-1));
}

float& WinFunc::operator[](const int n)
{
    return win(n);
}
