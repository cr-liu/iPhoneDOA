#ifndef MUSICLOC_HPP
#define MUSICLOC_HPP

#include <vector>
#include <eigen3/Eigen/Dense>

class MusicLoc
{
public:
    MusicLoc(const Eigen::MatrixXcf &R,
             const Eigen::MatrixXcf &A, int nSignal=1);
    int nSensor;
    int nDir;
    std::vector<float> spectrum;
    Eigen::VectorXf eigenValues;
    Eigen::MatrixXcf eigenVectors;
};

#endif // MUSICLOC_HPP
