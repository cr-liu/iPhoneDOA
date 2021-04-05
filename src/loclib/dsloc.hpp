#ifndef DSLOC_HPP
#define DSLOC_HPP

#include <vector>
#include <eigen3/Eigen/Dense>

class DSLoc
{
public:
    DSLoc(const Eigen::MatrixXcf &R, const Eigen::MatrixXcf &A);
    int nSensor;
    int nDir;
    std::vector<float> spectrum;

    std::vector<float> spectrumCosCart;
    std::vector<float> spectrumEucCart;
    std::vector<float> spectrumCosPolar;
    std::vector<float> spectrumEucPolar;
    void calc_cos_sim_cartesian(const Eigen::VectorXcf& crossR,
                      const Eigen::MatrixXcf& crossA);
    void calc_euc_dist_cartesian(const Eigen::VectorXcf& crossR,
                      const Eigen::MatrixXcf& crossA);
    void calc_cos_sim_polar(const Eigen::VectorXf& crossR,
                            const Eigen::MatrixXf& crossA);
    void calc_euc_dist_polar(const Eigen::VectorXf& crossR,
                             const Eigen::MatrixXf& crossA);
};

#endif // DSLOC_HPP
