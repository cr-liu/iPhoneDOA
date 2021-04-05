#ifndef CORRMAT_HPP
#define CORRMAT_HPP

#include <vector>
#include <fftw3.h>
#include <eigen3/Eigen/Dense>

class CorrMat
{
public:
    CorrMat(const int nft, fftwf_plan p = nullptr);
    ~CorrMat();
    void calc_corr(const Eigen::MatrixXf x);

    bool upperTriangular2Vec = false;
    float* in=nullptr;
    fftwf_complex *out = nullptr;
    fftwf_plan ftPlanR2C = nullptr;
    int nFFT;
    int nCh;
    int sndLen;
    Eigen::ArrayXf win;
    std::vector<Eigen::MatrixXcf> autocorrX; // freq x nCh x nCh
    std::vector<Eigen::VectorXcf> crossCorrX; // freq x C(nCh, 2)
    std::vector<Eigen::VectorXf> crossCorrXAng; // freq x C(nCh, 2)
    bool planCreated = false;
};

#endif // CORRMAT_HPP
