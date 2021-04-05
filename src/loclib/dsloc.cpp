#include "dsloc.hpp"

#include <iostream>

using namespace Eigen;

DSLoc::DSLoc(const Eigen::MatrixXcf &R, const Eigen::MatrixXcf &A) :
    nSensor(A.cols()),
    nDir(A.rows())
{
    // DS Spatial Spectrum
    spectrum.resize(nDir, 0.0);
    for(int i = 0; i < nDir; ++i) {
/*
        VectorXcd tmp;
        tmp.noalias() = R * a.row(i).transpose();
        for(m=0; m<numberOfSensor; ++m){
            spectrum[i] += (conj(a(i,m))*tmp[m]).real();
        }
*/
        std::complex<float> tmp;
        tmp = (A.row(i) * R).dot(A.row(i));
        spectrum[i] = tmp.real() / (float)nSensor;
    }
}

void DSLoc::calc_cos_sim_cartesian(const Eigen::VectorXcf &crossR,
                                   const Eigen::MatrixXcf &crossA) {
    assert(crossR.size() == crossA.cols());

    spectrumCosCart.resize(nDir, 0.0);
    for(int i = 0; i < nDir; ++i) {
        spectrumCosCart[i] = crossR.conjugate().dot(crossA.row(i)).real();
    }
}

void DSLoc::calc_euc_dist_cartesian(const Eigen::VectorXcf &crossR,
                         const Eigen::MatrixXcf &crossA) {
    assert(crossR.size() == crossA.cols());

    spectrumEucCart.resize(nDir, 0.0);
    for(int i = 0; i < nDir; ++i) {
        for(int j = 0; j < crossA.cols(); ++j) {
            spectrumEucCart[i] += (crossR[j].real() - crossA(i, j).real()) *
                    (crossR[j].real() - crossA(i, j).real());
            spectrumEucCart[i] += (crossR[j].imag() + crossA(i, j).imag()) *
                    (crossR[j].imag() + crossA(i, j).imag());
        }
    }
}

void DSLoc::calc_cos_sim_polar(const Eigen::VectorXf &crossR,
                               const Eigen::MatrixXf &crossA) {
    assert(crossR.size() == crossA.cols());

    spectrumCosPolar.resize(nDir, 0.0);
    for(int i = 0; i < nDir; ++i) {
        spectrumCosPolar[i] = crossR.dot(crossA.row(i));
    }
}

void DSLoc::calc_euc_dist_polar(const Eigen::VectorXf &crossR,
                                const Eigen::MatrixXf &crossA) {
    assert(crossR.size() == crossA.cols());

    spectrumEucPolar.resize(nDir, 0.0);
    for(int i = 0; i < nDir; ++i) {
        for(int j = 0; j < crossA.cols(); ++j) {
            spectrumEucPolar[i] += (crossR[j] - crossA(i, j)) * (crossR[j] - crossA(i, j));
        }
    }
}
