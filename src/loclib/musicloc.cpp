#include "musicloc.hpp"

#include <iostream>

using namespace Eigen;

MusicLoc::MusicLoc(const Eigen::MatrixXcf &R,
                   const Eigen::MatrixXcf &A, int nSignal) :
    nSensor(A.cols()),
    nDir(A.rows())
{
    // MUSIC Spatial Spectrum
    spectrum.resize(nDir, 0.0);

    SelfAdjointEigenSolver<MatrixXcf> eigenSolver(R);
    if(eigenSolver.info() != Success) {
        std::cout << "Eigen solver failed!" << std::endl;
        return;
    }

    eigenValues = eigenSolver.eigenvalues();
    eigenVectors = eigenSolver.eigenvectors();

    for(int i=0; i<nDir; ++i){
        for(int n=0; n<nSensor-nSignal; ++n){
            std::complex<double> d;
            d = A.row(i).conjugate().dot(eigenVectors.col(n));
            /*
            for(m=0; m<M; ++m){
                d += conj(A[i][m])*EigenVector[m][n];
            }
            */
            spectrum[i] += (conj(d) * d).real();
        }
        spectrum[i] = 1.0 / spectrum[i];
    }
}
