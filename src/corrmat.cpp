#include "corrmat.hpp"

#include "winfunc.hpp"

using namespace Eigen;

CorrMat::CorrMat(const int nft, fftwf_plan p) :
    nFFT(nft)
{
    int nFtBin(nFFT/2+1);

    in = (float*)fftwf_malloc(sizeof(float) * nFFT);
    out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nFtBin));

    if(p==nullptr) {
        ftPlanR2C = fftwf_plan_dft_r2c_1d(nFFT, in, out, FFTW_MEASURE);
        planCreated = true;
    } else {
        ftPlanR2C = p;
    }

    WinFunc w(nFFT); // hamming window
    win = std::move(w.win);
}

CorrMat::~CorrMat() {
    fftwf_free(in);
    fftwf_free(out);
    if(planCreated) {
        fftwf_cleanup();
        fftwf_destroy_plan(ftPlanR2C);
    }
}

void CorrMat::calc_corr(const Eigen::MatrixXf x) {
    // x: input sound; ch x samples
    nCh = x.rows();
    sndLen = x.cols();
    int nChunk(floor(sndLen / nFFT)), nFtBin(nFFT / 2 + 1),
            nPair(nCh * (nCh - 1) / 2);
    autocorrX.resize(nFtBin);
    for(Eigen::MatrixXcf& autocorrPerFreq : autocorrX){
        autocorrPerFreq.setZero(nCh, nCh);
    }

    Matrix<std::complex<float>, Dynamic, Dynamic, RowMajor> XX;
    XX.resize(nCh, nFtBin);

    for(int iChunk = 0; iChunk < nChunk; ++iChunk){
        for(int iCh = 0; iCh < nCh; ++iCh){
            for(int i = 0; i < nFFT; ++i){
                in[i] = x(iCh, iChunk*nFFT+i) * win[i];
            }
            fftwf_execute_dft_r2c(ftPlanR2C, in, out);
//            copy(out, out+nFtBin, XX.row(iCh).data());
            memcpy(XX.row(iCh).data(), out, nFtBin*sizeof(fftwf_complex));
        }
        // Correlation Matrix
        for(int iBin=0; iBin<nFtBin; ++iBin){
            autocorrX[iBin] += XX.col(iBin) * XX.col(iBin).adjoint();
        }
    }

    if(upperTriangular2Vec) {
        crossCorrX.resize(nFtBin);
        for(Eigen::VectorXcf& crossCorrPerFreq : crossCorrX) {
            crossCorrPerFreq.setZero(nPair);
        }

        crossCorrXAng.resize(nFtBin);
        for(Eigen::VectorXf& crossCorrPerFreqAng : crossCorrXAng) {
            crossCorrPerFreqAng.setZero(nPair);
        }

        for(int iBin = 0; iBin < nFtBin; ++iBin) {
            int iPair(0);
            for(int iCh1 = 0; iCh1 < nCh - 1; ++iCh1) {
                for(int iCh2 = iCh1 + 1; iCh2 < nCh; ++iCh2) {
                    assert(iPair < nPair);
                    crossCorrX[iBin][iPair] = autocorrX[iBin](iCh1, iCh2)
                            / abs(autocorrX[iBin](iCh1, iCh2));
                    crossCorrXAng[iBin][iPair] =
                            atan2(- crossCorrX[iBin][iPair].imag(),
                                  crossCorrX[iBin][iPair].real());
                    iPair++;
                }
            }
        }
    }
}
