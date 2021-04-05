#include <iostream>
#include <fstream>

#include "wavfile.hpp"
#include "corrmat.hpp"
#include "loclib/arraymanifold.hpp"
#include "loclib/dsloc.hpp"
#include "loclib/musicloc.hpp"

using namespace std;

int main()
{
    int nFFT = 2048; // 64

    // calcultate manifold vectors
//    ArrayManifold arrManif("iphonexsmax.txt", nFFT);
//    arrManif.generate_spherical_mesh();
//    arrManif.generate_a_mat();

    // measured manifold vectors
    ArrayManifold arrManif(2048);
    arrManif.read_ir_files("../../data/ir/iphonexsmax/ir-xsmax-16k/ir000.wav",
                           360, 1);

    arrManif.calc_cross_a();

//    std::vector<Eigen::MatrixXcf> aMat = arrManif.downsample_a_mat(32);
//    std::vector<Eigen::MatrixXcf> aVecCross = arrManif.downsample_a_vec_cross(32);
//    std::vector<Eigen::MatrixXf> aVecCrossAng = arrManif.downsample_a_vec_cross_ang(32);

    std::string filename("hello10degxsmax16k.wav");
    WavFile f(filename);

    CorrMat corr(nFFT);
    corr.upperTriangular2Vec = true;

    std::fstream resFile;
    resFile.open("localization-results.csv", std::fstream::out);

    int nSec(0);
    while(!f.eof) {
        std::vector<float> spatialSpec(arrManif.nDir, 0.0);
        Eigen::MatrixXf sndBlock = f.next_n_samples(16000);
        corr.calc_corr(sndBlock);

        std::vector<float> spatialSpecCosCart(arrManif.nDir, 0.0);
        std::vector<float> spatialSpecEucCart(arrManif.nDir, 0.0);
        std::vector<float> spatialSpecCosPolar(arrManif.nDir, 0.0);
        std::vector<float> spatialSpecEucPolar(arrManif.nDir, 0.0);

        for(int i = 1; i < nFFT / 2 + 1; ++i) {
            DSLoc locDs(corr.autocorrX[i], arrManif.aMat[i]);
            std::transform(spatialSpec.begin(), spatialSpec.end(),
                           locDs.spectrum.begin(), spatialSpec.begin(),
                           std::plus<float>());

//            locDs.calc_cos_sim_cartesian(corr.crossCorrX[i], arrManif.aVecCross[i]);
//            std::transform(spatialSpecCosCart.begin(), spatialSpecCosCart.end(),
//                           locDs.spectrumCosCart.begin(), spatialSpecCosCart.begin(),
//                           std::plus<float>());

            MusicLoc locMusic(corr.autocorrX[i], arrManif.aMat[i], 1);
            std::transform(spatialSpecCosCart.begin(), spatialSpecCosCart.end(),
                           locMusic.spectrum.begin(), spatialSpecCosCart.begin(),
                           std::plus<float>());

            locDs.calc_euc_dist_cartesian(corr.crossCorrX[i], arrManif.aVecCross[i]);
            std::transform(spatialSpecEucCart.begin(), spatialSpecEucCart.end(),
                           locDs.spectrumEucCart.begin(), spatialSpecEucCart.begin(),
                           std::plus<float>());

            locDs.calc_cos_sim_polar(corr.crossCorrXAng[i], arrManif.aVecCrossAng[i]);
            std::transform(spatialSpecCosPolar.begin(), spatialSpecCosPolar.end(),
                           locDs.spectrumCosPolar.begin(), spatialSpecCosPolar.begin(),
                           std::plus<float>());

            locDs.calc_euc_dist_polar(corr.crossCorrXAng[i], arrManif.aVecCrossAng[i]);
            std::transform(spatialSpecEucPolar.begin(), spatialSpecEucPolar.end(),
                           locDs.spectrumEucPolar.begin(), spatialSpecEucPolar.begin(),
                           std::plus<float>());
        }

        for(float v : spatialSpec) {
            resFile << v / 32 << ", ";
        }
        resFile << std::endl;

        for(float v : spatialSpecCosCart) {
            resFile << v / 32 << ", ";
        }
        resFile << std::endl;

        for(float v : spatialSpecEucCart) {
            resFile << v / 32 << ", ";
        }
        resFile << std::endl;

        for(float v : spatialSpecCosPolar) {
            resFile << v / 32 << ", ";
        }
        resFile << std::endl;

        for(float v : spatialSpecEucPolar) {
            resFile << v / 32 << ", ";
        }
        resFile << std::endl;

        resFile << std::endl;

        std::cout << nSec << " seconds finished." << std::endl;
        nSec++;
    }

    resFile.close();

    return 0;
}
