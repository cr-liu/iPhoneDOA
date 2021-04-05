#ifndef ARRAYMANIFOLD_HPP
#define ARRAYMANIFOLD_HPP

#include <vector>
#include <array>
#include <string>

#include <eigen3/Eigen/Dense>
#include <fftw3.h>

class ArrayParam
{
public:
    ArrayParam(std::string filename);
    ArrayParam(); // for recorded impulse response; arrayGeo will not init
    std::string arrayName = "";
    Eigen::MatrixXf arrayGeo; // nSensor x 3; x, y, z in meter
    int nSensor = 0;
};

class ArrayManifold
{
public:
    ArrayManifold(std::string geoFilename, int nft,
                  int sr=16000, float distance=0.0);
    ArrayManifold(int nft); // for recorded impulse response
    ~ArrayManifold();
    ArrayParam* param=nullptr;
    float distanceNear = 0;
    int sampRate;
    int nfft;
    int nDir;
    float sndVelocity = 343.0;
    std::vector<std::array<float, 2>> dirMesh; // azi & ele for every dir
    std::vector<Eigen::MatrixXcf> aMat; // freq x nDir x nSensor
    std::vector<Eigen::MatrixXcf> aVecCross; // freq x nDir x nPair
    std::vector<Eigen::MatrixXf> aVecCrossAng; // freq x nDir x nPair
    bool ftPlanCreated = false;
    float* in = nullptr;
    fftwf_complex* out = nullptr;
    fftwf_plan ftPlanR2C = nullptr;
    Eigen::ArrayXf win;

    void generate_spherical_mesh(float aziL=0.0, float aziH=360.0,
                                 float aziStep=1.0, float eleL=0.0,
                                 float eleH=0.0,  float eleStep=1.0);
    void generate_a_mat();
    void calc_cross_a();
    Eigen::MatrixXcf compute_target_a(const float azi, const float ele);
    void read_ir_files(const std::string filename, int angleEnd, int angleStep);
    void read_ir(const std::string filename, int iDir);
    void exe_fft(const Eigen::MatrixXf& data, int iDir);

    std::vector<Eigen::MatrixXcf> downsample_a_mat(int factor);
    std::vector<Eigen::MatrixXcf> downsample_a_vec_cross(int factor);
    std::vector<Eigen::MatrixXf> downsample_a_vec_cross_ang(int factor);
};

#endif // ARRAYMANIFOLD_HPP
