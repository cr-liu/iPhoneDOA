#include "arraymanifold.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <sndfile.hh>

#include "../winfunc.hpp"

ArrayParam::ArrayParam(std::string filename)
{
    std::fstream geoFile;
    geoFile.open(filename, std::fstream::in);
    if(!geoFile.is_open()) {
        std::cout<< "failed open file: " << filename << std::endl;
        return;
    }

    std::string lineStr;
    int iRow(0);
    while(std::getline(geoFile, lineStr)) {
        if(lineStr.length() < 2) {
            continue;
        }
        std::stringstream ss(lineStr);
        if(arrayName == "") {
            ss >> arrayName;
            continue;
        } else if(nSensor == 0 ) {
            ss >> nSensor;
            arrayGeo.resize(nSensor, 3);
            iRow = 0;
            continue;
        }
        for(int j=0; j<3; ++j) {
            float tmpF;
            ss >> tmpF;
            arrayGeo(iRow, j) = tmpF;
        }
        if(++iRow == arrayGeo.rows()) {
            break;
        }
    }

    // std::cout << arrayGeo << std::endl;
    geoFile.close();
}

ArrayParam::ArrayParam()
{
}

ArrayManifold::ArrayManifold(std::string geoFilename,int nft,
                             int sr, float distance) :
    distanceNear(distance),
    sampRate(sr),
    nfft(nft)
{
    param = new ArrayParam(geoFilename);
}

ArrayManifold::ArrayManifold(int nft) :
    nfft(nft)
{
    param = new ArrayParam();
    aMat.resize(nfft / 2 + 1);
}


ArrayManifold::~ArrayManifold() {
    delete param;
    if(ftPlanCreated) {
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
        fftwf_destroy_plan(ftPlanR2C);
    }
}

void ArrayManifold::generate_spherical_mesh(float aziL, float aziH, float aziStep,
                                            float eleL, float eleH, float eleStep) {
    dirMesh.reserve(std::max(1, (int)ceil((aziH - aziL) / aziStep)) *
                    std::max(1, (int)ceil((eleH - eleL) / eleStep)));
    int dirCnt(0), nPoint, i;
    float aziStepHat;
    float aziLRad, aziHRad, aziStepRad, aux;

    aziLRad = aziL * M_PI / 180;
    aziHRad = aziH * M_PI / 180;
    aziStepRad = aziStep * M_PI / 180;

    dirCnt = 0;
    for(int elevation=(int)eleL; elevation<=(int)eleH; elevation+=(int)eleStep) {
        aux = (aziHRad - aziLRad) * cos(elevation) / aziStepRad;
        nPoint = ceil(aux - 0.0001);
        aziStepHat = (aziHRad - aziLRad) / nPoint;

        for (i=0; i<nPoint; i++) {
            std::array<float, 2> coord;
            coord[0] = aziLRad + aziStepHat*i;
            coord[1] = (float)elevation * M_PI / 180;
            dirMesh.push_back(std::move(coord));
            dirCnt++;
        }
    }
    nDir = dirCnt;
    std::cout << "num of dirs: " << nDir << std::endl;
}

void ArrayManifold::generate_a_mat() {
    aMat.resize(nfft / 2 + 1);
    for(Eigen::MatrixXcf& aVecs : aMat) {
        aVecs.resize(nDir, param->nSensor);
    }

    for(int i=0; i<nDir; ++i) {
        Eigen::MatrixXcf aVecs =
                    compute_target_a(dirMesh[i][0], dirMesh[i][1]);
        for(int j=0; j<nfft/2+1; ++j) {
            aMat[j].row(i) = aVecs.row(j);
        }
    }
}

Eigen::MatrixXcf ArrayManifold::compute_target_a(const float azi, const float ele) {
    Eigen::MatrixXcf targetA;
    targetA.resize(nfft/2+1, param->nSensor);
    float omega, tau, dist;
    int i;
    Eigen::Array3f source;
    Eigen::Matrix3f pm; // projection matrix; onto line between origin to source
    Eigen::Array3f pv; // projection of sensor vector onto line from origin to source

    float ce(cos(ele));
    if(distanceNear>0) { // near mode
        source[0] = distanceNear * cos(azi) * ce;
        source[1] = distanceNear * sin(azi) * ce;
        source[2] = distanceNear * sin(ele);
    }else{ // far mode
        source[0] = cos(azi) * ce;
        source[1] = sin(azi) * ce;
        source[2] = sin(ele);
    }

    float normfactor;
    normfactor = source.square().sum();
    pm = (source.matrix() * source.matrix().transpose()) / normfactor;

    for(int iFreq=0; iFreq<nfft/2+1; ++iFreq) {
        omega = 2 * M_PI * sampRate * iFreq / nfft;
        if (distanceNear==0){ // far mode
            for (i=0; i<param->nSensor; ++i) {
                pv.matrix() = pm * param->arrayGeo.row(i).transpose();
                normfactor = sqrt(pv.square().sum());
                // (signed) length of the projection
                if((pv * source).sum() > 0) {
                    tau = -normfactor / sndVelocity;
                } else {
                    tau = normfactor / sndVelocity;
                }
                targetA(iFreq, i) = std::complex<float>(
                            cos(omega*tau), -sin(omega*tau));
            }
        } else { // near mode
            for (i=0; i<param->nSensor; i++) {
                dist = sqrt((param->arrayGeo.row(i).transpose().array()
                             - source).square().sum());
                tau = (dist - distanceNear) / sndVelocity;
                targetA(iFreq, i) = std::complex<float>(
                            1/dist*cos(omega*tau), -1/dist*sin(omega*tau));
            }
        }
    }
    return targetA;
}

void ArrayManifold::calc_cross_a() {
    int nPair(param->nSensor * (param->nSensor - 1) / 2), nFreq(nfft / 2 + 1);

    aVecCross.resize(nFreq);
    for(Eigen::MatrixXcf& aVecCorr : aVecCross) {
        aVecCorr.resize(nDir, nPair);
    }

    aVecCrossAng.resize(nFreq);
    for(Eigen::MatrixXf& aVecCorrAng : aVecCrossAng) {
        aVecCorrAng.resize(nDir, nPair);
    }

    for(int iDir = 0; iDir < nDir; ++iDir) {
        for(int iFreq = 0; iFreq < nFreq; ++iFreq) {
            int iPair(0);
            for(int iCh1 = 0; iCh1 < param->nSensor - 1; ++iCh1) {
                for(int iCh2 = iCh1 + 1; iCh2 < param->nSensor; ++iCh2) {
                    assert(iPair < nPair);
                    aVecCross[iFreq](iDir, iPair) = aMat[iFreq](iDir, iCh1) *
                            std::conj(aMat[iFreq](iDir, iCh2));
                    aVecCrossAng[iFreq](iDir, iPair) =
                            atan2(aVecCross[iFreq](iDir, iPair).imag(),
                            aVecCross[iFreq](iDir, iPair).real());
                    iPair++;
                }
            }
        }
    }
}

void ArrayManifold::read_ir_files(const std::string filename, int angleEnd,
                                  int angleStep) {
    size_t dotPos = filename.find_last_of(".");
    std::string preStr, numStr, postStr;
    preStr = filename.substr(0, dotPos - 3);
    numStr = filename.substr(dotPos - 3, 3);
    postStr = filename.substr(dotPos);

    int degree;
    try {
        degree = stoi(numStr);
    }  catch (...) {
        std::cout << "Filename is not end with 3 digit number!" << std::endl;
    }

    in = (float*)fftwf_malloc(sizeof(float) * nfft);
    out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));
    ftPlanR2C = fftwf_plan_dft_r2c_1d(nfft, in, out, FFTW_MEASURE);
    ftPlanCreated = true;

    nDir = (angleEnd - degree) / angleStep;

    int iDir(0);
    for(int deg = degree; deg < angleEnd; deg += angleStep) {
        numStr = std::to_string(deg);
        while(numStr.length() < 3) {
            numStr = "0" + numStr;
        }
        read_ir(preStr + numStr + postStr, iDir);
        iDir++;
    }
}

void ArrayManifold::read_ir(const std::string filename, int iDir) {
    SndfileHandle file(filename);
    if(file.channels() == 0) {
        std::cout << filename + " not found!" << std::endl;
        return;
    }
    int fileCh = file.channels(), fileSr = file.samplerate();
    if(param->nSensor < 2) {
        for(Eigen::MatrixXcf& aVecs : aMat) {
            param->nSensor = fileCh;
            sampRate = fileSr;
            aVecs.resize(nDir, fileCh);
        }
    }

    assert(fileCh == param->nSensor);

    Eigen::MatrixXf snd;
    snd.resize(fileCh, file.frames());
    file.read(snd.data(), file.channels() * file.frames());

    exe_fft(snd, iDir);
}

void ArrayManifold::exe_fft(const Eigen::MatrixXf& data, int iDir) {
    WinFunc w(nfft);
    win = std::move(w.win);

    for(int iCh = 0; iCh < data.rows(); ++iCh) {
        for(int i = 0; i < nfft; ++i) {
            in[i] = data(iCh, i) * win[i];
        }
        fftwf_execute_dft_r2c(ftPlanR2C, in, out);
        for(int j = 0; j < (nfft / 2 + 1); ++j) {
            float norm = sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]);
            aMat[j](iDir, iCh) = std::complex<float>(out[j][0] / norm,
                    -out[j][1] / norm);
        }
    }
}

std::vector<Eigen::MatrixXcf> ArrayManifold::downsample_a_mat(int factor) {
    assert(factor & (factor - 1) == 0);

    std::vector<Eigen::MatrixXcf> res;
    for(uint i = 0; i < aMat.size(); i += factor) {
        res.push_back(aMat[i]);
    }

    return res;
}

std::vector<Eigen::MatrixXcf> ArrayManifold::downsample_a_vec_cross(int factor) {
    assert(factor & (factor - 1) == 0);

    std::vector<Eigen::MatrixXcf> res;
    for(uint i = 0; i < aMat.size(); i += factor) {
        res.push_back(aVecCross[i]);
    }

    return res;
}

std::vector<Eigen::MatrixXf> ArrayManifold::downsample_a_vec_cross_ang(int factor) {
    assert(factor & (factor - 1) == 0);

    std::vector<Eigen::MatrixXf> res;
    for(uint i = 0; i < aMat.size(); i += factor) {
        res.push_back(aVecCrossAng[i]);
    }

    return res;
}
