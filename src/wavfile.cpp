#include "wavfile.hpp"

#include <iostream>
#include <sndfile.hh>

WavFile::WavFile(const std::string filename)
{
    read_wav_file(filename);
}

WavFile::WavFile(const std::string filename, int n) :
    nCh(n)
{
    std::size_t dotPos = filename.find_last_of(".");
    std::string preStr, numStr, postStr, wavFilename;
    preStr = filename.substr(0, dotPos - 2);
    numStr = filename.substr(dotPos - 2, 2);
    postStr = filename.substr(dotPos);

    int fileNum;
    try {
        fileNum = std::stoi(numStr);
    }  catch (...) {
        std::cout << "The last 2 chars of filename are not numbers!" << std::endl;
        return;
    }

    for(int i=0; i<n; ++i) {
        numStr = std::to_string(fileNum + i);
        if(numStr.length() < 2) {
            numStr = "0" + numStr;
        }
        wavFilename = preStr + numStr + postStr;
        std::cout << wavFilename << std::endl;
        read_wav_file(wavFilename, i);
    }
}

void WavFile::read_wav_file(const std::string filename, int chIndex) {
    SndfileHandle file;
    file = SndfileHandle(filename);

    // check nullptr caused crash! check channels instead
    if(file.channels() == 0) {
        std::cout << filename + " not found!" << std::endl;
        return;
    }
    int fileCh = file.channels(), fileSr = file.samplerate(), fileFmt = file.format();
    int nSample = file.frames();
    if(chIndex < 0) {
        nCh = fileCh;
        sr = fileSr;
        format = fileFmt;
        sndData.resize(nCh, nSample);
        file.read(sndData.data(), file.frames() * nCh);
    } else if(chIndex == 0 && fileCh == 1) {
        sr = fileSr;
        format = fileFmt;
        sndData.resize(nCh, nSample);
        Eigen::RowVectorXf snd(sndData.cols());
        file.read(snd.data(), file.frames());
        sndData.row(chIndex) = snd;
    } else if(fileCh == 1 && fileSr == sr && fileFmt == format) {
        Eigen::RowVectorXf snd(sndData.cols());
        file.read(snd.data(), file.frames());
        sndData.row(chIndex) = snd;
    } else {
        std::cout << "Expecting 1ch files with same format!" << std::endl;
    }

    // test
    //    Eigen::VectorXd ch1 = sndData.row(0);
    //    file = SndfileHandle("test.wav", SFM_WRITE, format, 1, sr);
    //    file.write(ch1.data(), ch1.size());
}

Eigen::MatrixXf WavFile::next_n_samples(int n) {
    std::cout << currIndex << "th sample of " << sndData.cols() << std::endl;
    if(currIndex + n < sndData.cols()) {
        currIndex += n;
        return sndData.block(0, currIndex - n, nCh, n);
    } else {
        eof = true;
        return sndData.block(0, currIndex, nCh, sndData.cols()-currIndex);
    }
}
