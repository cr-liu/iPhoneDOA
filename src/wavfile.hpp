#ifndef WAVFILE_HPP
#define WAVFILE_HPP

#include <string>

#include <eigen3/Eigen/Dense>

class WavFile
{
public:
    WavFile(const std::string filename);
    WavFile(const std::string filename, int n); // read multiple 1ch files
    Eigen::MatrixXf next_n_samples(int n);

    int nCh;
    int sr = 0;
    int format;
    Eigen::MatrixXf sndData; // ch x samples
    bool eof=false;

private:
    void read_wav_file(const std::string filename, int chIndex=-1);
    int currIndex = 0; // used as iterator
};

#endif // WAVFILE_HPP
