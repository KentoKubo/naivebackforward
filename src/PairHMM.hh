#ifndef PAIRHMM_HH
#define PAIRHMM_HH

#include "Sequence.hh"
#include "PairHMMPar.hh"

namespace naive{
    class PairHMM{
        private:
        PairHMMPar par;
        void dumpMatrix(const std::string fileName, const std::vector<std::vector<float>>& matrix);

        public:
        float forward(const std::vector<uchar>& seq1, const std::vector<uchar>& seq2);
        float backward(const std::vector<uchar>& seq1, const std::vector<uchar>& seq2);
    };
}

#endif