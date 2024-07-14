#ifndef PAIRHMM_HH
#define PAIRHMM_HH

#include "Sequence.hh"
#include "PairHMMPar.hh"

namespace naive{
    class PairHMM{
        private:
        PairHMMPar par;

        public:
        float forward(const std::vector<uchar>& seq1, const std::vector<uchar>& seq2);
        float backward(const std::vector<uchar>& seq1, const std::vector<uchar>& seq2);
    };
}

#endif