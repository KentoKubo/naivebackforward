#ifndef PAIRHMM_HH
#define PAIRHMM_HH

#include "Sequence.hh"
#include "PairHMMPar.hh"

namespace naive{
    class PairHMM{
        public:
        float forward(const Sequence& seq1, const Sequence& seq2);
        float backward(const Sequence& seq1, const Sequence& seq2);
    };
}

#endif