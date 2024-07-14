#include "PairHMM.hh"
#include "PairHMMPar.hh"

float naive::PairHMM::forward(const std::vector<uchar>& seq1, const std::vector<uchar>& seq2){
    // initialize
    llong seq1Len = seq1.size();
    llong seq2Len = seq2.size();

    // all f*(i,-1), f*(-1,j) = 0
    std::vector<std::vector<float>> match(seq1Len+2, std::vector<float>(seq2Len+2, xlog(0.0)));
    std::vector<std::vector<float>> insert1(seq1Len+2, std::vector<float>(seq2Len+2, xlog(0.0)));
    std::vector<std::vector<float>> insert2(seq1Len+2, std::vector<float>(seq2Len+2, xlog(0.0)));
    
    // fM(0,0) = 1, fI1(0,0) = fI2(0,0) = 0
    match[1][1] = xlog(1.0);
    
    // recursion of i = 1 or j = 1
    for(llong i = 2; i < seq1Len+2; i++){
        match[i][1] = par.getLogEmitProb(seq1[i-2], 0, State::MATCH) 
            + xlogsumexp(
                par.getLogTransProb(State::MATCH, State::MATCH) + match[i-1][0], 
                par.getLogTransProb(State::INS1, State::MATCH) + insert1[i-1][0], 
                par.getLogTransProb(State::INS2, State::MATCH) + insert2[i-1][0]
                );
        insert1[i][1] = par.getLogEmitProb(seq1[i-2], 0, State::INS1) // emit prob is 0 because seq1 is a gap in this state
            + xlogsumexp(
                par.getLogTransProb(State::MATCH, State::INS1) + match[i][0],
                par.getLogTransProb(State::INS1, State::INS1) + insert1[i][0]
                );
        insert2[i][1] = par.getLogEmitProb(seq1[i-2], 0, State::INS2) 
            + xlogsumexp(
                par.getLogTransProb(State::MATCH, State::INS2) + match[i-1][1], 
                par.getLogTransProb(State::INS2, State::INS2) + insert2[i-1][1]
                );
    }
    
    for(llong j = 2; j < seq2Len+2; j++){
        match[1][j] = par.getLogEmitProb(0, seq2[j-2], State::MATCH) 
            + xlogsumexp(
                par.getLogTransProb(State::MATCH, State::MATCH) + match[0][j-1], 
                par.getLogTransProb(State::INS1, State::MATCH) + insert1[0][j-1], 
                par.getLogTransProb(State::INS2, State::MATCH) + insert2[0][j-1]
                );
        insert1[1][j] = par.getLogEmitProb(0, seq2[j-2], State::INS1) 
            + xlogsumexp(
                par.getLogTransProb(State::MATCH, State::INS1) + match[1][j-1],
                par.getLogTransProb(State::INS1, State::INS1) + insert1[1][j-1]
                );
        insert2[1][j] = par.getLogEmitProb(0, seq2[j-2], State::INS2) // emit prob is 0 because seq2 is a gap in this state
            + xlogsumexp(
                par.getLogTransProb(State::MATCH, State::INS2) + match[0][j], 
                par.getLogTransProb(State::INS2, State::INS2) + insert2[0][j]
                );
    }

    // recursion else
    for (llong i = 2; i < seq1Len+2; i++){
        for (llong j = 2; j < seq2Len+2; j++){
            match[i][j] = par.getLogEmitProb(seq1[i-2], seq2[j-2], State::MATCH) 
                + xlogsumexp(
                    par.getLogTransProb(State::MATCH, State::MATCH) + match[i-1][j-1], 
                    par.getLogTransProb(State::INS1, State::MATCH) + insert1[i-1][j-1], 
                    par.getLogTransProb(State::INS2, State::MATCH) + insert2[i-1][j-1]
                    );
            insert1[i][j] = par.getLogEmitProb(seq1[i-2], seq2[j-2], State::INS1) 
                + xlogsumexp(
                    par.getLogTransProb(State::MATCH, State::INS1) + match[i][j-1],
                    par.getLogTransProb(State::INS1, State::INS1) + insert1[i][j-1]
                    );
            insert2[i][j] = par.getLogEmitProb(seq1[i-2], seq2[j-2], State::INS2) 
                + xlogsumexp(
                    par.getLogTransProb(State::MATCH, State::INS2) + match[i-1][j], 
                    par.getLogTransProb(State::INS2, State::INS2) + insert2[i-1][j]
                    );
        }
    }

    // termination
    float score = xlogsumexp(match[seq1Len+1][seq2Len+1], insert1[seq1Len+1][seq2Len+1], insert2[seq1Len+1][seq2Len+1]);
}