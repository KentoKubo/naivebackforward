#include <iostream>
#include <iomanip>

#include "Sequence.hh"
#include "PairHMM.hh"

int main(int argc, char** argv){
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta file including 2 DNA sequences>" << std::endl;
        return 1;
    }
    
    std::cout << "INFO: reading fasta file: " << argv[1] << std::endl;
    
    naive::MultiSequences multiSequences;
    multiSequences.readFasta(argv[1]);
    
    naive::PairHMM pairHMM;
    float forwardScore = pairHMM.forward(multiSequences.getSequence(0).getEncodedSeq(), multiSequences.getSequence(1).getEncodedSeq());
    
    std::cout << std::fixed << std::setprecision(7) << "INFO: forward score: " << forwardScore << std::endl;
}