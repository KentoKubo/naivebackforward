#include <iostream>

#include "Sequence.hh"

int main(int argc, char** argv){
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta file including 2 DNA sequences>" << std::endl;
        return 1;
    }
    
    std::cout << "INFO: reading fasta file: " << argv[1] << std::endl;
    
    naive::MultiSequences multiSequences;
    multiSequences.readFasta(argv[1]);
    
    naive::Sequence seq1 = multiSequences.getSequence(0);
    std::cout << "INFO: seq1 len" << seq1.getEncodedSeq().size() << std::endl;
    
}