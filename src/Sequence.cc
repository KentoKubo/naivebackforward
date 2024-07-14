#include <fstream>
#include <algorithm>
#include <stdexcept>
#include "Sequence.hh"

naive::Sequence::Sequence(const uint& seqId, const std::string& seqName, const std::string& sequence){
    this->seqId = seqId;
    this->seqName = seqName;
    this->encodedSeq.reserve(sequence.size());

    for(const char &c : sequence){
        this->encodedSeq.emplace_back(naive::CharCoder::encode(c));
    }
}

void naive::MultiSequences::readFasta(const std::string& filePath){
    int seqCtr = 0;
    std::string seqName = "";
    std::string seqStr = "";

    std::string line = "";
    std::ifstream ifs(filePath);

    if(ifs.fail()){
        throw std::runtime_error("File reading error: cannot open file -> " + filePath);
    }

    while (getline(ifs, line)){
        if(line.length() == 0) continue;
        if(line[0] == '>'){
            if(seqStr.empty()) {
                seqName = line.substr(1);
                continue;
            }
            std::transform(seqStr.cbegin(), seqStr.cend(), seqStr.begin(), toupper);
            multiSequences.emplace_back(Sequence(seqCtr, seqName, seqStr));
            ++seqCtr;
            seqName = line.substr(1);
            seqStr = "";
        }
        else {
            seqStr += line;
        }
    }
    if(!seqStr.empty()){
        std::transform(seqStr.cbegin(), seqStr.cend(), seqStr.begin(), toupper);
        multiSequences.emplace_back(Sequence(seqCtr, seqName, seqStr));
    }

    seqNum = seqCtr + 1;
}
