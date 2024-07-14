#ifndef SEQUENCE_HH
#define SEQUENCE_HH

#include <string>
#include <vector>
#include <cassert>

namespace naive{
    using uchar = unsigned char;
    using uint = unsigned int;
    using ulong = unsigned long;
    using llong = long long int;

    class CharCoder{
        public:
        static inline uchar encode(const char& c){
            switch(c){
            case 'N':
                return 0;
            case 'A':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 3;
            case 'T':
                return 4;
            case 'U':
                return 4;
            default:
                return 5;
            }
        }

        static inline char decode(const uchar& u){
            switch(u){
            case 0:
                return 'N';
            case 1:
                return 'A';
            case 2:
                return 'C';
            case 3:
                return 'G';
            case 4:
                return 'U';
            default:
                return '*';
            }
        }
    };

    class Sequence{
        private:
        uint seqId;
        std::string seqName;
        std::vector<uchar> encodedSeq;

        public:
        Sequence(const uint& seqId, const std::string& seqName, const std::string& sequence);

        const uint getSeqId() const {
            return seqId;
        }

        const std::string& getSeqName() const {
            return seqName;
        }
        
        const std::vector<uchar>& getEncodedSeq() const {
            return encodedSeq;
        }
    };

    class MultiSequences{
        private:
        uint seqNum;
        std::vector<Sequence> multiSequences;

        public:
        const uint getSeqNum() const {
            return seqNum;
        }

        const Sequence& getSequence(const uint& seqId) const {
            return multiSequences[seqId];
        }

        void readFasta(const std::string& filePath);
    };
}

#endif
