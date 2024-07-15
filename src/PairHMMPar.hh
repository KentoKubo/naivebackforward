#ifndef _PAIRHMM_PAR_HH_
#define _PAIRHMM_PAR_HH_

#include "Utils.hh"
#include "Sequence.hh"

namespace naive{
	const int N_STATES = 3;
	const int N_OUTPUTS = 25;
	const int N_CHARS = 5;
	
	// INS1 is the insertion for seq1, so seq1 is a gap and seq2 is a char
	// INS2 is the insertion for seq2, so seq2 is a gap and seq1 is a char
    enum struct State{ INS2, INS1, MATCH };

	class PairHMMPar{
		private:
		double TransProbs[N_STATES][N_STATES] = {
			{0.666439, 0.041319, 0.292242}, // INS2
			{0.041319, 0.666439, 0.292242}, // INS1
			{0.022666, 0.022666, 0.954668}  // MATCH
		};
		double EmitProbs[N_OUTPUTS][N_STATES] = {
			{0.000000, 0.000000, 1.000000}, // ..
			{0.000000, 0.211509, 0.000000}, // .A
			{0.000000, 0.257349, 0.000000}, // .C
			{0.000000, 0.271398, 0.000000}, // .G
			{0.000000, 0.259744, 0.000000}, // .U
			{0.211509, 0.000000, 0.000000}, // A.
			{0.000000, 0.000000, 0.134009}, // AA
			{0.000000, 0.000000, 0.027164}, // AC
			{0.000000, 0.000000, 0.049659}, // AG
			{0.000000, 0.000000, 0.028825}, // AU
			{0.257349, 0.000000, 0.000000}, // C.
			{0.000000, 0.000000, 0.027164}, // CA
			{0.000000, 0.000000, 0.140242}, // CC
			{0.000000, 0.000000, 0.037862}, // CG
			{0.000000, 0.000000, 0.047735}, // CU
			{0.271398, 0.000000, 0.000000}, // G.
			{0.000000, 0.000000, 0.049659}, // GA
			{0.000000, 0.000000, 0.037862}, // GC
			{0.000000, 0.000000, 0.178863}, // GG
			{0.000000, 0.000000, 0.032351}, // GU
			{0.259744, 0.000000, 0.000000}, // U.
			{0.000000, 0.000000, 0.028825}, // UA
			{0.000000, 0.000000, 0.047735}, // UC
			{0.000000, 0.000000, 0.032351}, // UG
			{0.000000, 0.000000, 0.099694}, // UU
		};
		
		public:
		double getLogTransProb(State from, State to){
			return xlog(TransProbs[static_cast<int>(from)][static_cast<int>(to)]);
		}
		double getLogEmitProb(const uchar& c1, const uchar& c2, State state){
			return xlog(EmitProbs[c1 * N_CHARS + c2][static_cast<int>(state)]);
		}
	};
}

#endif