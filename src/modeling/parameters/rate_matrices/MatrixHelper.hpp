#ifndef MATRIX_HELPER_HPP
#define MATRIX_HELPER_HPP
#include <array>

/**
 * @brief A namespace for miscellaneous functions/values that all codon rate matrices need. We have
 * written these so that they can either be evaluated at compile time or at startup.
 */
namespace MatrixHelper {
    static constexpr std::array<int, 61> aaMap = {8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 19, 19, 15, 15, 15, 15, 1, 18, 1, 9, 4, 9, 4};
    static constexpr std::array<const char*, 61> codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

    /**
     * @brief The pairs of coordinates in a codon matrix that are exactly 1 nucleotide away from each other
     */
    static const std::array<std::pair<int,int>, 263> validPairs = [] {
        std::array<std::pair<int,int>, 263> s;
        int currentIndex = 0;
        for(int i = 0; i < 61; i++){
            for(int j = i + 1; j < 61; j++){
                int mismatch = 0;
                bool isTransition = false;
                for(int k = 0; k < 3; k++){
                    if(codons[i][k] != codons[j][k]){
                        mismatch++;
                        if(mismatch > 1){
                            break;
                        }
                    }
                }
                if(mismatch == 1){
                    s[currentIndex] = {i, j};
                    currentIndex++;
                }
            }
        }

        return s;
    }();

    /**
     * @brief The pairs of coordinates in a codon matrix that are a transition mutation (as opposed to a transversion)
     */
    static const std::array<std::pair<int,int>, 89> transitionPairs = [] {
        std::array<std::pair<int,int>, 89> s;
        int currentIndex = 0;
        for(int i = 0; i < 61; i++){
            for(int j = i + 1; j < 61; j++){
                int mismatch = 0;
                bool isTransition = false;
                for(int k = 0; k < 3; k++){
                    if(codons[i][k] != codons[j][k]){
                        mismatch++;
                        if(mismatch > 1){
                            break;
                        }
                        if((codons[i][k] == 'A' && codons[j][k] == 'G') || (codons[i][k] == 'G' && codons[j][k] == 'A') || 
                        (codons[i][k] == 'T' && codons[j][k] == 'C') || (codons[i][k] == 'C' && codons[j][k] == 'T'))
                            isTransition = true;
                    }
                }
                if(mismatch == 1){
                    if(isTransition){
                        s[currentIndex] = {i, j};
                        currentIndex++;
                    }
                }
            }
        }
        return s;
    }();

    /**
     * @brief The pairs of coordinates in a codon matrix that are nonsynonymous
     */
    static const std::array<std::pair<int,int>, 67> synonymousPairs = [] {
        std::array<std::pair<int,int>, 67> s;
        int currentIndex = 0;
        for(int i = 0; i < 61; i++){
            for(int j = i + 1; j < 61; j++){
                int mismatch = 0;
                bool isTransition = false;
                for(int k = 0; k < 3; k++){
                    if(codons[i][k] != codons[j][k]){
                        mismatch++;
                        if(mismatch > 1){
                            break;
                        }
                    }
                }
                if(mismatch == 1){
                    if(aaMap[i] == aaMap[j]){
                        s[currentIndex] = {i, j};
                        currentIndex++;
                    }
                }
            }
        }
        return s;
    }();

    /**
     * @brief The pairs of coordinates in a codon matrix that are synonymous
     */
    static const std::array<std::pair<int,int>, 197> nonsynonymousPairs = [] {
        std::array<std::pair<int,int>, 197> s;
        int currentIndex = 0;
        for(int i = 0; i < 61; i++){
            for(int j = i + 1; j < 61; j++){
                int mismatch = 0;
                bool isTransition = false;
                for(int k = 0; k < 3; k++){
                    if(codons[i][k] != codons[j][k]){
                        mismatch++;
                        if(mismatch > 1){
                            break;
                        }
                    }
                }
                if(mismatch == 1){
                    if(aaMap[i] != aaMap[j]){
                        s[currentIndex] = {i, j};
                        currentIndex++;
                    }
                }
            }
        }
        return s;
    }();
}

#endif