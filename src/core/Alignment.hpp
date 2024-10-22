#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP
#include <string>
#include "ncl/nxsmultiformat.h"

class Alignment{
    public:
                                    Alignment(void) = delete;
                                    Alignment(std::string fn);
                                    ~Alignment();
        unsigned long long int      getCharCode(int i, int j) {return matrix[i][j];}
        unsigned long long int**    getMatrix() {return matrix;}
        int                         getNumChar() {return numChar;}
        int                         getNumTaxa() {return numTaxa;}
        int                         getStateSpace() {return stateSpace;}
        std::vector<std::string>    getTaxaNames() {return taxaNames;}
        std::vector<double>         getStateFrequencies() {return frequencies;}
    private:
        unsigned long long int**    matrix;
        int                         numTaxa;
        int                         numChar;
        int                         stateSpace;
        void                        readCodonData(NxsCharactersBlock* charBlock);
        std::vector<std::string>    taxaNames;
        std::vector<double>         frequencies;
};

#endif