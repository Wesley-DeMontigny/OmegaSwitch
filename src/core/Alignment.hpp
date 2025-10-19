#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP
#include <string>
#include "ncl/nxsmultiformat.h"

/**
 * @brief Uses the nexus class library to load a NEXUS file into a format that respects ambiguous
 * characters. To be used to initialize a ConditionalLikelihood object. Theoretically, this could just
 * take in an alignment and convert codons to integers representing an ID, but we want to account for 
 * all possible ambiguity codes. In a perfect world this would be able to identify potential codons with
 * the middle nucleotide as a gap (due to sequencing error). Right now we go a little bit overkill by setting
 * bit flags in an unsigned long long int matrix.
 */
class Alignment{
    public:
                                    Alignment(void) = delete;
                                    Alignment(std::string fn);                                  // Costructor from a NEXUS file
                                    Alignment(int* siteMatrix, int numChar, int numTaxa);       // Constructor from a simulated site matrix
                                    ~Alignment();                                               // Destructor

        int                         getNumChar() {return numChar;}                              // Returns the number of characters (sites) in the alignment
        int                         getNumTaxa() {return numTaxa;}                              // Returns the number of ta
        int                         getStateSpace() {return stateSpace;}                        // Returns the size of the state space (61)
        std::vector<std::string>    getTaxaNames() {return taxaNames;}                          // Returns the vector of taxa names from the alignment
        unsigned long long int      getCharCode(int i, int j) {return matrix[i][j];}            // Returns entry ij
        unsigned long long int**    getMatrix() {return matrix;}                                // Returns the whole character matrix
    private:
        void                        readCodonData(NxsCharactersBlock* charBlock);               // Read the codon data from the NEXUS character block into the matrix
        
        int                         numChar;                                                    // The number of characters (sites) in the alignment
        int                         numTaxa;                                                    // The number of taxa in the alignment
        int                         stateSpace;                                                 // The size of the state space (61)
        std::vector<std::string>    taxaNames;                                                  // The taxa names in the alignment
        unsigned long long int**    matrix;                                                     // The matrix containing the bit flags for each entry in the alignment
};

#endif