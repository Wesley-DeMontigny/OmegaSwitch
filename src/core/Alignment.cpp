#include "Alignment.hpp"
#include "Msg.hpp"
#include "ncl/nxscharactersblock.cpp"

/**
 * @brief Construct a new alignment from the provided NEXUS file.
 * 
 */
Alignment::Alignment(std::string fn) {
    
    MultiFormatReader nexusReader;
    const char* fileName = fn.c_str();
    nexusReader.ReadFilepath(fileName, MultiFormatReader::NEXUS_FORMAT);

    size_t numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    if(numTaxaBlocks > 1)
        Msg::error("Too many taxa blocks (> 1)");

    for(size_t tBlock = 0; tBlock < numTaxaBlocks; tBlock++){
        NxsTaxaBlock* taxaBlock = nexusReader.GetTaxaBlock(tBlock);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        const unsigned numCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        const unsigned numUnalignedCharBlocks = nexusReader.GetNumUnalignedBlocks(taxaBlock);
        
        if(numUnalignedCharBlocks > 0)
            Msg::error("No unaligned data allowed!");
        if(numCharBlocks > 1)
            Msg::error("Too many char blocks (> 1)");

        for(size_t cBlock = 0; cBlock < numCharBlocks; cBlock++){
            NxsCharactersBlock* charBlock = nexusReader.GetCharactersBlock(taxaBlock, cBlock);
            std::string charBlockTitle = charBlock->GetTitle();
            stateSpace = 61;
            readCodonData(charBlock);
        }
        
    }
}

/**
 * @brief Construct a new alignment from a siteMatrix. This kind of matrix consits of entries
 * with values 0-60 to indicate the codon at that taxon/site pair. The matrix should be indexed
 * so that matrix[t*numChar + c] gives the character code for taxon t and site c
 */
Alignment::Alignment(int* siteMatrix, int nC, int nT) : numChar(nC), numTaxa(nT) {

    matrix = new std::bitset<61>*[numTaxa];

    matrix[0] = new std::bitset<61>[numTaxa*numChar];
    for(int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1] + numChar;
    for(int i = 0; i < numTaxa; i++)
        for(int j = 0; j < numChar; j++)
            matrix[i][j] = std::bitset<61>{};

    for(int i = 0; i < numTaxa; i++){
        for(int j = 0; j < numChar; j++){
            int charType = siteMatrix[j*numTaxa +i] % 61;
            if(charType < 0 || charType >= 61)
                Msg::error("Invalid character state from simulation!");
            matrix[i][j][charType] = 1;
        }
    }

    std::cout << "Initialized alignment from simulation" << std::endl;
}

/**
 * @brief Destructor
 * 
 */
Alignment::~Alignment(){
    delete [] matrix[0];
    delete [] matrix;
}

/**
 * @brief Takes in a NEXUS character block from the NEXUS class library and constructs a codon matrix. Each
 * bit of each entry of the ULL encodes whether or not that site can emit that character. This allows us to consider
 * ambiguous sites in an easy way without collecting a list of what ambiguous sites can represent. 
 */
void Alignment::readCodonData(NxsCharactersBlock* charBlock){
    numTaxa = charBlock->GetNumActiveTaxa();
    numChar = charBlock->GetNumActiveChar();
    if(numChar % 3 != 0)
        Msg::error("The number of characters must be a multiple of three");
    numChar = (int)numChar/3;
    CodonRecodingStruct codonDict = getCodonRecodingStruct(NxsGeneticCodesEnum::NXS_GCODE_STANDARD); //Hard code this for now

    matrix = new std::bitset<61>*[numTaxa];

    matrix[0] = new std::bitset<61>[numTaxa*numChar];
    for(int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1] + numChar;
    for(int i = 0; i < numTaxa; i++)
        for(int j = 0; j < numChar; j++)
            matrix[i][j] = std::bitset<61>{};

    char* state = new char[4];

    for(int i = 0; i < numTaxa; i++){
        taxaNames.push_back(charBlock->GetTaxonLabel(i));
        for(int j = 0; j < numChar; j++){
            state[0] = charBlock->GetState(i, j*3);
            state[1] = charBlock->GetState(i, j*3 + 1);
            state[2] = charBlock->GetState(i, j*3 + 2);
            state[3] = '\0';
            for(int k = 0, len = 61; k < len; k++){
                if(codonDict.codonStrings[k] == state){
                    matrix[i][j][k] = 1;
                    break;
                }

                if(k == len-1) { // We did not find something here - anything could have emitted it
                    matrix[i][j].set();
                }
            }
        }
    }
    delete [] state;
}