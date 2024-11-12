#include "Alignment.hpp"
#include "Msg.hpp"
#include "ncl/nxscharactersblock.cpp"

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
            for(int i = 0; i < stateSpace; i++)
                frequencies.push_back(0);
            readCodonData(charBlock);
        }
        
    }
}

Alignment::~Alignment(){
    delete [] matrix[0];
    delete [] matrix;
}

void Alignment::readCodonData(NxsCharactersBlock* charBlock){
    numTaxa = charBlock->GetNumActiveTaxa();
    numChar = charBlock->GetNumActiveChar();
    if(numChar % 3 != 0)
        Msg::error("The number of characters must be a multiple of three");
    numChar = (int)numChar/3;
    CodonRecodingStruct codonDict = getCodonRecodingStruct(NxsGeneticCodesEnum::NXS_GCODE_STANDARD); //Hard code this for now

    matrix = new unsigned long long int*[numTaxa];

    matrix[0] = new unsigned long long int[numTaxa*numChar];
    for(int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1] + numChar;
    for(int i = 0; i < numTaxa; i++)
        for(int j = 0; j < numChar; j++)
            matrix[i][j] = 0;

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
                    matrix[i][j] = 1ULL << k;
                    frequencies[k] += 1;
                    break;
                }

                if(k == len-1) { // We did not find something here...
                    matrix[i][j] = (1ULL << 61) - 1;
                    //Msg::error("Found an unsupported codon (" + std::string(state) + ") at "  + to_string(i) + "," + to_string(j));             
                }
            }
        }
    }
    delete [] state;

    double total = 0.0;
    for(double v : frequencies){
        total += v + 1;
    }

    for(double& v : frequencies){
        v = (v + 1)/total;
    }
}