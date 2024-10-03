#include "CodonMultiMatrix.hpp"
#include "core/Matrix.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"

CodonMultiMatrix::CodonMultiMatrix(DirichletProcessPrior* o, BasicParameter<double>* k, std::vector<BasicParameter<double>*> pi) : 
                                   currentQMatrix(61, 61, 0.0), dNdS(o), 
                                   oldQMatrix(61, 61, 0.0), stationaryDist(pi), transitionTransversionRatio(k) {
    
    std::vector<int> aaMap = {8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 19, 19, 15, 15, 15, 15, 1, 18, 1, 9, 4, 9, 4};  
    std::vector<char*> codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

    double ttR = transitionTransversionRatio->getValue();

    // Because of the complicated nature of this matrix, we need to classify each of the positions in the matrix;
    for(int i = 0; i < 61; i++){
        for(int j = i + 1; j < 61; j++){
            int mismatch = 0;
            for(int k = 0; k < 3; k++){
                if(codons[i][k] != codons[j][k]){
                    mismatch++;
                    if(mismatch > 1){
                        break;
                    }
                    if((codons[i][k] == 'A' && codons[j][k] == 'G') || (codons[i][k] == 'G' && codons[j][k] == 'A') || 
                       (codons[i][k] == 'T' && codons[j][k] == 'C') || (codons[i][k] == 'C' && codons[j][k] == 'T'))
                        transition.insert(std::make_pair(i, j));
                }
            }
            if(mismatch == 1){
                valid.insert(std::make_pair(i, j));
                if(aaMap[i] != aaMap[j])
                    nonsynonymous.insert(std::make_pair(i, j));
            }
        }
        returnStationary.push_back(stationaryDist[i]->getValue());
    }

    for(auto coord : valid){
        currentQMatrix(coord.first, coord.second) = returnStationary[coord.first];
        currentQMatrix(coord.second, coord.first) = returnStationary[coord.second];
    }
    for(auto coord : transition){
        currentQMatrix(coord.first, coord.second) *= ttR;
        currentQMatrix(coord.second, coord.first) *= ttR;
    }
    
    oldQMatrix = currentQMatrix;
    dirty();
}

void CodonMultiMatrix::accept() {
    dNdS->accept();
    dNdS->clean();
    transitionTransversionRatio->accept();
    transitionTransversionRatio->clean();

    for(BasicParameter<double>* f : stationaryDist){
        f->accept();
        f->clean();
    }

    oldQMatrix = currentQMatrix;
}

void CodonMultiMatrix::reject() {
    dNdS->reject();
    dNdS->clean();
    transitionTransversionRatio->reject();
    transitionTransversionRatio->clean();

    for(BasicParameter<double>* f : stationaryDist){
        f->reject();
        f->clean();
    }

    currentQMatrix = oldQMatrix;
}

void CodonMultiMatrix::regenerate() {
    dNdS->regenerate();
    transitionTransversionRatio->regenerate();

    if(dNdS->isDirty() || transitionTransversionRatio->isDirty())
        this->dirty();

    bool dirtyStationary =false;
    for(BasicParameter<double>* f : stationaryDist){
        f->regenerate();

        if(f->isDirty()){
            this->dirty();
            dirtyStationary = true;
        }
    }

    if(this->isDirty()){
        double ttR = transitionTransversionRatio->getValue();
        if(dirtyStationary){
            returnStationary.clear();
            for(int i = 0; i < 61; i++){
                returnStationary.push_back(stationaryDist[i]->getValue());
            }
        }

        for(auto coord : valid){
            currentQMatrix(coord.first, coord.second) = returnStationary[coord.first];
            currentQMatrix(coord.second, coord.first) = returnStationary[coord.second];
        }
        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) *= ttR;
            currentQMatrix(coord.second, coord.first) *= ttR;
        }
    }
}

std::vector<double> CodonMultiMatrix::stationary(){
    return returnStationary;
}

Matrix<double> CodonMultiMatrix::Q(int i) {
    Matrix<double> returnMatrix = currentQMatrix;
    double ddR = dNdS->getCategoryValue(i);

    for(auto coord : nonsynonymous){
        returnMatrix(coord.first, coord.second) *= ddR;
        returnMatrix(coord.second, coord.first) *= ddR; 
    }

    for(int i = 0; i < 61; i++){
        double total = 0.0;
        for(int j = 0; j < 61; j++){
            if(j != i) total += returnMatrix(i , j);
        }
        returnMatrix(i , i) = total * -1;
    }
    return returnMatrix;
}
