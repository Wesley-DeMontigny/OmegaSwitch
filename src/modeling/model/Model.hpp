#ifndef MODEL_HPP
#define MODEL_HPP
#include <string>
#include<tuple>
#include<vector>

/**
 * @brief 
 * 
 */
class Model{
    public:
        virtual double                  lnLikelihood() = 0;                                     //
        virtual double                  lnPrior() = 0;                                          //
        virtual std::vector<double>     getTunableParameterRecord() const = 0;                  //
        virtual std::vector<double>     getTunableParameters() const = 0;                       //
        virtual void                    accept() = 0;                                           //
        virtual void                    printAcceptanceRates() = 0;                             //
        virtual void                    printTabular(int i) = 0;                                //
        virtual void                    regenerateLikelihood() = 0;                             //
        virtual void                    reject() = 0;                                           //
        virtual void                    setTunableParameters(const std::vector<double> & v) = 0;//
        virtual void                    tuneMoves() = 0;                                        //
        virtual void                    writeLogData(int i) = 0;                                //
        virtual void                    writeLogHeaders() = 0;                                  //
};

#endif