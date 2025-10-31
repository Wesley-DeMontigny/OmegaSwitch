#ifndef PARAMETER_HPP
#define PARAMETER_HPP
#include <vector>
#include <string>

/**
 * @brief The base class for all parameter types. We leave everything virtual except for the
 * things related to the "dirtyness" of a parameter (which indicates to the Model if it has
 * been updated).
 */
class Parameter {
    public:
        virtual void accept()=0;                // Accept the changes that have been made to the parameter
        virtual void reject()=0;                // Reject the changes that have been made to the parameter 
        virtual double lnPrior()=0;             // Return the log prior probability for the current parameter value
        virtual void tune()=0;                  // Tune MCMC moves that operate on this parameter
        void clean() {dirtyFlag = false;}       // Mark the parameter has having not been updated
        void dirty() {dirtyFlag = true;}        // Mark the parameter as having been updated
        bool isDirty() const {return dirtyFlag;}// Check if the parameter has been updated
    private:
        bool dirtyFlag;                         // The flag indicating if this parameter has been updating
};

#endif