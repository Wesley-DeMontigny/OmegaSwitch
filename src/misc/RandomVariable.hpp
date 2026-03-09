#ifndef RANDOM_VARIABLE_HPP
#define RANDOM_VARIABLE_HPP

#include <cstdint>
#include <vector>

enum {
    N = 624,
    M = 397,
    R = 31,
    A = 0x9908B0DF,
    F = 1812433253,
    U = 11,
    S = 7,
    B = 0x9D2C5680,
    T = 15,
    C = 0xEFC60000,
    L = 18,
    MASK_LOWER = (1ull << R) - 1,
    MASK_UPPER = (1ull << R)
};

/**
 * @brief Mersenne-Twister random number generator
 */
class RandomVariable {

    public:
        static RandomVariable&  randomVariableInstance(void){               // Returns reference to RandomVariable object
                                static RandomVariable singleRandomVariable; // singleRandomVariable is initialized once
                                return singleRandomVariable;
        }

        double                  uniformRv(void);                            // Draw from a uniform distribution on the unit interval

    private:
                                RandomVariable(void);                       // Default constructor 
                                RandomVariable(RandomVariable& r);          // Copy constructor
                                RandomVariable(uint32_t seed);              // Initialize with a particular seed
        RandomVariable&         operator=(const RandomVariable&);           // Assignment operator
        uint32_t                extractU32(void);                           // Extracts a random integer from the generator's internal state
        void                    initialize(uint32_t seed);                  // Initialize with a particular seed
        void                    twist(void);                                // Performs the "twist" transformation on the internal state array
        uint16_t                index;                                      // Current position within the internal state array
        uint32_t                mt[N];                                      // Holds the state for producing pseudorandom values
};

#endif
