

#include "latticeBoltzmann_/descriptorInC.h"
#include "latticeBoltzmann_/roundOffPolicyCImpl.h"

//--------------------- Descriptor constructors ---------------------------

struct DescriptorInC constructD3Q19Descriptor()
{
    struct DescriptorInC descriptorCase=
    {
     ///the decriptor constants
    //enum  	{ d = 3, q = 19 };
    3,

    19,

    "D3Q19",

    (double)1 / (double) 3,

    1,

    (double)1 / (double)3,

    (double)3,

    /// The round-off policy funciton pointers
    SkordosFactor_DefaultRoundOffPolicy,

    rhoBar_DefaultRoundOffPolicy,

    fullRho_DefaultRoundOffPolicy,

    invRho_DefaultRoundOffPolicy,

    rhoMinus1_DefaultRoundOffPolicy,

    ///The external field constants
    0,

    0,

    0,

    0,

    { 0, 1, 1, 1, 2, 2, 2, 2, 2, 2,
             1, 1, 1, 2, 2, 2, 2, 2, 2, },

    {
        1.0/3.0,
        1.0/18.0, 1.0/18.0, 1.0/18.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0,
        1.0/18.0, 1.0/18.0, 1.0/18.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0,
    },

    {
        { 0, 0, 0},
        {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
        {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
        {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
        { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
        { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
        { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
    }

    };
    return descriptorCase;
}
