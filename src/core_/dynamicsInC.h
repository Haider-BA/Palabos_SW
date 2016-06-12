#ifndef DYNAMICSINC_H_INCLUDED
#define DYNAMICSINC_H_INCLUDED



#include "core_/globalDefsInC.h"
#include "latticeBoltzmann_/descriptorInC.h"

struct DynamicsInC
{
    /// Pointer to the function that implement the member "collide" in Class Dynamics
    void (*computeRhoBarJ) (struct DynamicsInC*, struct UniformSequence_const, double* , struct UniformSequence);

    /// Pointer to the function that implement the member "collide" in Class Dynamics
    void (*collide) (struct DynamicsInC* , struct UniformSequence, struct BlockStatisticsInC*);

    /// Pointer to the function that implement the member "collideExternal" in Class Dynamics
    void (*collideExternal) (struct DynamicsInC* , struct UniformSequence, double, struct UniformSequence_const, struct BlockStatisticsInC* ,double);

    /// Pointer to the function that implement the member "computeEquilibrium" in Class Dynamics
    double (*computeEquilibrium) (struct DynamicsInC *, int, double, struct UniformSequence_const, double, double);

    /// Pointer to the function that implement the member "computeEquilibria" in Class Dynamics
    void (*computeEquilibria) (struct DynamicsInC *, struct UniformSequence, double, struct UniformSequence_const, double, double);

    /// Pointer to background descriptor
    struct DescriptorInC* descriptor_;

    /// Most dynamics use the relaxation parameter omega, so store it in the Base structure
    double omega;
};


//----------------------- Dynamics Constructor --------------------------
void constructDynamics(struct DynamicsInC*, struct DescriptorInC*);

void constructBasicBulkDynamics(struct DynamicsInC* , struct DescriptorInC *, double);




#endif // DYNAMICSINC_H_INCLUDED
