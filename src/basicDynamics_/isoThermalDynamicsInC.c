

#include "basicDynamics_/isoThermalDynamicsInC.h"
#include "basicDynamics_/isoThermalDynamicsCImpl.h"


void constructBGKDynamics(struct DynamicsInC* dynamicsCase, struct DescriptorInC* _descriptor, double omega_)
{
    constructBasicBulkDynamics(dynamicsCase , _descriptor, omega_);
    (*dynamicsCase).collide = collide_BGKDynamics;
    (*dynamicsCase).collideExternal = collideExternal_BGKDynamics;
    (*dynamicsCase).computeEquilibrium = computeEquilibrium_BGKDynamics;
}
