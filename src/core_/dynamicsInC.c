

#include "core_/dynamicsInC.h"
#include "core_/dynamicsCImpl.h"
#include "core_/blockStatisticsInC.h"

//----------------------- Dynamics Constructor --------------------------
void constructDynamics(struct DynamicsInC* dynamicsCase, struct DescriptorInC* _descriptor)
{
    (*dynamicsCase).descriptor_=_descriptor;
    (*dynamicsCase).collideExternal=collideExternal_Dynanics;
    (*dynamicsCase).computeEquilibria=computeEquilibria_Dynamics;
}

void constructBasicBulkDynamics(struct DynamicsInC* dynamicsCase, struct DescriptorInC* _descriptor, double omega_)
{
    constructDynamics(dynamicsCase, _descriptor);
    (*dynamicsCase).computeRhoBarJ=computeRhoBarJ_BasicBulkDynamics;
    (*dynamicsCase).omega=omega_;
}
