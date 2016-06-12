#ifndef DYNAMICSCIMPL_H_INCLUDED
#define DYNAMICSCIMPL_H_INCLUDED

#include "core_/globalDefsInC.h"
#include "core_/blockStatisticsInC.h"
#include "core_/dynamicsInC.h"

//----------------------- Dynamics -------------------------------

void collideExternal_Dynanics(struct DynamicsInC* _this, struct UniformSequence cellPopulation, double rhoBar,
                              struct UniformSequence_const j, struct BlockStatisticsInC* statistics, double thetaBar);

void computeEquilibria_Dynamics (struct DynamicsInC*, struct UniformSequence, double, struct UniformSequence_const, double, double);

//--------------------------- BasicBulkDynamics --------------------------------
void computeRhoBarJ_BasicBulkDynamics (struct DynamicsInC*, struct UniformSequence_const, double*, struct UniformSequence);


#endif // DYNAMICSCIMPL_H_INCLUDED
