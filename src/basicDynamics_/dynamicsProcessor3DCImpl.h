#ifndef DYNAMICSPROCESSOR3DCIMPL_H_INCLUDED
#define DYNAMICSPROCESSOR3DCIMPL_H_INCLUDED


#include "core_/dynamicsInC.h"

struct DynamicsInC;
struct BlockStatisticsInC;

void ExternalRhoJcollideAndStream3DInC_collide(struct DynamicsInC* backGroundDym, void* cellGroup, int cellLength,
                                          const void* rhoBarGroup, const void* jGroup, int cellPopPos, int fieldsDomain[9],
                                          int basicLength, int offset1Array[3], int offset2Array[3], int domainCoords[6],
                                          struct BlockStatisticsInC* stat);

void ExternalRhoJcollideAndStream3DInC_bulkCollideAndStream(struct DynamicsInC* backGroundDym, void* cellGroup, int cellLength,
                                          const void* rhoBarGroup, const void* jGroup, int cellPopPos, int fieldsDomain[9],
                                          int basicLength, int offset1Array[3], int offset2Array[3], int domainCoords[6],
                                          struct BlockStatisticsInC* stat);

void ExternalRhoJcollideAndStream3DInC_boundaryStream(struct DynamicsInC* backGroundDym, void* cellGroup, int cellLength,
                                          int cellPopPos, int fieldsDomain[3], int basicLength, int domainCoords[6],
                                          int boundCoords[6]);


#endif // DYNAMICSPROCESSOR3DCIMPL_H_INCLUDED
