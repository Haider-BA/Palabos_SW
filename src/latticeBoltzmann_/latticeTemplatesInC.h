#ifndef LATTICETEMPLATESINC_H_INCLUDED
#define LATTICETEMPLATESINC_H_INCLUDED

#include "latticeBoltzmann_/descriptorInC.h"

void swapAndStream3DInC_D3Q19(struct DescriptorInC* descriptor_, void* cellGroup, int cellLength,
                                         int cellPopPos, int fieldsDomain[3], int basicLength,
                                         int iX, int iY, int iZ);

void swapAndStream3DInC(struct DescriptorInC* descriptor_, void* cellGroup, int cellLength,  int cellPopPos,
                   int fieldsDomain[3], int basicLength, int iX, int iY, int iZ);


#endif // LATTICETEMPLATESINC_H_INCLUDED
