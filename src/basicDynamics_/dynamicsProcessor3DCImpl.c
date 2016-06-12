

#include "core_/globalDefsInC.h"
#include "core_/cellCImpl.h"
#include "latticeBoltzmann_/latticeTemplatesInC.h"
#include "basicDynamics_/dynamicsProcessor3DCImpl.h"

#include "stdio.h"

void ExternalRhoJcollideAndStream3DInC_collide(struct DynamicsInC* backGroundDym, void* cellGroup, int cellLength,
                                          const void* rhoBarGroup, const void* jGroup, int cellPopPos, int fieldsDomain[9],
                                          int basicLength, int offset1Array[3], int offset2Array[3], int domainCoords[6],
                                          struct BlockStatisticsInC* stat)
{
    for (int iX=domainCoords[0]; iX<=domainCoords[1]; ++iX) {
        for (int iY=domainCoords[2]; iY<=domainCoords[3]; ++iY) {
            for (int iZ=domainCoords[4]; iZ<=domainCoords[5]; ++iZ) {
                struct UniformSequence cellPopulations = {(*(*backGroundDym).descriptor_).q, basicLength,
                                                   cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                                              +iY*cellLength*fieldsDomain[2]
                                                              +iZ*cellLength)+cellPopPos};
                double* rhoBar =
                    (double*) (rhoBarGroup+((iX+offset1Array[0])*basicLength*fieldsDomain[4]*fieldsDomain[5]
                                          +(iY+offset1Array[1])*basicLength*fieldsDomain[5]
                                          +(iZ+offset1Array[2])*basicLength));

                int jds =  (*(*backGroundDym).descriptor_).d;
                struct UniformSequence_const jPops =
                    {jds, basicLength,
                     jGroup+((iX+offset2Array[0])*basicLength*jds*fieldsDomain[7]*fieldsDomain[8]
                     +(iY+offset2Array[1])*basicLength*jds*fieldsDomain[8]
                     +(iZ+offset2Array[2])*basicLength*jds)};

                (*backGroundDym).collideExternal(backGroundDym, cellPopulations, *rhoBar, jPops, stat, 0.0);
                cellRevert(cellPopulations);
            }
        }
    }
}

void ExternalRhoJcollideAndStream3DInC_bulkCollideAndStream(struct DynamicsInC* backGroundDym, void* cellGroup, int cellLength,
                                          const void* rhoBarGroup, const void* jGroup, int cellPopPos, int fieldsDomain[9],
                                          int basicLength, int offset1Array[3], int offset2Array[3], int domainCoords[6],
                                          struct BlockStatisticsInC* stat)
{
    int latticeDomain[]={fieldsDomain[0],fieldsDomain[1],fieldsDomain[2]};
    //printf("\nThe lattice domain is %dx%dx%d\n", latticeDomain[0], latticeDomain[1], latticeDomain[2]);
    for (int iX=domainCoords[0]; iX<=domainCoords[1]; ++iX) {
        for (int iY=domainCoords[2]; iY<=domainCoords[3]; ++iY) {
            for (int iZ=domainCoords[4]; iZ<=domainCoords[5]; ++iZ) {
                struct UniformSequence cellPopulations = {(*(*backGroundDym).descriptor_).q, basicLength,
                                                   cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                                              +iY*cellLength*fieldsDomain[2]
                                                              +iZ*cellLength)+cellPopPos};
                double* rhoBar =
                    (double*) (rhoBarGroup+((iX+offset1Array[0])*basicLength*fieldsDomain[4]*fieldsDomain[5]
                                          +(iY+offset1Array[1])*basicLength*fieldsDomain[5]
                                          +(iZ+offset1Array[2])*basicLength));

                int jds =  (*(*backGroundDym).descriptor_).d;
                struct UniformSequence_const jPops =
                    {jds, basicLength,
                     jGroup+((iX+offset2Array[0])*basicLength*jds*fieldsDomain[7]*fieldsDomain[8]
                     +(iY+offset2Array[1])*basicLength*jds*fieldsDomain[8]
                     +(iZ+offset2Array[2])*basicLength*jds)};

                (*backGroundDym).collideExternal(backGroundDym, cellPopulations, *rhoBar, jPops, stat, 0.0);
                //printf("**********flag2.2.1.1**********");
                swapAndStream3DInC((*backGroundDym).descriptor_, cellGroup, cellLength, cellPopPos, latticeDomain, basicLength, iX, iY, iZ);
            }
        }
    }
}




void ExternalRhoJcollideAndStream3DInC_boundaryStream(struct DynamicsInC* backGroundDym, void* cellGroup, int cellLength,
                                          int cellPopPos, int fieldsDomain[3], int basicLength, int domainCoords[6],
                                          int boundCoords[6])
{
    double swapTmp;
    for (int iX=domainCoords[0]; iX<=domainCoords[1]; ++iX) {
        for (int iY=domainCoords[2]; iY<=domainCoords[3]; ++iY) {
            for (int iZ=domainCoords[4]; iZ<=domainCoords[5]; ++iZ) {
                for (int iPop=1; iPop<=(*(*backGroundDym).descriptor_).q/2; ++iPop) {
                    int nextX = iX + (*(*backGroundDym).descriptor_).c[iPop][0];
                    int nextY = iY + (*(*backGroundDym).descriptor_).c[iPop][1];
                    int nextZ = iZ + (*(*backGroundDym).descriptor_).c[iPop][2];
                    if ( nextX>=boundCoords[0] && nextX<=boundCoords[1] && nextY>=boundCoords[2] &&
                         nextY<=boundCoords[3] && nextZ>=boundCoords[4] && nextZ<=boundCoords[5] )
                    {
                        //swap(
                        swapTmp=
                             *((double*)(cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                        +iY*cellLength*fieldsDomain[2]
                                        +iZ*cellLength)+cellPopPos
                                        +basicLength*(iPop+(*(*backGroundDym).descriptor_).q/2)));


                        *((double*)(cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +iY*cellLength*fieldsDomain[2]
                                    +iZ*cellLength)+cellPopPos
                                    +basicLength*(iPop+(*(*backGroundDym).descriptor_).q/2)))
                                    =
                            *((double*)(cellGroup+(nextX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                        +nextY*cellLength*fieldsDomain[2]
                                        +nextZ*cellLength)+cellPopPos
                                        +basicLength*iPop));

                        *((double*)(cellGroup+(nextX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +nextY*cellLength*fieldsDomain[2]
                                    +nextZ*cellLength)+cellPopPos
                                    +basicLength*iPop))   =    swapTmp;
                         //   );
                    }
                }
            }
        }
    }
}
