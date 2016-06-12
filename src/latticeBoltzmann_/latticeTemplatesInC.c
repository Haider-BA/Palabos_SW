

#include "latticeBoltzmann_/latticeTemplatesInC.h"
#include "stdio.h"



void swapAndStream3DInC_D3Q19(struct DescriptorInC* descriptor_, void* cellGroup, int cellLength,
                                         int cellPopPos, int fieldsDomain[3], int basicLength,
                                         int iX, int iY, int iZ)
{
    double fTmp;
    //printf("**********flag2.2.1.1.1.1**********\n");
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX-1, iY,   iZ,   1, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX,   iY-1, iZ,   2, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX,   iY  , iZ-1, 3, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX-1, iY-1, iZ,   4, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX-1, iY+1, iZ,   5, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX-1, iY  , iZ-1, 6, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX-1, iY  , iZ+1, 7, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX  , iY-1, iZ-1, 8, &fTmp);
    swapAndStreamCellInC_D3Q19(cellGroup, cellLength, cellPopPos, fieldsDomain, basicLength,
                               iX, iY, iZ, iX  , iY-1, iZ+1, 9, &fTmp);
}

void swapAndStreamCellInC_D3Q19(
      void* cellGroup, int cellLength, int cellPopPos, int fieldsDomain[3], int basicLength,
      int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, double* fTmp )
{
    //printf("**********flag2.2.1.1.1.1.1**********\n");
    //printf("%d", fieldsDomain[0]);
    //printf("\n**********flag2.2.1.1.1.1.2**********");
    *fTmp   =  *((double*) (cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +iY*cellLength*fieldsDomain[2]
                                    +iZ*cellLength)+cellPopPos+iPop*basicLength));





    *((double*) (cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +iY*cellLength*fieldsDomain[2]
                                    +iZ*cellLength)+cellPopPos+iPop*basicLength))
    =
    *((double*) (cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +iY*cellLength*fieldsDomain[2]
                                    +iZ*cellLength)+cellPopPos+(iPop+9)*basicLength));
    //printf("**********flag2.2.1.1.1.1.3**********");

    *((double*) (cellGroup+(iX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +iY*cellLength*fieldsDomain[2]
                                    +iZ*cellLength)+cellPopPos+(iPop+9)*basicLength))
    =
    *((double*) (cellGroup+(nX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +nY*cellLength*fieldsDomain[2]
                                    +nZ*cellLength)+cellPopPos+iPop*basicLength));
    //printf("**********flag2.2.1.1.1.1.4**********");

    *((double*) (cellGroup+(nX*cellLength*fieldsDomain[1]*fieldsDomain[2]
                                    +nY*cellLength*fieldsDomain[2]
                                    +nZ*cellLength)+cellPopPos+iPop*basicLength))
    =
    *fTmp;
    //printf("**********flag2.2.1.1.1.1.5**********");
}


void swapAndStream3DInC(struct DescriptorInC* descriptor_, void* cellGroup, int cellLength,  int cellPopPos,
                   int fieldsDomain[3], int basicLength, int iX, int iY, int iZ)
{
    if((*descriptor_).d==3)
        switch((*descriptor_).q)
    {
        case 13:
            {
                break;
            }
        case 15:
            {
                break;
            }
        case 19:
            {
                //printf("**********flag2.2.1.1.1**********");
                swapAndStream3DInC_D3Q19(descriptor_, cellGroup, cellLength,
                                         cellPopPos, fieldsDomain, basicLength,
                                         iX, iY, iZ);
                break;
            }
        case 27:
            {
                break;
            }
        case 121:
            {
                break;
            }
        default:;
    }
    else if((*descriptor_).d==2)
        switch((*descriptor_).q)
    {
        case 5:
            {
                break;
            }
        case 9:
            {
                break;
            }
        case 37:
            {
                break;
            }
        default:
            {
                break;
            }
    }
    else ;///Output An error!!
}

