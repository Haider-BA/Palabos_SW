

#ifndef DYNAMICSPROCESSOR3DINC_HH_INCLUDED
#define DYNAMICSPROCESSOR3DINC_HH_INCLUDED

#include "basicDynamics_/dynamicsProcessor3DInC.h"
extern "C"{
    #include "basicDynamics_/dynamicsProcessor3DCImpl.h"
}

using namespace plb;
/* ************* Class ExternalRhoJcollideAndStream3D ******************* */

template<typename T, template<typename U> class Descriptor>
ExternalRhoJcollideAndStream3DInC<T,Descriptor>::ExternalRhoJcollideAndStream3DInC(DynamicsInC& backGroundDym_)
:backGroundDym(backGroundDym_)
{
    double refer;
    T TemplpateType;
    PLB_PRECONDITION(typeid(TemplpateType)==typeid(refer));
}


template<typename T, template<typename U> class Descriptor>
void ExternalRhoJcollideAndStream3DInC<T,Descriptor>::collide (
        BlockLattice3D<T,Descriptor>& lattice, Box3D const& domain,
        ScalarField3D<T> const& rhoBarField, Dot3D const& offset1,
        TensorField3D<T,3> const& jField, Dot3D const& offset2, BlockStatisticsInC& stat )
{
    void* cellGroup= &(lattice.get(0, 0, 0));
    const void* rhoBarGroup= &(rhoBarField.get(0, 0, 0));
    const void* jGroup= &(jField.get(0, 0, 0));
    ///Relative position of the cell populations array
    int cellPopPos=int((&((lattice.get(0, 0, 0))[0]))-
                       (double*)(&(lattice.get(0, 0, 0))));
    int cellLength = sizeof(lattice.get(0, 0, 0));
    int basicLength = sizeof(double);
    int fieldsDomain[]={lattice.getNx(), lattice.getNy(), lattice.getNz(),
                        rhoBarField.getNx(), rhoBarField.getNy(), rhoBarField.getNz(),
                        jField.getNx(), jField.getNy(), jField.getNz()};
    int domainCoords[]={domain.x0, domain.x1,
                        domain.y0, domain.y1,
                        domain.z0, domain.z1};
    int offset1Array[] = {offset1.x, offset1.y, offset1.z};
    int offset2Array[] = {offset2.x, offset2.y, offset2.z};

    ExternalRhoJcollideAndStream3DInC_collide(&backGroundDym, cellGroup, cellLength,
                                              rhoBarGroup, jGroup, cellPopPos, fieldsDomain,
                                              basicLength, offset1Array, offset2Array, domainCoords,
                                              &stat);
}

template<typename T, template<typename U> class Descriptor>
void ExternalRhoJcollideAndStream3DInC<T,Descriptor>::bulkCollideAndStream (
        BlockLattice3D<T,Descriptor>& lattice, Box3D const& domain,
        ScalarField3D<T> const& rhoBarField, Dot3D const& offset1,
        TensorField3D<T,3> const& jField, Dot3D const& offset2, BlockStatisticsInC& stat )
{
    void* cellGroup= &(lattice.get(0, 0, 0));
    const void* rhoBarGroup= &(rhoBarField.get(0, 0, 0));
    const void* jGroup= &(jField.get(0, 0, 0));
    ///Relative position of the cell populations array
    int cellPopPos=int((&((lattice.get(0, 0, 0))[0]))-
                       (double*)(&(lattice.get(0, 0, 0))));
    int cellLength = sizeof(lattice.get(0, 0, 0));
    int basicLength = sizeof(double);
    int fieldsDomain[]={lattice.getNx(), lattice.getNy(), lattice.getNz(),
                        rhoBarField.getNx(), rhoBarField.getNy(), rhoBarField.getNz(),
                        jField.getNx(), jField.getNy(), jField.getNz()};
    int domainCoords[]={domain.x0, domain.x1,
                        domain.y0, domain.y1,
                        domain.z0, domain.z1};
    int offset1Array[] = {offset1.x, offset1.y, offset1.z};
    int offset2Array[] = {offset2.x, offset2.y, offset2.z};

    //pcout<<"**********flag2.2.1**********"<<std::endl;

    ExternalRhoJcollideAndStream3DInC_bulkCollideAndStream(&backGroundDym, cellGroup, cellLength,
                                              rhoBarGroup, jGroup, cellPopPos, fieldsDomain,
                                              basicLength, offset1Array, offset2Array, domainCoords,
                                              &stat);

}


template<typename T, template<typename U> class Descriptor>
void ExternalRhoJcollideAndStream3DInC<T,Descriptor>::boundaryStream (
        BlockLattice3D<T,Descriptor>& lattice,
        Box3D const& bound, Box3D const& domain )
{
    // Make sure domain is contained within bound
    PLB_PRECONDITION( contained(domain, bound) );

    void* cellGroup= &(lattice.get(0, 0, 0));
    ///Relative position of the cell populations array
    int cellPopPos=int((&((lattice.get(0, 0, 0))[0]))-
                       (double*)(&(lattice.get(0, 0, 0))));
    int cellLength = sizeof(lattice.get(0, 0, 0));
    int basicLength = sizeof(double);
    int fieldsDomain[]={lattice.getNx(), lattice.getNy(), lattice.getNz()};
    int domainCoords[]={domain.x0, domain.x1,
                        domain.y0, domain.y1,
                        domain.z0, domain.z1};
    int boundCoords[]={bound.x0, bound.x1,
                       bound.y0, bound.y1,
                       bound.z0, bound.z1};

    ExternalRhoJcollideAndStream3DInC_boundaryStream(&backGroundDym, cellGroup, cellLength,
                                          cellPopPos, fieldsDomain, basicLength, domainCoords,
                                          boundCoords);


}

template<typename T, template<typename U> class Descriptor>
void ExternalRhoJcollideAndStream3DInC<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    global::timer("collideAndStream").start();
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]);
    ScalarField3D<T> const& rhoBarField =
        dynamic_cast<ScalarField3D<T> const&>(*atomicBlocks[1]);
    TensorField3D<T,3> const& jField =
        dynamic_cast<TensorField3D<T,3> const&>(*atomicBlocks[2]);

    struct BlockStatisticsInC stat;//BlockStatistics& stat = lattice.getInternalStatistics();

    static const plint vicinity = Descriptor<T>::vicinity;
    Box3D extDomain(domain.enlarge(vicinity));

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);

    global::profiler().start("collStream");
    global::profiler().increment("collStreamCells", extDomain.nCells());

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    //pcout<<"**********flag2.1**********"<<std::endl;
    collide(lattice,
            Box3D(extDomain.x0,extDomain.x0+vicinity-1,
                  extDomain.y0,extDomain.y1,
                  extDomain.z0,extDomain.z1),
            rhoBarField, offset1, jField, offset2, stat);
    collide(lattice,
            Box3D(extDomain.x1-vicinity+1,extDomain.x1,
                  extDomain.y0,extDomain.y1,
                  extDomain.z0,extDomain.z1),
            rhoBarField, offset1, jField, offset2, stat);
    collide(lattice,
            Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                  extDomain.y0,extDomain.y0+vicinity-1,
                  extDomain.z0,extDomain.z1),
            rhoBarField, offset1, jField, offset2, stat);
    collide(lattice,
            Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                  extDomain.y1-vicinity+1,extDomain.y1,
                  extDomain.z0,extDomain.z1),
            rhoBarField, offset1, jField, offset2, stat);
    collide(lattice,
            Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                  extDomain.y0+vicinity,extDomain.y1-vicinity,
                  extDomain.z0,extDomain.z0+vicinity-1),
            rhoBarField, offset1, jField, offset2, stat);
    collide(lattice,
            Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                  extDomain.y0+vicinity,extDomain.y1-vicinity,
                  extDomain.z1-vicinity+1,extDomain.z1),
            rhoBarField, offset1, jField, offset2, stat);

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    //pcout<<"**********flag2.2**********"<<std::endl;
    bulkCollideAndStream(lattice,
                         Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                               extDomain.y0+vicinity,extDomain.y1-vicinity,
                               extDomain.z0+vicinity,extDomain.z1-vicinity),
                         rhoBarField, offset1, jField, offset2, stat);

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    //pcout<<"**********flag2.3**********"<<std::endl;
    boundaryStream(lattice, extDomain, Box3D(extDomain.x0,extDomain.x0+vicinity-1,
                                             extDomain.y0,extDomain.y1,
                                             extDomain.z0,extDomain.z1) );
    //pcout<<"**********flag2.4**********"<<std::endl;
    boundaryStream(lattice, extDomain, Box3D(extDomain.x1-vicinity+1,extDomain.x1,
                                             extDomain.y0,extDomain.y1,
                                             extDomain.z0,extDomain.z1) );
    boundaryStream(lattice, extDomain, Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                                             extDomain.y0,extDomain.y0+vicinity-1,
                                             extDomain.z0,extDomain.z1) );
    boundaryStream(lattice, extDomain, Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                                             extDomain.y1-vicinity+1,extDomain.y1,
                                             extDomain.z0,extDomain.z1) );
    boundaryStream(lattice, extDomain, Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                                             extDomain.y0+vicinity,extDomain.y1-vicinity,
                                             extDomain.z0,extDomain.z0+vicinity-1) );
    boundaryStream(lattice, extDomain, Box3D(extDomain.x0+vicinity,extDomain.x1-vicinity,
                                             extDomain.y0+vicinity,extDomain.y1-vicinity,
                                             extDomain.z1-vicinity+1,extDomain.z1) );
    //pcout<<"**********flag2.5**********"<<std::endl;
    global::profiler().stop("collStream");
    global::timer("collideAndStream").stop();
    //pcout<<"**********flag2.6**********"<<std::endl;
}


template<typename T, template<typename U> class Descriptor>
ExternalRhoJcollideAndStream3DInC<T,Descriptor>*
    ExternalRhoJcollideAndStream3DInC<T,Descriptor>::clone() const
{
    return new ExternalRhoJcollideAndStream3DInC<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ExternalRhoJcollideAndStream3DInC<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}


#endif // DYNAMICSPROCESSOR3DINC_HH_INCLUDED
