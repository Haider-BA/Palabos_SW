#ifndef DESCRIPTORINC_H_INCLUDED
#define DESCRIPTORINC_H_INCLUDED


///------------------------ The Descriptor structure --------------------
struct DescriptorInC
{
    ///the decriptor constants
    const int d;

    const int q;

    const char 	name [50];

    const double invD;

    const int 	vicinity;

    const double 	cs2 ;//= (double)1 / (double)3;

    const double	invCs2 ;//= (double)3;

    /// The round-off policy funciton pointers
    int 	(*SkordosFactor) ();

    double 	(*rhoBar) (double);

    double 	(*fullRho) (double);

    double 	(*invRho) (double);

    double 	(*rhoMinus1) (double);

    ///The external field constants
    const int 	numScalars;

    const int 	numSpecies;

    const int 	forceBeginsAt;

    const int 	sizeOfForce;

    ///Since the arrays in C can't be determined dynamically, the d and q is
    ///replaced with constants 3 and 200 to retain enough space for any
    ///descriptor.
    const int   cNormSqr [200];

    const double        t [200];

    const int  c [200][3];

};

//--------------------- Descriptor constructors ---------------------------

struct DescriptorInC constructD3Q19Descriptor();



#endif // DESCRIPTORINC_H_INCLUDED
