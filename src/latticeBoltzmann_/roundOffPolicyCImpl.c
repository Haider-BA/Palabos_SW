


#include "latticeBoltzmann_/roundOffPolicyCImpl.h"



//------------------- DefaultRoundOffPolicy -----------------------
int SkordosFactor_DefaultRoundOffPolicy ()
{
    return 1;
}

double rhoBar_DefaultRoundOffPolicy (double rho)
{
    return rho-(double)1;
}

double fullRho_DefaultRoundOffPolicy (double rhoBar)
{
    return rhoBar+(double)1;
}

double invRho_DefaultRoundOffPolicy (double rhoBar)
{
    return (double)1/(rhoBar+(double)1);
}

double rhoMinus1_DefaultRoundOffPolicy (double rhoBar)
{
    return rhoBar;
}




