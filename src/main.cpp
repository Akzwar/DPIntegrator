#include <stdio.h>
#include "dpintegrator.h"

int main(int argc, char** argv)
{
	NormalSphereFuncModel* model = new NormalSphereFuncModel;
	DPIntegrator integr( model, 0, 100, 10e-18 );
	//while( integr.getT() < integr.getTk() )
	for(int i=0;i<200; i++)
	{
		integr.PhaseVect().print();
		integr.NextStep();	
	}
	//printf("%f",(double)LegandrSin(4,0,1));
	return 0;
}
