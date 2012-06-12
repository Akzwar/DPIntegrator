#include <stdio.h>
#include "dpintegrator.h"

int main(int argc, char** argv)
{
	DotMassesModel* model = new DotMassesModel;
	DPIntegrator integr( model, 0, 100, 10e-8 );
	while( integr.getT() < integr.getTk() )
	{
	//	printf("%f\n",(double)integr.getStep());	
	integr.PhaseVect().print();
		integr.NextStep();	
	}
	return 0;
}
