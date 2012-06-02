#include <stdio.h>
#include "dpintegrator.h"

int main(int argc, char** argv)
{
	KeplerModel* model = new KeplerModel;
	DPIntegrator integr( model, 0, 100, 10e-8 );
	for(int i = 0; i<=100; i++)
	{	
		integr.NextStep();	
		printf("%f\n",(double)(integr.PhaseVect())[0]);
	}
	delete model;
	return 0;
}
