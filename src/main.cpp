#include <stdio.h>
#include "dpintegrator.h"

int main(int argc, char** argv)
{
	ShapingFilterExp* model = new ShapingFilterExp(5,1);
	DPIntegrator integr( model, 0, 100, 10e-8 );
	while( integr.getT() < integr.getTk() )
	//for(int i=0;i<200; i++)
	{
		Vect Out = integr.PhaseVect();
		Out.print();
		integr.NextStep();	
	}
	//printf("%f",(double)LegandrSin(4,0,1));
	return 0;
}
