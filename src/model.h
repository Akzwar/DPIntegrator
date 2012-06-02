#include "../../VectTest/src/Vect.h"
#include <math.h>
#define mu  2.9591220829566038481357324012661e-4
#include <stdio.h>
class Model 
{
 protected:
	Vect IV;
 public:
	Model(){}
	virtual const Vect getRight( const Vect& Vt, long double t ){printf("Fuuuuuuuuuuuuuuuu\n");}
	Vect InitVect(){ return IV; } 
};

class KeplerModel : public Model
{
 public:
	KeplerModel() : Model()
	{
		IV = Vect( 6 );
		IV[0] = 0.81397168111;//X
		IV[1] = -0.010143699142;//dX
		IV[2] = 0.5225171672;//...
		IV[3] = 0.012874795077;
		IV[4] = 0.22658533987;
		IV[5] = 0.0055828817986;
	}	
	virtual const Vect getRight( const Vect& Vt, long double t )
	{
		
		Vect V = Vect(Vt);
		Vect retVect( V.size() );
		long double R;
		
		R = sqrt( V[0] * V[0] + V[2] * V[2] + V[4] * V[4] );
		retVect[0] = V[1];
		retVect[2] = V[3];
		retVect[4] = V[5];
		retVect[1] = -mu * V[0]/( R * R * R );
		retVect[3] = -mu * V[2]/( R * R * R );
		retVect[5] = -mu * V[4]/( R * R * R );
		return retVect;
	}
};
