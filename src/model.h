#include "../../Quat_Test/src/Quat.h"
#include <math.h>
#define mu  2.9591220829566038481357324012661e-4
#include <stdio.h>
class Model 
{
 protected:
	Vect IV;
 public:
	Model(){}
	virtual const Vect getRight( const Vect& Vt, long double t ){}
	Vect InitVect(){ return IV; }
};

class KeplerModel : public Model
{
 public:
	KeplerModel() : Model()
	{
		IV = Vect( 6 );
		IV[ V_X ] = 0.11601490913916648627;
		IV[ V_Y ] = -0.92660555364038517604;
		IV[ V_Z ] = -0.40180627760698804496;
		IV[ V_dX ] = 0.01681162005220228976;
		IV[ V_dY ] = 0.00174313168798203152;
		IV[ V_dZ ] = 0.00075597376713614610; 
	}	
	virtual const Vect getRight( const Vect& Vt, long double t )
	{
		
		Vect V = Vect(Vt);
		Vect retVect( V.size() );
		long double R;
		
		R = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] + V[V_Z] * V[V_Z] );
		retVect[V_X] = V[V_dX];
		retVect[V_Y] = V[V_dY];
		retVect[V_Z] = V[V_dZ];
		retVect[V_dX] = -mu * V[V_X]/( R * R * R );
		retVect[V_dY] = -mu * V[V_Y]/( R * R * R );
		retVect[V_dZ] = -mu * V[V_Z]/( R * R * R );
		return retVect;
	}
};
