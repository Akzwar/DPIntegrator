#include "../../Quat_Test/src/Quat.h"
#include <math.h>
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
 protected:
	long double mu;
 public:
	KeplerModel() : Model()
	{  
		mu = 2.9591220829566038481357324012661e-4;

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
		long double R = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] + V[V_Z] * V[V_Z] );
		retVect[V_X] = V[V_dX];
		retVect[V_Y] = V[V_dY];
		retVect[V_Z] = V[V_dZ];
		retVect[V_dX] = -mu * V[V_X]/( R * R * R );
		retVect[V_dY] = -mu * V[V_Y]/( R * R * R );
		retVect[V_dZ] = -mu * V[V_Z]/( R * R * R );
		return retVect;
	}
};

struct MassDot
{
	double mass,X,Y,Z;
};

class DotMassesModel : public KeplerModel
{
 private:
	MassDot MD[8];
	long double Emass;
	long double aem;
 public:
	DotMassesModel() : KeplerModel()
	{	
		mu = 398600.436 * 10e9;
		IV[ V_X ] = -5000;
		IV[ V_Y ] = 0;
		IV[ V_Z ] = 0;
		IV[ V_dX ] = 0;
		IV[ V_dY ] = 0;
		IV[ V_dZ ] = 0;
		Emass = 5.9742*10e24;
		aem = 149597870.66;
		MD[0].mass =  10652.73; MD[0].X = -4200.917; MD[0].Y = 3067.194;  MD[0].Z = -217.754;
		MD[1].mass =  -4160.24; MD[1].X = 856.705;   MD[1].Y = 5487.255;  MD[1].Z = -599.402;
		MD[2].mass =  41357.43; MD[2].X = 2843.858;  MD[2].Y = -1911.723; MD[2].Z = 3151.231;
		MD[3].mass = -48929.27; MD[3].X = 1469.622;  MD[3].Y = -3260.508; MD[3].Z = 2157.722;
		MD[4].mass =   3856.77; MD[4].X = 2825.441;  MD[4].Y = 2630.083;  MD[4].Z = -4148.659;
		MD[5].mass = -50642.30; MD[5].X = -1259.093; MD[5].Y = 271.749;   MD[5].Z = -3997.510;
		MD[6].mass =   6602.59; MD[6].X = -5413.541; MD[6].Y = 531.582;	  MD[6].Z = -1948.403;
		MD[7].mass = -19965.52; MD[7].X = 184.686;   MD[7].Y = 4435.039;  MD[7].Z = 3698.865;
		for( int i = 0; i < 8; i++ )
		{
			//MD[i].X /= aem;
			//MD[i].Y /= aem;
			//MD[i].Z /= aem;
		}
	}
	virtual const Vect getRight( const Vect& Vt, long double t )
	{
		Vect V( Vt );
		Vect retVect( V.size() );
		Vect summ( V.size() );
		long double dX, dY, dZ, R, eps, X, Y, Z;
		for( int i = 0; i < 8; i++ )
		{
			dX = MD[i].X - V[V_X];
			dY = MD[i].Y - V[V_Y];
			dZ = MD[i].Z - V[V_Z];
			//dX = sqrt(dX*dX);
			//dY = sqrt(dY*dY);
			//dZ = sqrt(dZ*dZ);
			R = sqrt( dX * dX + dY * dY + dZ * dZ );
			eps = MD[i].mass * 10e10 / Emass;
			X = -eps * dX / pow( R, 3.0/2.0 );
			Y = -eps * dY / pow( R, 3.0/2.0 );
			Z = -eps * dZ / pow( R, 3.0/2.0 );
			//printf("%f,%f,%f\n",(double)X*10e10, (double)Y, (double)Z);
			summ += Vect( X, Y, Z ) ; 
		}
		summ = summ * mu;
		//summ.print();
		retVect[V_X] = V[V_dX];
		retVect[V_Y] = V[V_dY];
		retVect[V_Z] = V[V_dZ];
		retVect[V_dX] = summ[V_X];
		retVect[V_dY] = summ[V_Y];
		retVect[V_dZ] = summ[V_Z];
		return retVect;	
	}
};
