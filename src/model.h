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

	
long double LegandrSin( int n, int m, long double Phi )
{
	double sm;
	if( (m - 1) == 0 )
		sm = 0.5;
	else 
		sm = 1;
	if( n < m )
		return 0;
	if( n == m && m == 0 )
		return 1;
	if( n == m && m!=0 )
		return LegandrSin( n-1, m-1, Phi ) * cos(Phi) * sqrt( (2*n+1)/(float)(2*n*sm) );
	return LegandrSin( n-1, m, Phi ) * sin(Phi) * sqrt( (4*n*n - 1)/(n*n-m*m) ) - LegandrSin( n-2, m, Phi ) * sqrt((( (n-1)*(n-1) )-m*m)*(2*n+1)/((n*n-m*m)*(2*n-3)));
}

long double dLegandrSin( int n, int m, long double Phi )
{
	double sm;
	if( m == 0 )
		sm = 0.5;
	else 
		sm = 1;
	return -( m * tan(Phi) * LegandrSin( n, m, Phi ) - sqrt( sm*(n-m)*(n+m+1) ) * LegandrSin( n, m+1, Phi ) );
}

class NormalSphereFuncModel : public Model
{
 private:
	static const double ae = 6378.136;
	static const long double C20 = -484165 * 10e-9;
	static const long double C40 = 790.3 * 10e-9;	
	long double mu;
	long double P20,P40;
 public:
	NormalSphereFuncModel() : Model()
	{
		mu = 398600.436*10e9;
		IV = Vect( 6 );
		IV[ V_X ] = 0.11601490913916648627*11000;
		IV[ V_Y ] = -0.92660555364038517604*11000;
		IV[ V_Z ] = -0.40180627760698804496*11000;
		IV[ V_dX ] = 0.11681162005220228976*11000;
		IV[ V_dY ] = 0.00174313168798203152*11000;
		IV[ V_dZ ] = 0.00075597376713614610*11000; 

	}
	virtual const Vect getRight( const Vect& Vt, long double t )
	{
		Vect V( Vt );
		Vect retVect( V.size() );
		Vect dV0( V.size() );
		long double Ro, Phi, Lambda, dPhi, dRo, dLambda, R;
		long double r0 = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] );
		R = V.Length();
		Ro = R;
		Phi = atan( V[V_Z]/r0 );
		Lambda = atan( V[V_Y] / V[V_X] );
		//printf("SC=");
		//Vect(Ro,Phi,Lambda).print();
		//dRo = -mu / ( Ro * Ro ) - 3 * C20 * mu * ae * ae * LegandrSin( 2, 0 , Phi ) / pow( Ro, 4 ) - 5 * C40 * mu * pow(ae,4) * LegandrSin( 4,0,Phi ) / pow( Ro, 6 );
		//dRo = -(3*C20*ae*ae*LegandrSin(2,0,Phi)*Ro*Ro + 5*C40*pow(ae,4)*LegandrSin(4,0,Phi)+pow(Ro,4))/pow(Ro,6);
		//dPhi = mu/Ro * (C20 * pow(ae/Ro,2)*dLegandrSin(2,0,Phi)+C40*pow(ae/Ro,4)*dLegandrSin(4,0,Phi));
		dRo = -mu/Ro/ae*(3*pow(ae/Ro,3)*C20*LegandrSin(2,0,Phi)+5*pow(ae/Ro,5)*C40*LegandrSin(4,0,Phi));
		dPhi = -mu/Ro/ae*(pow(ae/Ro,3)*C20*dLegandrSin(2,0,Phi)+pow(ae/Ro,5)*C40*dLegandrSin(4,0,Phi));
		printf("-----%f\n",(double)dRo);
		dLambda = 0;
		dV0 = Vect(dRo, dPhi, dLambda);
		//dV0.print();
		Matrix MdV0;
		MdV0.addRow(dV0);
		MdV0.Transpose();
		Matrix MP(3,3);
		MP[0][0] = V[V_X] / Ro; MP[0][1] = -V[V_X] * V[V_Z] / ( Ro * r0 ); MP[0][2] = -V[V_Y] / r0;
		MP[1][0] = V[V_Y] / Ro; MP[1][1] = -V[V_Y] * V[V_Z] / ( Ro * r0 ); MP[1][2] = V[V_X] / r0;
		MP[2][0] = V[V_Z] / Ro; MP[2][1] = -r0 / Ro;			   MP[2][2] = 0;
		Matrix res = MP * MdV0;
		retVect[V_X] = V[ V_dX ];
		retVect[V_Y] = V[ V_dY ];
		retVect[V_Z] = V[ V_dZ ];
		retVect[V_dX] = res[0][0];
		retVect[V_dY] = res[1][0];
		retVect[V_dZ] = res[2][0];
		return retVect;
	}
};

class AnomalSphereFuncModel : public Model
{
 private:
	double C[37][37],S[37][37];
	long double mu;
 public:
	AnomalSphereFuncModel() : Model()
	{
		C[1][0] = 0; C[1][1] = 0; 
		C[2][0] = -484164.95;	C[2][1] = 0.05; C[2][2] = 2438.76;
		C[3][0] = 957.16; C[3][1] = ; C[3][2] = ; C[3][3] = ;
                mu = 398600.436*10e9;
                IV = Vect( 6 );
                IV[ V_X ] = 0.11601490913916648627*11000;
                IV[ V_Y ] = -0.92660555364038517604*11000;
                IV[ V_Z ] = -0.40180627760698804496*11000;
                IV[ V_dX ] = 0.11681162005220228976*11000;
                IV[ V_dY ] = 0.00174313168798203152*11000;
                IV[ V_dZ ] = 0.00075597376713614610*11000;
	}
	virtual const getRight( const Vect& V, long double t )
	{
		int N = 35;
		long double summ1 = 0;
		long double summ2 = 0;
		Vect V(Vt);
		Vect retVect( V.size() );
		long double Ro, Phi, Lambda, dPhi, dRo, dLambda;
		long double r0 = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] );
		Ro = V.Length();
		Phi = atan( V[V_Z] / r0 );
		Lambda = atan( V[V_Y] / V[V_X] );
		for( int n = 0; n < N; n++ )
		{
			summ2 = 0;
			for( int m = 0; m <= n; m++ )
				summ2 += (C[n][m]*cos(m*Lambda) + S[n][m]*sin(m*Lambda))*LegandrSin(n,m,Phi);
			summ1 += (n+1)*pow(ae/Ro,n+1)*summ2;
		}
		dRo = -mu * summ1 / (ae*Ro);
		summ1 = 0;
		for( int n = 0; n < N; n++ )
		{
			summ2 = 0;
			for( int m = 0; m <= n; m++ )
				summ2 += (C[n][m]*cos(m*Lambda) + S[n][m]*sin(m*Lambda))*dLegandrSin(n,m,Phi);
			summ1 += pow(ae/Ro,n+1)*summ2;
		}
		dPhi = -mu * summ1 / (ae*Ro);
		summ1 = 0;
		for( int n = 0; n < N; n++ )
		{
			summ2 = 0;
			for( int m = 0; m <= n; m++ )
				summ2 += (-C[n][m]*sin(m*Lambda) + S[n][m]*cos(m*Lambda))*m*LegandrSin(n,m,Phi);
			summ1 += pow(ae/Ro,n+1)*summ2;
		}
		dLambda = -mu * summ1 / (ae*Ro*cos(Phi));			
		Matrix MP(3,3);
                MP[0][0] = V[V_X] / Ro; MP[0][1] = -V[V_X] * V[V_Z] / ( Ro * r0 ); MP[0][2] = -V[V_Y] / r0;
                MP[1][0] = V[V_Y] / Ro; MP[1][1] = -V[V_Y] * V[V_Z] / ( Ro * r0 ); MP[1][2] = V[V_X] / r0;
                MP[2][0] = V[V_Z] / Ro; MP[2][1] = -r0 / Ro;                       MP[2][2] = 0;
		Matrix gS;
		gs.AddRow(Vect(dRo,dPhi,dLambda));
		gs.Transpose();
		Matrix res = MP * gS;
		retVect = V * (-mu/V.Length()) + Vect(res[0][0],res[1][0],res[2][0]);
		return retVect;
	}	
};
