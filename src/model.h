#include "../../Quat_Test/src/Quat.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#pragma once
template<class DD> class Model 
{
 protected:
	Vect<DD> IV;
    DD mu;
 public:
	Model(){}
    Model(DD p, DD e, DD i, DD Omega, DD omega, DD Theta)
    {
        mu = 3.986004418e5 * 3600 * 3600;
        IV = Vect<DD>( 6 );
        DD R = p / ( 1 + e * cos( Theta ));
        DD u = Theta + omega;

        IV[V_X] = R * ( cos( u ) * cos( Omega ) - sin( u ) * sin( Omega ) * cos( i ) );
        IV[V_Y] = R * ( cos( u ) * sin( Omega ) - sin( u ) * cos( Omega ) * cos( i ) );
        IV[V_Z] = R * sin( u ) * sin( i );

        DD Vr = sqrt( mu / p ) * e * sin( Theta ) / R;
        DD Vn = sqrt( mu / p ) * ( 1 + e * cos( Theta) );

        IV[V_dX] = IV[V_X] * Vr + Vn * ( -sin( u ) * cos( Omega ) 
                                        -cos( u ) * sin( Omega ) * cos( i ) );
        IV[V_dY] = IV[V_Y] * Vr + Vn * ( -sin( u ) * sin( Omega ) 
                                        +cos( u ) * cos( Omega ) * cos( i ) );
        IV[V_dZ] = IV[V_Z] * Vr + Vn *cos( u ) * sin( i );
    } 
	virtual const Vect<DD> getRight( const Vect<DD>& Vt, DD t ){}
	Vect<DD> InitVect(){ return IV; }
};

template<class DD> class KeplerModel : public Model<DD>
{
 public:
    KeplerModel() : Model<DD>(){}
	KeplerModel( double R ) : Model<DD>()
	{  
		this->mu = 2.9591220829566038481357324012661e-4;
		this->IV = Vect<DD>( 6 );
		this->IV[ V_X ] = 0.11601490913916648627 * R;
		this->IV[ V_Y ] = -0.92660555364038517604 * R;
		this->IV[ V_Z ] = -0.40180627760698804496 * R;
		this->IV[ V_dX ] = 0.01681162005220228976 * R;
		this->IV[ V_dY ] = 0.00174313168798203152 * R;
		this->IV[ V_dZ ] = 0.00075597376713614610 * R; 
	}	
    KeplerModel(DD p, DD e, DD i, DD Omega, DD omega, DD Theta) : Model<DD>( p,e,i,Omega,omega,Theta ){}

	virtual const Vect<DD> getRight( const Vect<DD>& Vt, DD t )
	{
		
		Vect<DD> V(Vt);
		Vect<DD> retVect( V.size() );
		DD R = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] + V[V_Z] * V[V_Z] );
		retVect[V_X] = V[V_dX];
		retVect[V_Y] = V[V_dY];
		retVect[V_Z] = V[V_dZ];
		retVect[V_dX] = -this->mu * V[V_X]/( R * R * R );
		retVect[V_dY] = -this->mu * V[V_Y]/( R * R * R );
		retVect[V_dZ] = -this->mu * V[V_Z]/( R * R * R );
		return retVect;
	}
};

struct MassDot
{
	double mass,X,Y,Z;
};

template<class DD>
class DotMassesModel : public Model<DD>
{
 private:
	MassDot MD[8];
	DD Emass;
	DD aem;
 public:
	DotMassesModel() : Model<DD>()
	{	
		Emass = 5.9742*10e23;
		MD[0].mass =  10652.73; MD[0].X = -4200.917; MD[0].Y = 3067.194;  MD[0].Z = -217.754;
		MD[1].mass =  -4160.24; MD[1].X = 856.705;   MD[1].Y = 5487.255;  MD[1].Z = -599.402;
		MD[2].mass =  41357.43; MD[2].X = 2843.858;  MD[2].Y = -1911.723; MD[2].Z = 3151.231;
		MD[3].mass = -48929.27; MD[3].X = 1469.622;  MD[3].Y = -3260.508; MD[3].Z = 2157.722;
		MD[4].mass =   3856.77; MD[4].X = 2825.441;  MD[4].Y = 2630.083;  MD[4].Z = -4148.659;
		MD[5].mass = -50642.30; MD[5].X = -1259.093; MD[5].Y = 271.749;   MD[5].Z = -3997.510;
		MD[6].mass =   6602.59; MD[6].X = -5413.541; MD[6].Y = 531.582;	  MD[6].Z = -1948.403;
		MD[7].mass = -19965.52; MD[7].X = 184.686;   MD[7].Y = 4435.039;  MD[7].Z = 3698.865;
        for( int i = 0; i < 8; i++)
            MD[i].mass *= 10e10;
	}
    DotMassesModel(DD p, DD e, DD i, DD Omega, DD omega, DD Theta) : Model<DD>( p,e,i,Omega,omega,Theta ) 
	{
		Emass = 5.9742*10e23;
		MD[0].mass =  10652.73; MD[0].X = -4200.917; MD[0].Y = 3067.194;  MD[0].Z = -217.754;
		MD[1].mass =  -4160.24; MD[1].X = 856.705;   MD[1].Y = 5487.255;  MD[1].Z = -599.402;
		MD[2].mass =  41357.43; MD[2].X = 2843.858;  MD[2].Y = -1911.723; MD[2].Z = 3151.231;
		MD[3].mass = -48929.27; MD[3].X = 1469.622;  MD[3].Y = -3260.508; MD[3].Z = 2157.722;
		MD[4].mass =   3856.77; MD[4].X = 2825.441;  MD[4].Y = 2630.083;  MD[4].Z = -4148.659;
		MD[5].mass = -50642.30; MD[5].X = -1259.093; MD[5].Y = 271.749;   MD[5].Z = -3997.510;
		MD[6].mass =   6602.59; MD[6].X = -5413.541; MD[6].Y = 531.582;	  MD[6].Z = -1948.403;
		MD[7].mass = -19965.52; MD[7].X = 184.686;   MD[7].Y = 4435.039;  MD[7].Z = 3698.865;
        for( int i = 0; i < 8; i++)
            MD[i].mass *= 10e10;
	}

	virtual const Vect<DD> getRight( const Vect<DD>& Vt, DD t )
	{
		Vect<DD> V( Vt );
		Vect<DD> retVect( V.size() );
		Vect<DD> summ( V.size() );
		DD dX, dY, dZ, R, eps, X, Y, Z;
		for( int i = 0; i < 8; i++ )
		{
			dX = MD[i].X - V[V_X];
			dY = MD[i].Y - V[V_Y];
			dZ = MD[i].Z - V[V_Z];
			R =  dX * dX + dY * dY + dZ * dZ ;
			eps = MD[i].mass / Emass;
			X = -eps * dX * 0.5 / pow(R,1/3.0);
			Y = -eps * dY * 0.5 / pow(R,1/3.0);
			Z = -eps * dZ * 0.5 / pow(R,1/3.0);
			summ += Vect<DD>( X, Y, Z ) ; 
		}
		summ = summ * this->mu;
		retVect[V_X] = V[V_dX];
		retVect[V_Y] = V[V_dY];
		retVect[V_Z] = V[V_dZ];
		retVect[V_dX] = summ[V_X];
		retVect[V_dY] = summ[V_Y];
		retVect[V_dZ] = summ[V_Z];
		return retVect;	
	}
};

	
double LegandrSin( int n, int m, double Phi )
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
		return LegandrSin( n-1, m-1, Phi ) * cos(Phi) * sqrt( (2*n+1)/(double)(2*n)*1/sm );
	return LegandrSin( n-1, m, Phi ) * sin(Phi) * sqrt( (4*n*n - 1)/(double)(n*n-m*m) ) - LegandrSin( n-2, m, Phi ) * sqrt((( (n-1)*(n-1) )-m*m)*(2*n+1)/(double)((n*n-m*m)*(2*n-3)));
}

double dLegandrSin( int n, int m, double Phi )
{
	double sm;
	if( m == 0 )
		sm = 0.5;
	else 
		sm = 1;
	return -( m * tan(Phi) * LegandrSin( n, m, Phi ) - sqrt( sm*(n-m)*(n+m+1) ) * LegandrSin( n, m+1, Phi ) );
}

template<typename DD>
class NormalSphereFuncModel : public Model<DD>
{
 private:
	static const double ae = 6378.136;
	static const long double C20 = -484165 * 10e-8;
	static const long double C40 = 790.3 * 10e-8;	
 public:
	NormalSphereFuncModel() : Model<DD>()
	{
	}
    NormalSphereFuncModel(DD p, DD e, DD i, DD Omega, DD omega, DD Theta) : Model<DD>( p,e,i,Omega,omega,Theta )
    {
    }

	virtual const Vect<DD> getRight( const Vect<DD>& Vt, DD t )
	{
		Vect<DD> V( Vt );
		Vect<DD> retVect( V.size() );
		Vect<DD> dV0( V.size() );
		DD Ro, Phi, Lambda, dPhi, dRo, dLambda;
		DD r0 = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] );
		Ro = V.Length();
		Phi = atan2( V[V_Z], r0 );
		Lambda = atan2( V[V_Y], V[V_X] );
        dRo = -this->mu/Ro/ae*( 3*pow(ae/Ro,3)*C20*LegandrSin(2,0,Phi)
                                + 5*pow(ae/Ro,5)*C40*LegandrSin(2,0,Phi)+1);
        dPhi = -this->mu/Ro/ae*( pow(ae/Ro,3)*C20*dLegandrSin(2,0,Phi)
                                + pow(ae/Ro,5)*C40*dLegandrSin(4,0,Phi));
		dLambda = 0;
		dV0 = Vect<DD>(dRo, dPhi, dLambda);
		Matrix<DD> MdV0;
		MdV0.addRow(dV0);
		MdV0.Transpose();
		Matrix<DD> MP(3,3);
		MP[0][0] = V[V_X] / Ro; MP[0][1] = -V[V_X] * V[V_Z] / ( Ro * r0 ); MP[0][2] = -V[V_Y] / r0;
		MP[1][0] = V[V_Y] / Ro; MP[1][1] = -V[V_Y] * V[V_Z] / ( Ro * r0 ); MP[1][2] = V[V_X] / r0;
		MP[2][0] = V[V_Z] / Ro; MP[2][1] = -r0 / Ro;                       MP[2][2] = 0;
		Matrix<DD> res = MP * MdV0;
		retVect[V_X] = V[ V_dX ];
		retVect[V_Y] = V[ V_dY ];
		retVect[V_Z] = V[ V_dZ ];
		retVect[V_dX] = res[0][0];
		retVect[V_dY] = res[1][0];
		retVect[V_dZ] = res[2][0];
		return retVect;
	}
};

template<class DD>
class AnomalSphereFuncModel : public Model<DD>
{
 private:
	double C[7][7],S[7][7];
	static const double ae = 6378.136;
 public:
	AnomalSphereFuncModel() : Model<DD>()
	{
        for(int i = 0; i < 7; i++)
		    for( int j =0; j<7;j++ )
			{
				C[i][j] = 0;
				S[i][j] = 0;
			}
		C[1][0] = 0; C[1][1] = 0; 
		C[2][0] = -484164.95;	C[2][1] = 0.05; C[2][2] = 2438.76;
		C[3][0] = 957.16; C[3][1] = 2032.83; C[3][2] = 906.86; C[3][3] = 715.99;
		C[4][0] = 543.52; C[4][1] = -531.34; C[4][2] = 350.04; C[4][3] = 989.28; C[4][4] = -193.25;
		C[5][0] = 68.72;  C[5][1] = -67.36;  C[5][2] = 651.94; C[5][3] = -451.45; C[5][4] = -293.19;
		C[5][5] = 191.57; 
		C[6][0] = -143.45; C[6][1] = -77.66; C[6][2] = 38.02; C[6][3] = 56.65; C[6][4] = -90.89;
		C[6][5] = -268.71; C[6][6] = 14.97;
		S[1][1] = 0.0; S[2][1] = 0.01; S[2][2] = -1399.68; S[3][1] = 249.03; S[3][2] = -619.29;
		S[3][3] = 1400.74; S[4][1] = -471.95; S[4][2] = 651.66; S[4][3] = -197.67; S[4][4] = 298.76;
		S[5][1] = -84.56; S[5][2] = -326.28; S[5][3] = -201.09; S[5][4] = 54.38; S[5][5] = -669.62;
		S[6][1] = 32.21; S[6][2] = -364.79; S[6][3] = 2.64; S[6][4] = -465.16; S[6][5] = -535.01;
		S[6][6] = -243.90;
	}
    AnomalSphereFuncModel(DD p, DD e, DD i, DD Omega, DD omega, DD Theta) : Model<DD>( p,e,i,Omega,omega,Theta )
    {
        for(int i = 0; i < 7; i++)
		    for( int j =0; j<7;j++ )

			{
				C[i][j] = 0;
				S[i][j] = 0;
			}
		C[1][0] = 0; C[1][1] = 0; 
		C[2][0] = -484164.95;	C[2][1] = 0.05; C[2][2] = 2438.76;
		C[3][0] = 957.16; C[3][1] = 2032.83; C[3][2] = 906.86; C[3][3] = 715.99;
		C[4][0] = 543.52; C[4][1] = -531.34; C[4][2] = 350.04; C[4][3] = 989.28; C[4][4] = -193.25;
		C[5][0] = 68.72;  C[5][1] = -67.36;  C[5][2] = 651.94; C[5][3] = -451.45; C[5][4] = -293.19;
		C[5][5] = 191.57; 
		C[6][0] = -143.45; C[6][1] = -77.66; C[6][2] = 38.02; C[6][3] = 56.65; C[6][4] = -90.89;
		C[6][5] = -268.71; C[6][6] = 14.97;
		S[1][1] = 0.0; S[2][1] = 0.01; S[2][2] = -1399.68; S[3][1] = 249.03; S[3][2] = -619.29;
		S[3][3] = 1400.74; S[4][1] = -471.95; S[4][2] = 651.66; S[4][3] = -197.67; S[4][4] = 298.76;
		S[5][1] = -84.56; S[5][2] = -326.28; S[5][3] = -201.09; S[5][4] = 54.38; S[5][5] = -669.62;
		S[6][1] = 32.21; S[6][2] = -364.79; S[6][3] = 2.64; S[6][4] = -465.16; S[6][5] = -535.01;
		S[6][6] = -243.90;

    }
	virtual const Vect<DD> getRight( const Vect<DD>& Vt, DD t )
	{
		int N = 7;
		DD summ1 = 0;
		DD summ2 = 0;
		Vect<DD> V(Vt);
		Vect<DD> retVect( V.size() );
		DD Ro, Phi, Lambda, dPhi, dRo, dLambda;
		DD r0 = sqrt( V[V_X] * V[V_X] + V[V_Y] * V[V_Y] );
		Ro = V.Length();
        printf("%f\n",Ro);
		Phi = atan2( V[V_Z], r0 );
		Lambda = atan2( V[V_Y], V[V_X] );
		for( int n = 2; n < N; n++ )
		{
			summ2 = 0;
			for( int m = 0; m <= n; m++ )
				summ2 += (C[n][m]*cos(m*Lambda) + S[n][m]*sin(m*Lambda))*LegandrSin(n,m,Phi);
			summ1 += (n+1)*pow(ae/Ro,n+1)*summ2;
		}
		dRo = -this->mu * summ1 / (ae*Ro);
        printf("%f\n",dRo);
		summ1 = 0;
		for( int n = 2; n < N; n++ )
		{
			summ2 = 0;
			for( int m = 0; m <= n; m++ )
				summ2 += (C[n][m]*cos(m*Lambda) + S[n][m]*sin(m*Lambda))*dLegandrSin(n,m,Phi);
			summ1 += pow(ae/Ro,n+1)*summ2;
		}
		dPhi = -this->mu * summ1 / (ae*Ro);
		summ1 = 0;
		for( int n = 2; n < N; n++ )
		{
			summ2 = 0;

			for( int m = 0; m <= n; m++ )
				summ2 += (-C[n][m]*sin(m*Lambda) + S[n][m]*cos(m*Lambda))*m*LegandrSin(n,m,Phi);
			summ1 += pow(ae/Ro,n+1)*summ2;
		}
		dLambda = -this->mu * summ1 / (ae*Ro*cos(Phi));			
		Matrix<DD> MP(3,3);
                MP[0][0] = V[V_X] / Ro; MP[0][1] = -V[V_X] * V[V_Z] / ( Ro * r0 ); MP[0][2] = -V[V_Y] / r0;
                MP[1][0] = V[V_Y] / Ro; MP[1][1] = -V[V_Y] * V[V_Z] / ( Ro * r0 ); MP[1][2] = V[V_X] / r0;
                MP[2][0] = V[V_Z] / Ro; MP[2][1] = -r0 / Ro;                       MP[2][2] = 0;
		Matrix<DD> gS;
        Vect<double> dG0(dRo, dPhi, dLambda);
        //dG0.print();
		gS.addRow(dG0);
		gS.Transpose();
		Matrix<DD> res = MP * gS;
		retVect[V_X] = V[V_dX];
		retVect[V_Y] = V[V_dY];
		retVect[V_Z] = V[V_dZ];
		retVect[V_dX] = V[V_X]*(-this->mu/V.Length())+res[0][0];
		retVect[V_dY] = V[V_Y]*(-this->mu/V.Length())+res[1][0]; 
		retVect[V_dZ] = V[V_Z]*(-this->mu/V.Length())+res[2][0];
		return retVect;
	}	
};

inline double randn( double mu, double sigma )
{
	static bool deviateAvailable = false;
	static float storedDeviate;
	double polar, rsquared, var1, var2;
	if(!deviateAvailable)
	{
		do
		{
			var1 = 2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2 = 2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared = var1*var1 + var2*var2;
		}
		while ( rsquared>=1.0 || rsquared == 0.0 );
		polar = sqrt( -2.0*log(rsquared)/rsquared );
		storedDeviate = var1*polar;
		deviateAvailable = true;

		return var2*polar*sigma+mu;
	}
	else
	{
		deviateAvailable = false;
		return storedDeviate*sigma+mu;
	}

}

template<typename DD>
class ShapingFilterExp : public Model<DD>
{
 private: 
	DD D, lambda;
 public:
	ShapingFilterExp(DD _D, DD _lambda) : Model<DD>()
	{
		this->IV = Vect<DD>(2);
		this->IV[0] = 0;
		this->IV[1] = 0;
		D = _D;
		lambda = _lambda;
		srand(time(NULL));
	}
	virtual const Vect<DD> getRight( const Vect<DD>& Vt, DD t )
	{
		Vect<DD> V(Vt);
		Vect<DD> retVect(V.size());
		//long double nu = sqrt(exp(1/rand())*2*1/t)*sin(rand()*2*M_PI);
		DD nu = randn(0, sqrt(D));
		retVect[0] = V[1];
		retVect[1] = nu * sqrt(2*D*lambda)-V[0]*lambda;
		return retVect;
	}
};
