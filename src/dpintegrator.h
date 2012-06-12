#include "model.h"
#include <math.h>
class DPIntegrator
{
 private:
	long double t, t0, tk, step, h, Tout;
	Vect k[8];
	
	long double a[8][7];
	long double b[8];
	long double sb[8];
	long double b1[8];
	long double c[8];

	long double eps, eps_max, u;
	long double* keps;
	
	Vect OutVect;
	Vect CurrVect;
	bool issmallstep;
	Vect prevVect, tmpVect, epsVect;
 public:
	Model* M;
		
	DPIntegrator(){}
	~DPIntegrator()
	{
		if( keps != NULL )
			delete [] keps;
	}
	DPIntegrator( Model* _M, long double _t0, long double _tk, long double _eps_max, long double _h = 0 )
	{
		eps = 0;
		eps_max = _eps_max;
		long double v = 1;
		while ( ( 1 + v ) > 1 )
		{
			u = v;
			v = v / 2;
		}
		t0 = _t0;
		M = _M;
		tk = _tk;
		step = 10e-5;
		CurrVect = M->InitVect();
		if( _h == 0 )
			issmallstep = false;
		else
		{
			issmallstep = true;
			h = _h;
		}
		a[2][1]=1.0/5.0;
            	a[3][1]=3.0/40.0;
            	a[3][2]=9.0/40.0;
            	a[4][1]=44.0/45.0;
            	a[4][2]=-56.0/15.0;
            	a[4][3]=32.0/9.0;
            	a[5][1]=19372.0/6561.0;
            	a[5][2]=-25360.0/2187.0;
            	a[5][3]=64448.0/6561.0;
            	a[5][4]=-212.0/729.0;
            	a[6][1]=9017.0/3168.0;
            	a[6][2]=-355.0/33.0;
            	a[6][3]=46732.0/5247.0;
            	a[6][4]=49.0/176.0;
            	a[6][5]=-5103.0/18656.0;
            	a[7][1]=35.0/384.0;
            	a[7][2]=0;
            	a[7][3]=500.0/1113.0;
            	a[7][4]=125.0/192.0;
            	a[7][5]=-2187.0/6784.0;
            	a[7][6]=11.0/84.0;
            	b[1]=35.0/384.0;
            	b[2]=0;
            	b[3]=500.0/1113.0;
            	b[4]=125.0/192.0;
            	b[5]=-2187.0/6784.0;
            	b[6]=11.0/84.0;
            	b[7]=0;
            	b1[1]=5179.0/57600.0;
            	b1[2]=0;
            	b1[3]=7571.0/16695.0;
            	b1[4]=393.0/640.0;
            	b1[5]=-92097.0/339200.0;
            	b1[6]=187.0/2100.0;
            	b1[7]=1.0/40.0;
            	c[1]=0;
            	c[2]=1.0/5.0;
            	c[3]=3.0/10.0;
            	c[4]=4.0/5.0;
            	c[5]=8.0/9.0;
            	c[6]=1.0;
            	c[7]=1.0;
		Tout = step+1;
		keps = new long double[CurrVect.size()];
		t = t0;
	}
		
	void SmallStep()
	{	
		if(Tout > step )
		{
			BigStep();
			Tout=0;
		}
		double teta = h/step;
		sb[1]=teta*(1+teta*(-1337.0/480.0 + teta*(1039.0/360.0 + teta*(-1163.0/1152.0))));
		sb[2]=0.0;
		sb[3] = 100.0 * teta * teta * (1054.0/9275.0 + teta* ( -4682.0/27825.0 + teta*( 379.0/5565.0 ) ) )/3.0; 
		sb[4] = -5.0*teta*teta*( 27.0/40.0 + teta * (-9.0/5.0 + teta *(83.0/96.0)) )/2.0;
		sb[5] = 18225.0*teta*teta*(-3.0/250.0+teta*(22.0/375.0+teta*(-37.0/600.0)))/848.0;
		sb[6] = -22.0*teta*teta*(-3.0/10.0+teta*(29.0/30.0+teta*(-17.0/24.0)))/7.0;
		OutVect = (( k[1] * sb[1] ) +
                            ( k[2] * sb[2] ) +
                            ( k[3] * sb[3] ) +
                            ( k[4] * sb[4] ) +
                            ( k[5] * sb[5] ) +
                            ( k[6] * sb[6] ) +
                            ( k[7] * sb[7] )) 
                                * h;
       		Tout+=h;
	}
	
	void BigStep()
	{
		Vect summ;
		long double s;
		eps = eps_max + 1;
		double u2 = 2.0 * u / eps_max;
		while ( eps > eps_max )
		{
			tmpVect = CurrVect;
			epsVect = CurrVect;
			prevVect = CurrVect;
			
			for( int i = 1; i <= 7; i++ )
			{
				summ = Vect(tmpVect.size());
				for( int j = 1; j < i; j++ )
				{
					summ += k[j] * a[i][j];
				}
				k[i] = M->getRight( tmpVect + ( summ * step ), step ); 
			}
			summ = Vect(tmpVect.size());
			for( int i = 1; i <= 7; i++ )
				summ += k[i] * b1[i];
			epsVect += summ * step;
			summ = Vect(tmpVect.size());
			for( int i = 1; i <= 7; i++ )
				summ += k[i] * b[i];
			tmpVect += summ*step;
			
			s = 0;
			for( int i = 0; i < CurrVect.size(); i++ )
			{
				keps[i] = step * ( tmpVect[i] - epsVect[i] )/
						max( max( 1.0e-5, sqrt( tmpVect[i] * tmpVect[i] ) ) , max( sqrt( prevVect[i] * prevVect[i]  ), u2 ) );
				s += keps[i] * keps[i]; 
			}
			eps = sqrt( s / 2.0 );

			step = step / max( 0.1, min( 5.0, pow( eps / eps_max, 0.2 ) / 0.9 ) );
		}	
	t += step;
	CurrVect = tmpVect;
	}

	void NextStep()
	{
		if(issmallstep)
		SmallStep();
		else
		BigStep();
	}
	
	Vect PhaseVect()
	{
		if(issmallstep)
		return OutVect;
		return CurrVect;
	}
	long double getT()
	{return this->t;}
	long double getStep()
	{return this->step;}
	long double getTk()
	{return this->tk;}
};
