#include <stdio.h>

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <unistd.h>
#include <math.h>
#include "dpintegrator.h"
#include <vector>

#define M_3D 0
#define M_2D 1

int sign( double Value )
{
    if( Value < 0 )
        return -1;
    if( Value > 0 )
        return 1;
    if( Value == 0 )
        return 0;
}

int W=1280;
int H=1024;

using namespace std;

int window;
Model<double>* model;
DPIntegrator<double>* Integr;
Quat<double> rotVect,Q,Qm1;
Vect<double> PV;
int mode;
int  l = 0;

struct Gnomon
{
	Vect<double> V;
	double l;
};

class Globe
{
 private:
	vector<float> InitSphere;
	int SphereDotsCount;
	vector<float> Sphere;
 public:
	Gnomon IG,G;
	Globe(){}
	Globe(double R, double _l, double _alpha, double _beta)
	{
		G.V = Vect<double>(R*sin(_beta)*sin(_alpha),R*cos(_beta),R*sin(_beta)*cos(_alpha));
		G.l = _l;
		IG  = G;
		SphereDotsCount = 0;
		for( double alpha = 0; alpha <= M_PI; alpha += M_PI/16.0 )
			for( double beta = 0-0.0001; beta <=2* M_PI; beta += M_PI/16.0 )
			{
				Sphere.push_back(R * sin(beta) * sin(alpha));
				Sphere.push_back(R * cos(beta));
				Sphere.push_back(R * sin(beta) * cos(alpha));
				Sphere.push_back(R * sin(beta + M_PI/16.0) * sin(alpha + M_PI/16.0));
				Sphere.push_back(R * cos(beta + M_PI/16.0));
				Sphere.push_back(R * sin(beta + M_PI/16.0) * cos(alpha + M_PI/16.0));	
				SphereDotsCount +=2;
			}
		InitSphere = Sphere;
	}
	void Offset(double X, double Y, double Z)
	{
		Sphere.clear();
		G.V[V_X] = IG.V[V_X] + X;
		G.V[V_Y] = IG.V[V_Y] + Y;
		G.V[V_Z] = IG.V[V_Z] + Z; 
		for( int i = 0; i < SphereDotsCount*3; i+=3 )
		{		
			Sphere.push_back( InitSphere[i] + X );
			Sphere.push_back( InitSphere[i+1] + Y );
			Sphere.push_back( InitSphere[i+2] + Z );
		}
	}
	void Rotate(double t)
	{
		Quat<double> res,q,qm1,tmp;	
		double angle;
		angle =  t * 24.0 * 60.0 * 60.0 * 7.292115*10e-6 / 2 ;
		q = Quat <double> (cos(angle),Vect<double>(0,sin(angle),0));
		qm1 = Quat <double> (cos(angle), Vect<double>(0,-sin(angle),0));
		tmp = Quat <double> (0,Vect<double>(IG.V));
		res = q * tmp * qm1; 
		IG.V = res.V;
		for( int i = 0; i < SphereDotsCount*3; i+=3 )
		{
			tmp = Quat<double>(0,Vect<double>(InitSphere[i],InitSphere[i+1],InitSphere[i+2]));
			res = q * tmp * qm1;
			InitSphere[i] = res.V[V_X];
			InitSphere[i+1] = res.V[V_Y];
			InitSphere[i+2] =res.V[V_Z];
		}
	}
	void* ptr(){return &Sphere[0];}
	int Count(){return SphereDotsCount;}
};

Globe Earth;

void ReinitCamera()
{	
	gluLookAt(rotVect.V[V_X], rotVect.V[V_Y], rotVect.V[V_Z],
 		0,0,0, 
		0,1,0);
}

void Draw3D()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(70.0,W/(double)H,0.01,10000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    ReinitCamera();	

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glColor3f(0.0,0.0,1.0);
	glDrawArrays(GL_POINTS,0,Earth.Count());
	
    glColor3f(0.0,1.0,1.0);
    glBegin(GL_LINE_STRIP);
        for( int i = 0; i < l; i += 10 )
            glVertex3f( Integr->Result[i].Y * 10e-6, 
                        Integr->Result[i].Z * 10e-6,
                        Integr->Result[i].X * 10e-6 );
    glEnd();
}

void Draw2D()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0,W,H,0,0,1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); 
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

void DrawGLScene()
{
    if(mode == M_2D)
        Draw2D();
    if(mode == M_3D)
        Draw3D();
 
    glutSwapBuffers(); 	
}

bool IsDrag_Left;

int prevX, prevY;
void onMouseClick(int button, int state, int x, int y)
{
	if( button == GLUT_LEFT_BUTTON )
	{
		if( state == GLUT_DOWN )
		{
			IsDrag_Left = true;
			prevX = x;
			prevY = y;
		}
		if( state == GLUT_UP )
			IsDrag_Left = false;
	}
	/*if( button == GLUT_MIDDLE_BUTTON )
	{
		if( state == GLUT_DOWN )
		{
			IsDrag_Middle = true;
			prevX = x;
			prevY = y;
		}
		if( state == GLUT_UP )
			IsDrag_Middle = false;
	}*/
}

Vect<double> rotAxis;
void onMouseMove( int X, int Y )
{
    if( IsDrag_Left )
    {
    	int kx = 500;
    	int ky = 500;
    	int dX = X - prevX;
    	int dY = Y - prevY;
    		double Angle = -(double)dX / kx;
    		Q = Quat<double>(cos(Angle / 2.0 ),Vect<double>(0,sin(Angle / 2.0 ),0));
    		Qm1 = Quat<double>(cos(Angle / 2.0),Vect<double>(0,-sin(Angle / 2.0 ),0));
    		rotVect = Q * rotVect * Qm1;
    		Angle = -(double)dY / ky;
            rotAxis = Vect<double>(0,1,0) % rotVect.V;
    		Q = Quat<double>(cos(Angle / 2.0 ), rotAxis * (1 / rotAxis.Length()) * sin(Angle / 2.0) );
    		Qm1 = Quat<double>(cos(Angle / 2.0), rotAxis * (1 / rotAxis.Length()) * ( - sin(Angle / 2.0 ) ) );
    		rotVect = Q * rotVect * Qm1;
            if( fabs(rotVect.V[V_X]) <= 0.01 && fabs(rotVect.V[V_Z]) <= 0.01 )
            {
                rotVect.V[V_X] = sign(rotVect.V[V_X]) * 0.01;
                rotVect.V[V_Z] = sign(rotVect.V[V_Z]) * 0.01;
            }
    	prevX = X;
    	prevY = Y;
}

}

void keyPressed(unsigned char key, int x, int y)
{
	if( key == 's' )
		rotVect.V = rotVect.V * 1.1;
	if( key == 'w' )
		rotVect.V = rotVect.V * 0.9;	
    if( key == 'g' )
        if(mode == M_3D)
            mode = M_2D;
        else
            mode = M_3D;
}

void ReSizeGLScene(int Width, int Height)
{
	glViewport(0,0,Width,Height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70.0, Width/(double)Height, 0.01, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	ReinitCamera();
    W = Width;
    H = Height;
}

void InitGL(int Width, int Height)
{
	glClearColor( 0,0,0.07,0 );
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	ReSizeGLScene(Width, Height);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3,GL_FLOAT,0,Earth.ptr());
	glPointSize(3.0);

}

void onIdle(void)
{
    if( l < Integr->Result.size() )
    {
        l+=30;
    }
	glutPostRedisplay();
} 

void IntegrateModel()
{
    while(Integr->getT()<=Integr->getTk())
	    {
		    Integr->NextStep();
	    }
}

int main(int argc, char **argv)
{
	model = new KeplerModel<double>(10000,0.2,M_PI/4.0,M_PI/2.0,40/180.0*M_PI,0);
	Integr = new DPIntegrator<double>(model,0,3,10e-18);
	rotVect = Quat<double>(0,Vect<double>(0,0,0.3));
	PV = Vect<double>(3);	
	printf("GL inititialization start...\n");
	Earth = Globe(6378.137*10e-6,0,0,0);
    IntegrateModel();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
	glutInitWindowSize(W, H);
	glutInitWindowPosition(100, 100);

	window = glutCreateWindow("OpenGLTutor");
	glutDisplayFunc(&DrawGLScene);
	glutIdleFunc(&onIdle);
	glutReshapeFunc(&ReSizeGLScene);
	glutKeyboardFunc(&keyPressed);
	glutMouseFunc(onMouseClick);
	glutMotionFunc(onMouseMove);
	InitGL(W,H);
	printf("GL Init done...\n");
	glutMainLoop();
	delete Integr;
	delete model;
	return 0;
}
