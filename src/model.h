#include "../VectTest/src/Vect.h"

class Model 
{
 protected:
	Vect IV;
 public:
	Model(){}
	virtual Vect getRight( const Vect& V, long double t );
	Vect InitVect(){ return InitVect; } 
};
