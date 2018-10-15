// p2dp3d.h
#if !defined( _P2D_H_ )
#define _P2D_H_ 
class p2d {
   
public:

	p2d(double x1=100000.0, double y1=100000.0);

	void set_x (double x1) {x=x1;}

	double get_x () {return x;}

	void set_y (double y1) {y=y1;}

	double get_y () {return y;}

	virtual ~p2d(){
	//	cout << "Edge Object is deleted! "<<endl;
		//delete[] p2;
	}

protected:
	double x;
	double y;
	p2d* p2;

};

p2d::p2d(double x1, double y1) 
{
x=x1;
y=y1;
}

#endif // _P2D_H_ 

#if !defined( _P3D_H_ )
#define _P3D_H_ 

class p3d {

public:	

	p3d(double x1=100000.0, double y1=100000.0, double z1=100000.0);
	
	void set_x (double x1) {x=x1;}

	double get_x () {return x;}

	void set_y (double y1) {y=y1;}

	double get_y () {return y;}

    void set_z (double z1) {z=z1;}

	double get_z () {return z;}

	virtual ~p3d(){
//		delete[] p3;
	//	cout << "3d Object is deleted! "<<endl;
		}
protected:
	double x;
	double y;
	double z;
	p3d* p3;

};

p3d::p3d(double x1, double y1,double z1) 
{
x=x1;
y=y1;
z=z1;
}
#endif // _P3D_H_ 
