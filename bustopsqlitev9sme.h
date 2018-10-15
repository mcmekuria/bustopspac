#include <iostream>
#include <assert.h>
#include <cstdlib>
//#include <sqlite3.h>
#include <deque>
#include <map>
#include <cstring>
#include <cstddef>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <time.h>

extern "C"
{
#include <sys/types.h>
#include <sqlite3.h>
#include <spatialite.h>

}


using namespace std;

//using size_t;


extern const int inf=100000;  // vertices not yet visited are initialized with -1
//#define eidconst 10000000 
extern const int eidconst=10000000;  // reverse edges are  inserted using eid = original eid + eidconst

#define P(A) cout << #A << " : " << (A) << endl;

// Vector Geometry 
// dot product (3D) which allows vector operations in arguments
#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define norm(v)    sqrt(dot(v,v))  // norm = length of vector
#define d(u,v)     norm(u-v)       // distance = norm of difference


class vertexp;

class Point;

template <class T >
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
};

template<typename U>  
U from_String(const string& s) {
  std::istringstream is(s);
  U u;
  is >> u;
  return u;
}

class Ptxy
{
public:
    double x, y;

    Ptxy(double d = 0.0)
        : x(d), y(d) {}

    Ptxy(double x1, double y1)
        : x(x1), y(y1) {}

	void set_x (double x1) {x=x1;}

	double get_x () {return x;}

	void set_y (double y1) {y=y1;}

	double get_y () {return y;}


};

class LineSegment
{
public:
    Ptxy begin;
    Ptxy end;

    LineSegment(const Ptxy& begin1, const Ptxy& end1)
        : begin(begin1), end(end1) {}

    enum IntersectResult { PARALLEL, COINCIDENT, NOT_INTERESECTING, INTERESECTING };

    IntersectResult Intersect(const LineSegment& other_line, Ptxy& intersection)
    {
        double denom = ((other_line.end.y - other_line.begin.y)*(end.x - begin.x)) -
                      ((other_line.end.x - other_line.begin.x)*(end.y - begin.y));

        double nume_a = ((other_line.end.x - other_line.begin.x)*(begin.y - other_line.begin.y)) -
                       ((other_line.end.y - other_line.begin.y)*(begin.x - other_line.begin.x));

        double nume_b = ((end.x - begin.x)*(begin.y - other_line.begin.y)) -
                       ((end.y - begin.y)*(begin.x - other_line.begin.x));

        if(denom == 0.0)
        {
            if(nume_a == 0.0 && nume_b == 0.0)
            {
                return COINCIDENT;
            }
            return PARALLEL;
        }

        double ua = nume_a / denom;
        double ub = nume_b / denom;

        if(ua >= 0.0 && ua <= 1.0 && ub >= 0.0 && ub <= 1.0)
        {
            // Get the intersection point.
            intersection.x = begin.x + ua*(end.x - begin.x);
            intersection.y = begin.y + ua*(end.y - begin.y);

            return INTERESECTING;
        }

        return NOT_INTERESECTING;
    }
};

Ptxy& DoLineSegmentIntersection(const Ptxy& p0, const Ptxy& p1, const Ptxy& p2, const Ptxy& p3)
{
    LineSegment linesegment0(p0, p1);
    LineSegment linesegment1(p2, p3);

    Ptxy intersection;

    std::cout << "Line Segment 0: (" << p0.x << ", " << p0.y << ") to (" << p1.x << ", " << p1.y << ")\n"
              << "Line Segment 1: (" << p2.x << ", " << p2.y << ") to (" << p3.x << ", " << p3.y << ")\n";

    switch(linesegment0.Intersect(linesegment1, intersection))
    {
    case LineSegment::PARALLEL:
        std::cout << "The lines are parallel\n\n";
		break;
    case LineSegment::COINCIDENT:
        std::cout << "The lines are coincident\n\n";
        break;
    case LineSegment::NOT_INTERESECTING:
        std::cout << "The lines do not intersect\n\n";
        break;
    case LineSegment::INTERESECTING:
        std::cout << "The lines intersect at (" << intersection.x << ", " << intersection.y << ")\n\n";
        break;
    }
return intersection;
}


#ifndef SS_Common_H
#define SS_Common_H

enum boolean {FALSE=0, TRUE=1, ERROR=(-1)};

// Error codes
enum Error {
	Enot,	// no error
	Edim,	// error: dim of point invalid for operation
	Esum	// error: sum not affine (cooefs add to 1)
};

// utility macros
//#define	abs(x)		((x) >= 0 ? x : -(x));
//#define	min(x,y)	((x) < (y) ? (x) : (y));
//#define	max(x,y)	((x) > (y) ? (x) : (y));

#endif SS_Common_H


class edgev {
public:
	edgev(const edgev& eg )  {
	id = eg.id; // edge id
	esid = eg.esid; // edge object id
	orig = eg.orig; // origin id
	evp = eg.evp;  // edge pointer 
	lbl = eg.lbl;  // edge label
	toStop = eg.toStop;  // (on or off cost)
	ecost= eg.ecost;  // edge cost
	scost= eg.scost;  // edge vertex start cost
	tcost= eg.tcost;  // edge vertex end cost
	frid = eg.frid;  // edge from vertex
	toid = eg.toid;  // edge to vertex
	palong = eg.palong; // position along
	tway = eg.tway;  // edge to vertex
	dirn = eg.dirn;  // edge by default set to boundary
	slen = eg.slen; // edge length
	lts = eg.lts; // Level of traffic Stress
	stopon = eg.stopon; // Boarding stop
	stopoff = eg.stopoff;  // Alighting Stop
	efc = eg.efc; // edge feature class (from database,ArcGIS)
	orfc = eg.orfc; // origin feature class (Spatialite Database or ArcGIS)
	ornm = eg.ornm; // origin feature class (Spatialite Database or ArcGIS)
	enote = eg.enote; // edge note such as fixed , free 
	//  cout << "Vx[" << pid << "]" << endl;
	//    ++copycons;
  }
  edgev& operator=(const edgev& eg) {
    //cout << "(" << pid << ")=[" << vx.pid << "]" << endl;
	id = eg.id; // edge id
	esid = eg.esid; // edge object id
	orig = eg.orig; // origin id
	evp = eg.evp;  // edge pointer 
	lbl = eg.lbl;  // edge label
	toStop = eg.toStop;  // (on or off cost)
	ecost= eg.ecost;  // edge cost
	scost= eg.scost;  // edge vertex start cost
	tcost= eg.tcost;  // edge vertex end cost
	frid = eg.frid;  // edge from vertex
	toid = eg.toid;  // edge to vertex
	palong = eg.palong; // position along
	tway = eg.tway;  // edge to vertex
	dirn = eg.dirn;  // edge by default set to boundary
	slen = eg.slen; //length of edge
	lts = eg.lts;  // Level of Traffic Stress value
	stopon = eg.stopon; // Boarding stop
	stopoff = eg.stopoff;  // Alighting Stop
	efc = eg.efc; // edge feature class (from database,ArcGIS)
	orfc = eg.orfc; // origin feature class (Spatialite Database or ArcGIS)
	ornm = eg.ornm; // origin feature class (Spatialite Database or ArcGIS)
	enote = eg.enote; // edge note such as fixed , free 
	    return *this;
  }

	bool operator==(edgev eg);

  	long get_cid () {return id;}

    friend istream& operator>>( istream&, edgev& );
	
	enum EnumLbl { fixed=-1, free=0, labeled = 1 };

	edgev(long id1, long esid1,long evp1, short lbl1, double ecost1, double scost1, double tcost1,double slen1, 
		long frid1, long toid1, long orig1, double palong1, short tway1,short dirn1,short toStop1, int lts1, 
		int stopoff1,int stopon1,string efc1,string orfc1,string ormn1,string enote1 ) ; // constructor
 
	edgev() : id(0), esid(0), evp(0),  lbl(0),  ecost(inf),  scost(inf),tcost(inf),
		  slen(0),palong(0),frid(0), toid(0), orig(-1), tway(-1), dirn(0),toStop(0),
		  stopoff(-1),stopon(-1),lts(999),efc(""),orfc(""),ornm(""),enote("")
		  { 
		//cout << "ev[" << id << "]" << " Create " << " create " << endl; 
	}; 
    void ev( long eid1=0,long esid1=0,long evp1=0, short lbl1=0, double ecost1=inf, double scost1=inf,double tcost1=inf,
		double slen1=0, long frid1=0, long toid1=0, long orig1=-1, double palong1=0.0, short tway1=-1,short dirn1=0, 
		short toStop1=0, int lts1=999,int stopon1=-1,int stopoff1=-1,string efc1="",string orfc1="",string ornm1="",
		string enote1="")
	{
	id = eid1; // edge id
	esid = esid1; // edge object id
	orig = orig1; // origin id
	evp = evp1;  // edge pointer 
	lbl = lbl1;  // edge label
	toStop = toStop1;  // (on or off cost)
	ecost= ecost1;  // edge cost
	scost= scost1;  // edge vertex start cost
	tcost= tcost1;  // edge vertex end cost
	frid = frid1;  // edge from vertex
	toid = toid1;  // edge to vertex
	palong = palong1; // position along
	tway = tway1;  // edge to vertex
	dirn = dirn1;  // edge by default set to boundary
	slen = slen1;  // edge by length
	lts = lts1; // level of traffic stress
	stopon=stopon1;  // stop assignment for boardings
	stopoff=stopoff1; // stop assignment for alightings
	efc = efc1; // edge feature class (from database,ArcGIS)
	orfc = orfc1; // origin feature class (Spatialite Database or ArcGIS)
	ornm = ornm1; // origin feature class (Spatialite Database or ArcGIS)
	enote = enote1; // edge note such as fixed , free 
	}

	void set_id (long id1) {id=id1;} // edge id that matches the street object id

	long get_id () const {return id;} // edge id that matches the street object id

	void set_orig (long orig1) {orig=orig1;}

	long get_orig () const {return orig;}

	void set_esid (long esid1) {esid=esid1;}

	long get_esid () const {return esid;}

	void set_evp (long evp1) {evp=evp1;}

	long get_evp () const {return evp;}

	void set_lbl (short lbl1) {lbl=lbl1;}

	short get_lbl ()  const {return lbl;}

	void set_dirn (short dirn1) {dirn=dirn1;}

	short get_dirn ()  const {return dirn;}

	void set_cost (double ecost1) {ecost=ecost1;}

	double get_cost ()  const {return ecost;}

	void set_frid (long frid1) {frid=frid1;}

	long get_frid ()  const {return frid;}

	void set_toid (long toid1) {toid=toid1;}

	long get_toid ()  const {return toid;}

	void set_tway (short tway1) {tway=tway1;}

	long get_tway ()  const {return tway;}

	void set_tcost (double tcost1) {tcost=tcost1;}

	double get_tcost ()  const {return tcost;}

	void set_palong (double palong1) {palong=palong1;}

	double get_palong ()  const {return palong;}

	void set_scost (double scost1) {scost=scost1;}

	double get_scost ()  const {return scost;}

	void set_slen (double slen1) {slen=slen1;}

	double get_slen ()  const {return slen;}

	void set_efc (string efc1) {efc=efc1;}

	string get_efc ()  const {return efc;}

	void set_ornm (string ornm1) {ornm=ornm1;}

	string get_ornm ()  const {return ornm;}

	void set_stopon (int stopon1) {stopon=stopon1;}

	int get_stopon ()  const {return stopon;}

	void set_stopoff (int stopoff1) {stopoff=stopoff1;}

	int get_stopoff ()  const {return stopoff;}

	void set_lts (int lts1) {lts=lts1;}

	int get_lts ()  const {return lts;}

	void set_enote (string enote1) {enote=enote1;}

	string get_enote ()  const {return enote;}

	void set_toStop (short toStop1) {toStop=toStop1;}

	short get_toStop ()  const {return toStop;}

	void show_edgehdr(ostream& out)
    { 
		out << "eid" << "\t" << "frid" << "\t" << "toid" << "\t" 
		<< "orig" << "\t" << "ecost"<< "\t" << "scost"<< "\t" << 
		"tcost" << "\t" <<"Length" << "\t" << "evp" << "\t" << "lbl" << endl;
	};
 	void show_edge(ostream& out)
    { 
		out << get_id() << "\t" << get_frid() << "\t" << get_toid() << "\t" 
		<< get_orig()<< "\t" << get_lts() << "\t" << get_cost() << "\t" << get_scost() << "\t" 
		<<get_tcost() << "\t" <<get_slen() << "\t" << get_evp() << "\t" << get_lbl() << endl;
	};
   void show_edge_Dialog(void)
    { 
		cout << "E-Id " << get_id() << " V(i) " << get_frid() << "V(j) " << get_toid() << "Edge-Orig " << get_orig()<< " lts " << get_lts()  
			<< " Cost " << get_cost()<< " sCost " << get_scost()<< " tCost " << get_tcost() << " E-pred " << get_evp() << " lbl " << get_lbl() << endl;
	};

	void read_edgehdr(istream& in)
    { 
		in >>id>> evp>> efc>>frid>>stopoff>>toid>>stopon  
		>> ecost >> palong>> orfc>>orig >>ornm>>scost>>tcost>>lbl>>enote>>esid 
		>>slen;
	};

	void read_edge(istream& in)
    { 
		in >>id>> evp>> efc>>frid>>stopoff>>toid>>stopon  
		>> ecost >> palong>> orfc>>orig >>ornm>>scost>>tcost>>lbl>>enote>>esid 
		>>slen;
	};

		void edgev::serialize(ofstream& pEdge);
		void edgev::deserialize(ifstream& pEdge);
		void edgev::serializetext(ofstream& pEdge);
		void edgev::serializetexthdr(ofstream& pEdge);

protected:
//	static long create, assign, copycons, destroy, sid;
			 long id;  // edge id that matches the street object id
	         long esid; // edge serial number
			 long evp;  // edge pointer
	         short lbl; // label
	         short dirn; // forward = 1, reverse=-1, boundary=-2,unscanned=0;
			 short	toStop; // 1 for cost of travel to a stop (ons) 0 otherwise (offs) 
			 double ecost; // edge cost
         	 double scost; // start path cost
         	 double tcost; // total path cost
         	 double slen; // length of shape
         	 double palong; // total path cost
			 long frid;  // edge from id
			 long toid;  // edge to id
			 long orig;  // edge to id
			 short tway;  // two way or one way indicator
        	int stopon; //stop fo ons 
        	int stopoff; //stop for offs
        	int lts; //level of traffic stress
        	string efc; // edge feature class pointer
        	string orfc; //Origin feature class pointer
        	string ornm; //Origin Name pointer
        	string enote; //Edge Note pointer
};

void edgev::serialize(ofstream& pEdge)
{

 pEdge.write(reinterpret_cast<char *>(&id), sizeof(id));
 pEdge.write(reinterpret_cast<char *>(&esid), sizeof(esid));
 pEdge.write(reinterpret_cast<char *>(&evp), sizeof(evp));
 pEdge.write(reinterpret_cast<char *>(&lbl), sizeof(lbl));
 pEdge.write(reinterpret_cast<char *>(&toStop), sizeof(toStop));
 pEdge.write(reinterpret_cast<char *>(&dirn),sizeof(dirn));
 pEdge.write(reinterpret_cast<char *>(&ecost), sizeof(ecost));
 pEdge.write(reinterpret_cast<char *>(&scost), sizeof(scost));
 pEdge.write(reinterpret_cast<char *>(&tcost), sizeof(tcost));
 pEdge.write(reinterpret_cast<char *>(&slen), sizeof(slen));
 pEdge.write(reinterpret_cast<char *>(&palong),sizeof(palong));
 pEdge.write(reinterpret_cast<char *>(&frid),sizeof(frid));
 pEdge.write(reinterpret_cast<char *>(&toid),sizeof(toid));
 pEdge.write(reinterpret_cast<char *>(&orig),sizeof(orig));
 pEdge.write(reinterpret_cast<char *>(&tway),sizeof(tway));
 pEdge.write(reinterpret_cast<char *>(&stopon),sizeof(stopon));
 pEdge.write(reinterpret_cast<char *>(&stopoff),sizeof(stopoff));
 pEdge.write(reinterpret_cast<char *>(&lts),sizeof(lts));

 size_t sizet=efc.size();// store efc's length
 pEdge.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pEdge.write(efc.c_str(), sizet+1); // write final '\0' too
 sizet=orfc.size();// store orfc's length
 pEdge.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pEdge.write(orfc.c_str(), sizet+1); // write final '\0' too
 sizet=ornm.size();// store Origin Name's length
 pEdge.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pEdge.write(ornm.c_str(), sizet+1); // write final '\0' too
 sizet=enote.size();// store pacid's length
 pEdge.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pEdge.write(enote.c_str(), sizet+1); // write final '\0' too

}

void edgev::deserialize(ifstream& pEdge)
{
 int len=0;
 char *p=0;
 pEdge.read(reinterpret_cast<char *>(&id), sizeof(id));
 pEdge.read(reinterpret_cast<char *>(&esid), sizeof(esid));
 pEdge.read(reinterpret_cast<char *>(&evp), sizeof(evp));
 pEdge.read(reinterpret_cast<char *>(&lbl), sizeof(lbl));
 pEdge.read(reinterpret_cast<char *>(&toStop), sizeof(toStop));
 pEdge.read(reinterpret_cast<char *>(&dirn), sizeof(dirn));
 pEdge.read(reinterpret_cast<char *>(&ecost), sizeof(ecost));
 pEdge.read(reinterpret_cast<char *>(&scost), sizeof(scost));
 pEdge.read(reinterpret_cast<char *>(&tcost), sizeof(tcost));
 pEdge.read(reinterpret_cast<char *>(&slen), sizeof(slen));
 pEdge.read(reinterpret_cast<char *>(&palong),sizeof(palong));
 pEdge.read(reinterpret_cast<char *>(&frid),sizeof(frid));
 pEdge.read(reinterpret_cast<char *>(&toid),sizeof(toid));
 pEdge.read(reinterpret_cast<char *>(&orig),sizeof(orig));
 pEdge.read(reinterpret_cast<char *>(&tway),sizeof(tway));
 pEdge.read(reinterpret_cast<char *>(&stopon),sizeof(stopon));
 pEdge.read(reinterpret_cast<char *>(&stopoff),sizeof(stopoff));
 pEdge.read(reinterpret_cast<char *>(&lts),sizeof(lts));

 pEdge.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 pEdge.read(p, len+1); // copy string to temp, including '\0'
 efc=p; // copy temp to data member
 delete[] p;

 pEdge.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 pEdge.read(p, len+1); // copy string to temp, including '\0'
 orfc=p; // copy temp to data member
 delete[] p;
 
 pEdge.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 pEdge.read(p, len+1); // copy string to temp, including '\0'
 ornm=p; // copy temp to data member
 delete[] p;

 pEdge.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 pEdge.read(p, len+1); // copy string to temp, including '\0'
 enote=p; // copy temp to data member
 delete[] p;

}

edgev::edgev( long id1, long esid1=0,long evp1=0, short lbl1=0, double ecost1=0,double scost1=0,
			 double tcost1=0,double slen1=0,long frid1=-1, long toid1=-1,long orig1=-1, double palong1=0,
			 short tway1=-1,short dirn1=0,short toStop1=0,int lts1=999,int stopon1=-1,int stopoff1=-1,
			 string efc1="",string orfc1="",string ornm1="",string enote1="") : id(id1), esid(esid1),evp(evp1), 
			 lbl(lbl1), ecost(ecost1), scost(scost1), tcost(tcost1),slen(slen1),frid(frid1),toid(toid1),
			 orig(orig1), palong(palong),tway(tway1),dirn(dirn1),toStop(toStop1),lts(lts1),stopon(stopon1),
			 stopoff(stopoff1),efc(efc1),orfc(orfc1),ornm(ornm1),enote(enote1)  
{
  // cout << "ev[" << id << "]" << endl;
}

bool edgev::operator==(edgev ev2)
 {
   if(id!=ev2.id)
      return false;
   if(esid!=ev2.esid)
      return false;
   if(frid!=ev2.frid)
      return false;
   if(toid!=ev2.toid)
      return false;
   if(ecost!=ev2.ecost)
      return false;
   if(toStop!=ev2.toStop)
      return false;
   return true;
 }

void edgev::serializetext(ofstream& pEdge)
{
 pEdge <<id<<"\t"<<esid<<"\t"<<evp<<"\t"<<lbl<<"\t"<<toStop<<"\t"<<dirn<<"\t"<<
	 ecost<<"\t"<<scost<<"\t"<<tcost<<"\t"<<slen<<"\t"<<palong<<"\t"<<frid<<"\t"<<toid<<"\t"<<
 orig<<"\t"<<tway<<"\t"<<efc<<"\t"<<stopoff<<"\t"<<stopon<<"\t"<<orfc<<"\t"<<ornm<<"\t"<<enote<<endl;

}

void edgev::serializetexthdr(ofstream& pEdge)
{
 pEdge <<"id"<<"\t"<<"esid"<<"\t"<<"evp"<<"\t"<<"lbl"<<"\t"<<"toStop"<<"\t"<<"dirn"<<"\t"
	 <<"ecost"<<"\t"<<"scost"<<"\t"<<"tcost"<<"\t"<<"length"<<"\t"<<"palong"<<"\t"<<"frid"<<"\t"
	 <<"toid"<<"\t"<<"orig"<<"\t"<<"tway"<<"\t"<<"efc"<<"\t"<<"stopoff"<<"\t"<<"stopon"<<"\t"
	 <<"orfc"<<"\t"<<"ornm"<<"\t"<<"enote"<<endl;

}

    typedef vector<edgev> edgevect;


class inputfilelist {
public:

//scenetable, globalcosttable, periodtable, landusecodetable, stoptable, triptable,edgetable, vertextable, parceltable

	inputfilelist(string rtename , string dirname , string tblgcost1, string tblperiod1, string tbllanduse1, 
		string tblstop1, string tbltrip1,string tbledge1, string tblvertex1, string tblparcel1,
		string tblscene1, string filegcost1, string fileperiod1, string filelanduse1, 
		string filestop1,string filetrip1, string fileedge1, string filevertex1, string fileparcel1,
		string filescene1, bool historic1 , bool euclid1, bool pdRun1, bool parcel1 );

	inputfilelist( string rtename,string dirname , string tblgcost1, // global cost table
			 string tblperiod1, // period headway table
			  string tbllanduse1,  // landuse table 
			 string tblstop1,  // stop table
			 string tbltrip1,  // stop table
	         string tbledge1, // edge table
			 string tblvertex1, // vertex table
			string tblparcel1, // parcel table
			string tblscene1,  // Scene table
			 string inputDB1,  // input database file
			 string resultDB1,  // result database file
			bool historic1, // If Historic run is needed
			bool euclid1,  // if euclidean metric is to be used.
			bool pdRun1,  // Period Run
			bool parcel1,  // Parcel data
			int srid1,  // Spatial reference System used.
			int lts1,  // Spatial reference System used.
//			long buffer1,  // Buffer distance to be used as selection
			string filegcost1, // global cost file
			 string fileperiod1, // period headway file
			  string filelanduse1,  // landuse file 
			 string filestop1,  // stop file
			 string filetrip1,  // stop file
	         string fileedge1, // edge file
			 string filevertex1, // vertex file
			 string fileparcel1, // parcel file
			 string filescene1  // Scene file
			 );

	void set_rtename (string rtename1) {rtename=rtename1;}

	string get_rtename () {return rtename;}  // Route Name 

	void set_dirname (string dirname1) {dirname=dirname1;}

	string get_dirname () {return dirname;}  // route direction name 

	void set_filegcost (string filegcost1) {filegcost=filegcost1;}

	string get_filegcost () {return filegcost;}  // file for global cost 

	void set_fileperiod (string fileperiod1) {fileperiod=fileperiod1;}

	string get_fileperiod () {return fileperiod;}  // file for period data

	void set_filelanduse (string filelanduse1) {filelanduse=filelanduse1;}

	string get_filelanduse () {return filelanduse;}  // file for landuse coefficients 

	void set_filestop (string filestop1) {filestop=filestop1;}

	string get_filestop () {return filestop;}  // file for stop data 

	void set_filetrip (string filetrip1) {tbltrip=filetrip1;}

	string get_filetrip () {return filetrip;}  // trip table file  

	void set_fileedge (string fileedge1) {fileedge=fileedge1;}

	string get_fileedge () {return fileedge;}  // file for edge data 

	void set_filevertex (string filevertex1) {filevertex=filevertex1;}

	string get_filevertex () {return filevertex;}  // file for vertex data 

	void set_fileparcel (string fileparcel1) {fileparcel=fileparcel1;}

	string get_fileparcel () {return fileparcel;}  // file for parcel data 

	void set_tblgcost (string tblgcost1) {tblgcost=tblgcost1;}

	string get_tblgcost () {return tblgcost;}  // table global cost 

	void set_tblperiod (string tblperiod1) {tblperiod=tblperiod1;}

	string get_tblperiod () {return tblperiod;}  // table period data

	void set_tbllanduse (string tbllanduse1) {tbllanduse=tbllanduse1;}

	string get_tbllanduse () {return tbllanduse;}  // table landuse coefficients 

	void set_tblstop (string tblstop1) {tblstop=tblstop1;}

	string get_tblstop () {return tblstop;}  // table stop data 

	void set_tbltrip (string tbltrip1) {tbltrip=tbltrip1;}

	string get_tbltrip () {return tbltrip;}  // trip table data 

	void set_tbledge (string tbledge1) {tbledge=tbledge1;}

	string get_tbledge () {return tbledge;}  // table edge data 

	void set_tblvertex (string tblvertex1) {tblvertex=tblvertex1;}

	string get_tblvertex () {return tblvertex;}  // table vertex data 

	void set_tblparcel (string tblparcel1) {tblparcel=tblparcel1;}

	string get_tblparcel () {return tblparcel;}  // table parcel data 

	void set_inputDB (string inputDB1) {inputDB=inputDB1;}

	string get_inputDB () {return inputDB;}  // input database parcel data 

	void set_resultDB (string resultDB1) {resultDB=resultDB1;}

	string get_resultDB () {return resultDB;}  // result database parcel data 

	void set_tblscene (string tblscene1) {tblscene=tblscene1;}

	string get_tblscene () {return tblscene;}  // table scene data 

	void set_filescene (string filescene1) {filescene=filescene1;}

	string get_filescene () {return filescene;}  // file scene data 

	void set_historic (bool historic1) {historic=historic1;}

	bool get_historic () {return historic;}  // if historic run needs to be run first 

	void set_euclid (bool euclid1) {euclid=euclid1;}

	bool get_euclid () {return euclid;}  // if euclidean geometry is used for the Voronoi 

	bool get_parcel () {return parcel;}  // if parcel input data is povided 

	void set_parcel (bool parcel1) {parcel=parcel1;}

	bool get_pdRun () {return pdRun;}  // if historic run needs to be run first 

	void set_pdRun (bool pdRun1) {pdRun=pdRun1;}

	void set_srid (int srid1) {srid=srid1;}

	int get_srid () {return srid;}  // Spatial reference System ID as per EPSG 

	//void set_buffer (long buffer1) {buffer=buffer1;} // buffer distance used to select data  

	//long get_buffer () {return buffer;}  // buffer distance used to select data  

	void set_lts (int lts1) {lts=lts1;}

	int get_lts () {return lts;}  // Level of traffic Stress  

    void show_inputfiles(ostream& out)
    { 
		out <<endl
			<< "Route name " <<"\t" << get_rtename()<<endl
			<< "Direction " <<"\t" << get_dirname()<<endl
			<< "Input SQLite Database " <<"\t" << get_inputDB()<<endl
			<< "Result output SQLite Database " <<"\t" << get_resultDB()<<endl
			<< "Global Cost Table " <<"\t" << get_tblgcost()<<endl
			<< "Period Table " <<"\t" << get_tblperiod()<<endl
			<< "Land Use Code Table " <<"\t" << get_tbllanduse()<<endl
			<< "Stop Table " <<"\t" << get_tblstop()<<endl
			<< "Trip Table " <<"\t" << get_tbltrip()<<endl
			<< "Edge Table " <<"\t" << get_tbledge()<<endl
			<< "Vertex Table " <<"\t" << get_tblvertex()<<endl
			<< "Parcel Table " <<"\t" << get_tblparcel()<<endl
			<< "Scenario Table " <<"\t" << get_tblscene()<<endl
			<< "Input Spatial Database File Name  " <<"\t" << get_inputDB()<<endl
			<< "Result Spatial Database File Name  " <<"\t" << get_resultDB()<<endl
			<< "Global Cost File  " <<"\t" << get_filegcost()<<endl
			<< "Period File  " <<"\t" << get_fileperiod()<<endl
			<< "Land Use Code File  " <<"\t" << get_filelanduse()<<endl
			<< "Stop File  " <<"\t" << get_filestop()<<endl
			<< "trip File  " <<"\t" << get_filetrip()<<endl
			<< "Edge File  " <<"\t" << get_fileedge()<<endl
			<< "Vertex File  " <<"\t" << get_filevertex()<<endl
			<< "Parcel File  " <<"\t" << get_fileparcel()<<endl
			<< "Parcel Data  " <<"\t" << get_parcel()<<endl
			<< "Scenario File  " <<"\t" << get_filescene()<<endl
			<< "Historic Run is Needed " <<"\t" << get_historic()<<endl
			<< "Use Euclidian Geometry " <<"\t" << get_euclid()<<endl
			<< "Period Summary Run " <<"\t" << get_pdRun()<<endl
			<< "Use SRID " <<"\t" << get_srid()<<endl<<endl
			<< "Use LTS " <<"\t" << get_lts()<<endl<<endl;
	};


	virtual ~inputfilelist(){
		//	cout << "input file list Object is deleted! "<<endl;
	}

protected:
			string rtename; // route name
			string dirname; // route name
			string tblgcost; // global cost table
			string tblperiod; // period headway table
			string tbllanduse;  // landuse table 
			string tblstop;  // stop table
			string tbltrip;  // trip table
	        string tbledge; // edge table
			string tblvertex; // vertex table
			string tblparcel; // parcel table	         
			string tblscene; // scenario table
			string filegcost; // global cost file
			string fileperiod; // period headway file
			string filelanduse;  // landuse file 
			string filestop;  // stop file
			string filetrip;  // trip file
	        string fileedge; // edge file
			string filevertex; // vertex file
			string fileparcel; // parcel file	         
			string inputDB; // input database for all network and route tables	         
			string resultDB; // output results database for evaluation and analysis tables	         
			string filescene; // scenario file
			bool historic; // if historic run is to run first
			bool euclid; // if euclidean geometry from stop to parcel is to be used
			bool pdRun; // if Period Summary Run is only done
			bool parcel; // if Parcel data is provided
			int srid; // geometry Spatial Reference ID from EPSG
			int lts; // Level of traffic Stress for the Street network 
		//	long buffer; // buffer distance around route used to select data  
};

inputfilelist::inputfilelist(string rtename1 ="12",string dirname1 ="N" ,string tblgcost1="",string tblperiod1="", 
		string tbllanduse1="", string tblstop1="", string tbltrip1="",string tbledge1="", string tblvertex1="",
		 string tblparcel1="", string tblscene1="",string inputDB1="", string resultDB1="", 
		bool historic1 = true, bool euclid1 = false , bool pdRun1 = false, bool parcel1 = false,int srid1=2232, 
		int lts1=999,string filegcost1="", string fileperiod1="", string filelanduse1="", string filestop1="",
		string filetrip1="", string fileedge1="", string filevertex1="", string fileparcel1="", string filescene1="") 
{
	rtename=rtename1;
	dirname=dirname1;
	tblgcost=tblgcost1;
	tblperiod=tblperiod1;
	tbllanduse=tbllanduse1; 
	tblstop=tblstop1;
	tbltrip=tbltrip1;
	tbledge=tbledge1;
	tblvertex=tblvertex1;
	tblparcel=tblparcel1;
	tblscene=tblscene1;
	inputDB=inputDB1;
	resultDB=resultDB1;
	historic= historic1; 
	euclid=euclid1;
	pdRun=pdRun1;
	parcel=parcel1;
	srid=srid1;
	lts=lts1;
	//buffer=buffer1;
	filegcost=filegcost1;
	fileperiod=fileperiod1;
	filelanduse=filelanduse1; 
	filestop=filestop1;
	filetrip=filetrip1;
	fileedge=fileedge1;
	filevertex=filevertex1;
	fileparcel=fileparcel1;
	filescene=filescene1;
}



class globcost {
public:

//COSTWALK,COSTRIDE,UNITONTM,UNITOFFTM,MAXWLKDIST,PROPENSITY,FILESTEM,NOPERIODS,
//WALKSPD,FILEPATH

	globcost(double walkcost1,float unitontm1, float unitofftm1, short nopds1,short maxskip1,
		float maxwalkdist1,double ridecost1, float propensity1,float walkspd1,int dpdimension1, 
		double unitconv1);

	globcost(double walkcost1,  // vertex id
	         short nopds1, // period
	         short maxskip1, // Maximum no of stops to skip (m-1)
         	 float unitontm1,  // unit boarding time
             float unitofftm1,  // unit alighting time
			 float maxwalkdist1, // total cost
			 double ridecost1,  // ride cost 
			 float propensity1,  // propensity ratio
	         float walkspd1, // walk speed 
	         string  filestem1, // stem for file naming 
	         string  filepath1,
			 int dpdimension1,
			 double unitconv1) ; 
//: id(gc.id)
	globcost(const globcost& gc )  {
    //std::cout << "c[" << id << "]" << endl;
//    ++copycons;
  }
  globcost& operator=(const globcost& gc) {
    //cout << "(" << id << ")=[" << gc.id << "]" << endl;
    id = gc.id;
//    ++assign;
    return *this;
  }
  friend bool operator<(const globcost& gb, const globcost& gc ) {
    return gb.id < gc.id;
  }
  friend bool operator==(const globcost& gb,const globcost& gc) {
    return gb.id == gc.id;
  }
 virtual ~globcost() {
    //cout << "~[" << id << "]" << endl;
//    ++destroy;
  }


  friend ostream& operator<<(ostream& os, const globcost& gc) {
    return os << gc.id;
  }
  friend class gcostReport;


	void set_walkcost (double walkcost1) {walkcost=walkcost1;}

	double get_walkcost () {return walkcost;}

	void set_unitontm (float unitontm1) {unitontm=unitontm1;}

	float get_unitontm () {return unitontm;}  // unit on time 

	void set_unitofftm (float unitofftm1) {unitofftm=unitofftm1;}

	float get_unitofftm () {return unitofftm;}  // unit off time 

	void set_maxwalkdist (float maxwalkdist1) {maxwalkdist=maxwalkdist1;}

	float get_maxwalkdist () {return maxwalkdist;}  // unit off time 

	void set_propensity (float propensity1) {propensity=propensity1;}

	float get_propensity () {return propensity;}  // unit off time 

	void set_walkspd (float walkspd1) {walkspd=walkspd1;}

	float get_walkspd () {return walkspd;}  // unit off time 

	void set_filestem (string filestem1) {filestem=filestem1;}

	string get_filestem () {return filestem;}  // unit off time 

	void set_filepath (string filepath1) {filepath=filepath1;}

	string get_filepath () {return filepath;}  // unit off time 

	void set_nopds (short nopds1) {nopds=nopds1;}

	short get_nopds () {return nopds;}

	void set_maxskip (short maxskip1) {maxskip=maxskip1;}

	short get_maxskip () {return maxskip;}

	void set_ridecost (double ridecost1) {ridecost=ridecost1;}

	double get_ridecost () {return ridecost;}

	void set_dpdimension (int dpdimension1) {dpdimension=dpdimension1;}

	short get_dpdimension () {return dpdimension;}

	double get_unitconv () {return unitconv;}  // unit conversion factor 

	void set_unitconv (double unitconv1) {unitconv=unitconv1;}

	void show_globcost(ostream& out)
    { 
		out << get_walkcost() << "\t" << get_ridecost() <<"\t" << get_unitontm() << "\t" 
			<< get_unitofftm() << "\t" << get_maxwalkdist() << "\t" << get_propensity() << 
			"\t" << get_filestem()<<"\t" << get_filepath()<<"\t" << get_walkspd()<< 
			"\t" << get_nopds()<<"\t" << get_maxskip()<<"\t" <<get_dpdimension()<<"\t"<<get_unitconv()<<endl;
	};

    void show_globcosthdr(ostream& out)
    { 
		out << "walkcost" << "\t" << "ridecost" <<"\t" << "unitontm" << "\t" 
			<< "unitofftm" << "\t" << "maxwalkdist" << "\t" << "propensity" << 
			"\t" << "filestem"<<"\t" << "filepath"<<"\t" << "walkspd"<< 
			"\t" << "nopds"<<"\t" << "maxskip"<<"\t" <<"DPDimension"<<"\t"<<"unitconv"<<endl;
	};



private:
//	static long  create, assign, copycons, destroy;
	static long id;
	double walkcost;  // vertex id
	short nopds; // no of periods
    short maxskip; // Maximum no of stops to skip (m-1)
    float unitontm;  // unit boarding time
    float unitofftm;  // unit alighting time
	float maxwalkdist; // total cost
	double ridecost;  // ride cost 
	float propensity;  // propensity ratio
	float walkspd; // walk speed 
	string  filestem; // stem for file naming 
	string  filepath; // path name for output files
	int dpdimension;
	double unitconv;
};

globcost::globcost(double walkcost1=10,float unitontm1=2, float unitofftm1=1.5, short nopds1=2,	short maxskip1=2,  
		float maxwalkdist1=800,double ridecost1=5, float propensity1=0.3333,float walkspd1=1.2,int dpdimension1=5,
		double unitconv1=1.0)
		{
//    cout << "d[" << id << "]" << endl;
		id++; 
		walkcost = walkcost1;
		unitontm = unitontm1;
		unitofftm= unitofftm1;
		ridecost = ridecost1;
		nopds = nopds1;
		maxwalkdist= maxwalkdist1;
		propensity= propensity1;
		walkspd= walkspd1;
	    maxskip=maxskip1; // Maximum no of stops to skip (m-1)
		dpdimension=dpdimension1;
		unitconv=unitconv1;
}

class pdhdway {
public:

//FLDPERD,FLDBEGT,FLDHDWAY,FLDPDLEN,COSTOPER,Include

	pdhdway(short pdId1, string pdkey1,int numTrips1,
		double begtm1, float pdlen1,double hdway1,float opercost1, short include1);

	pdhdway( short pdId1, // period
			 string pdkey1, // string key 	
			 int numTrip1, // number of trips 	
			 float begtm1, // total cost in 24 hour clock
			 float pdlen1,  // pdlen in hours
			 double hdway1,  // headway in minutes 
	         float opercost1, // cost of operating a bus per hour
	         short include1 // include int he analysis
			 );

	void set_begtm (float begtm1) {begtm=begtm1;}

	float get_begtm () {return begtm;}  // unit off time 

	void set_pdlen (float pdlen1) {pdlen=pdlen1;}

	float get_pdlen () {return pdlen;}  // unit off time 

	void set_opercost (float opercost1) {opercost=opercost1;}

	float get_opercost () {return opercost;}  // unit off time 

	void set_pdId (short pdId1) {pdId=pdId1;}

	short get_pdId () {return pdId;}

	void set_pdKey (string pdKey1) {pdKey=pdKey1;}

	string get_pdKey () {return pdKey;}

	void set_numTrips (int numTrips1) {numTrips=numTrips1;}

	short get_numTrips () {return numTrips;}

	void set_hdway (double hdway1) {hdway=hdway1;}

	double get_hdway () {return hdway;}

	void set_include (short include1) {include=include1;}

	short get_include () {return include;}


	
	void show_pdhdway(ostream& out)
    { 
		out << get_pdId()<< "\t" << get_pdKey()<<"\t" << get_begtm()<< "\t" << get_pdlen() 
			<< "\t" << get_numTrips()<< "\t" << get_hdway() << "\t" << get_opercost()<< "\t" << get_include()<<endl;
	};

    void show_pdhdwayhdr(ostream& out)
    { 
		out << "pdId"<< "\t" << "PdKey"<<"\t" << "begtm"<<"\t"<< "pdlen" <<"\t"<< "NumTrips"
			<< "\t"  "hdway" << "\t" << "opercost"<< "\t" << "Include"<<endl;
	};

    void calc_hdway(float tmFact=60.0)
    { 
		hdway = pdlen*tmFact/numTrips;
	};


	virtual ~pdhdway(){
		//	cout << "Vertex Object is deleted! "<<endl;
			//delete[] pv;
//		    plist.clear();
	}

protected:
	         short pdId; // period Id
			 string pdKey; // pdKey  
			 int numTrips;  // number of trips in a period used to calculate headway if needed  
			 float begtm; // Beginning of Period time
			 double hdway;  // headway 
			 float pdlen;  // period length 
	         float opercost; // Operating Cost 
	         short include; // include in trip run 
};

pdhdway::pdhdway(short pdId1=3,string pdKey1="", int numTrips1=0, float begtm1=3.5,
				 float pdlen1=2.5,double hdway1=5, float opercost1=80.0,short include1=1) 
{
	pdhdway::hdway = hdway1;
	pdhdway::pdId = pdId1;
	pdhdway::pdKey = pdKey1;
	pdhdway::numTrips = numTrips1;
	pdhdway::begtm= begtm1;
	pdhdway::pdlen= pdlen1;
	pdhdway::opercost= opercost1;
	pdhdway::include= include1;
}

class lucodes {
public:

//PTYPE,DESC_,LU,,KEYFLD,ONCOEF,OFFCOEF


	lucodes(string pType1,string Desc1,string LUC1,string KeyProp1, float OnCoeff1,float OffCoeff1,int PdOveRide1 );

	lucodes( string pType1, // Property Type
			 string LUC1, // Land Use Code
			 string KeyProp1,  // Key Property of parcel to be used with coeff's 
			 float OnCoeff1,  // On Coeff
	         float OffCoeff1, // Off Coeff
	         int PdOveRide1 // Overide Period Prod. Attractio AM/PM zeroing of residential and commercial Land Use 
			 );

	void set_pType (string pType1) {pType=pType1;}

	string get_pType () {return pType;}  // Property Type 

	void set_LUC (string LUC1) {LUC=LUC1;}

	string get_LUC () {return LUC;}  // Land Use Code 

	void set_KeyProp (string KeyProp1) {KeyProp=KeyProp1;}

	string get_KeyProp () {return KeyProp;}  // Key Parcel Property used for prod/attract 

	void set_Desc (string Desc1) {Desc=Desc1;}

	string get_Desc () {return Desc;}  // Property Type description

	void set_OnCoeff (float OnCoeff1) {OnCoeff=OnCoeff1;}

	float get_OnCoeff () {return OnCoeff;}  // unit off time 

	void set_OffCoeff (float OffCoeff1) {OffCoeff=OffCoeff1;}

	float get_OffCoeff () {return OffCoeff;}  // unit off time 

	void set_PdOveRide (int PdOveRide1) {PdOveRide = PdOveRide1;} // Period OveRide

	int get_PdOveRide () {return PdOveRide;}  // Period OveRide

	void show_lucodes(ostream& out)
    { 
		out << get_pType() <<"\t" << get_Desc() << "\t" << get_LUC() << "\t" << 
			get_KeyProp()<< "\t" << get_OnCoeff()<< "\t" << get_OffCoeff()<<endl;
	};

    void show_lucodeshdr(ostream& out)
    { 
		out << "pType" <<"\t" << "Desc" << "\t" << "LUC" << "\t" << 
			"KeyProp"<< "\t" << "OnCoeff"<< "\t" << "OffCoeff"<<"PdOveRide"<<endl;
	};


	virtual ~lucodes(){
		//	cout << "Vertex Object is deleted! "<<endl;
			//delete[] pv;
//		    plist.clear();
	}

protected:
			string pType; // Property Type
			string Desc; // Property Type Description
			string LUC; // Land Use Code
			string KeyProp;  // General Property of parcel to be used with coeff's 
			float OnCoeff;  // On Coeff
	        float OffCoeff; // Off Coeff
			int PdOveRide; // Overide Period Prod. Attractio AM/PM zeroing of residential and commercial Land Use 
};

lucodes::lucodes(string pType1="0",string Desc1="",string LUC1="",
				 string KeyProp1="", float OnCoeff1=0.0,float OffCoeff1=0.0, int PdOveRide1=0) 
{
	lucodes::pType = pType1;
	lucodes::Desc = Desc1;
	lucodes::LUC = LUC1;
	lucodes::KeyProp= KeyProp1;
	lucodes::OnCoeff= OnCoeff1;
	lucodes::OffCoeff= OffCoeff1;	
	lucodes::PdOveRide= PdOveRide1;
}




#ifndef ORIGINVX_H //Walk, Ride, Operating cost
#define ORIGINVX_H

class originvx {
public:
    
//OBJECTID,ORDER,LABEL,STOPNAME,CUM_DIST,CUM_TIME,UNDCRTM,CUM_TIME_P,HISTORIC,INBOUND,HISTONS,HISTOFFS,	
//HISTDEPVOL,ONS,OFFS,DEPVOL,PROBSTOP,DEPDELAY,ARRDELAY,DELAY,PVAL,AVAL,Include,CUM_TIME2	
//WkTmOns,WkTmOffs,External_,pid,fid,OperCost,TCost,External,Eliminated,StopLoc
    
	originvx(long pid1,int pidp, short lbl1, double cost1, long fid1, long orig1);

	originvx(long pid1,int pidp, short lbl1, double cost1, long fid1, long orig1, double tcost1);

//	void insert_plist (list<originvx*> plist1, originvx* pxp1) { plist1.insert(pxp1); plist = plist1;}

//	std::list<originvx*> get_plist () {return plist;}

	void set_id (long pid1) {pid=pid1;}

	long get_id () {return pid;}

	void set_idp (int pidp1) {pidp=pidp1;}

	int get_idp () {return pidp;}  // predecessor

	void set_lbl (short lbl1) {lbl=lbl1;}

	short get_lbl () {return lbl;}

	void set_cost (double cost1) {cost=cost1;}

	double get_cost () {return cost;}

	void set_tcost (double tcost1) {tcost=tcost1;}

	double get_tcost () {return tcost;}

	void set_fid (long fid1) {fid=fid1;}

	long get_fid () {return fid;}

	void set_orig (long orig1) {orig=orig1;}

	long get_orig () {return orig;}

    void show_vertof(ostream& out)
    { 
		out << get_id() << "\t" << get_idp() << "\t" << get_orig() << "\t" 
		  << get_cost() << "\t" << get_tcost() << "\t" << get_lbl()<< get_fid() << endl;
	};

    void show_verthdr(ostream& out)
    { 
		out  << "vid" << "\t" << "pid" << "\t" << "orig" << "\t" 
		 << "vcost" << "\t" << "tcost" << "\t" << "lbl"<<"\t"<<"fid" << endl;
	};

	void show_vertex(void)
    { 
		cout << "Vertex Id : " << get_id() << endl;
		cout << "Pred Id: " << get_fid() << endl;
		cout << "Vertex Feature Id: " << get_fid() << endl;
		cout << "Vertex Origin : " << get_orig() << endl;
		cout << "Vertex Cost : " << get_cost() << endl;
		cout << "Vertex Total Cost : " << get_tcost() << endl;
		cout << "Vertex lbl : " << get_lbl() << endl;
	};

void originvx::serialize(ofstream& pVx);
void originvx::deserialize(ifstream& pVx);


	virtual ~originvx(){
		//	cout << "Vertex Object is deleted! "<<endl;
			//delete[] pv;
//		    plist.clear();
	}

protected:
	unsigned long pid;  // vertex id
	         int pidp; // predecessor
//			 std::list<originvx*> plist; // predecessor list with ties
	         short lbl; // label
         	 double cost;  // cost at vertex
         	 double tcost; // total cost
			 long fid;  //  
			 long orig;  // Origin for vertex least cost path
};

void originvx::serialize(ofstream& pVx)
{
 pVx.write(reinterpret_cast<char *>(&pid), sizeof(pid));
 pVx.write(reinterpret_cast<char *>(&pidp), sizeof(pidp));
 pVx.write(reinterpret_cast<char *>(&lbl), sizeof(lbl));
 pVx.write(reinterpret_cast<char *>(&cost), sizeof(cost));
 pVx.write(reinterpret_cast<char *>(&tcost), sizeof(tcost));
 pVx.write(reinterpret_cast<char *>(&fid), sizeof(fid));
 pVx.write(reinterpret_cast<char *>(&orig),sizeof(orig));
}

void originvx::deserialize(ifstream& pVx)
{
 pVx.read(reinterpret_cast<char *>(&pid), sizeof(pid));
 pVx.read(reinterpret_cast<char *>(&pidp), sizeof(pidp));
 pVx.read(reinterpret_cast<char *>(&lbl), sizeof(lbl));
 pVx.read(reinterpret_cast<char *>(&cost), sizeof(cost));
 pVx.read(reinterpret_cast<char *>(&tcost), sizeof(tcost));
 pVx.read(reinterpret_cast<char *>(&fid), sizeof(fid));
 pVx.read(reinterpret_cast<char *>(&orig),sizeof(orig));

}

originvx::originvx(long pid1=0, int pidp=0, short lbl1=0, double cost1=100000,  long fid1=0, long orig1=0,double tcost1=-1) 
{
	originvx::pid = pid1;
	originvx::lbl = lbl1;
	originvx::cost= cost1;
	originvx::fid = fid1;
	originvx::orig = orig1;
	originvx::tcost= tcost1;
}
#endif // _ORIGINVX_H_ 


class tstop {
public:
//OBJECTID,ORDER,LABEL,STOPNAME,CUM_DIST,CUM_TIME,UNDCRTM,CUM_TIME_P,HISTORIC,INBOUND,HISTONS,HISTOFFS,	
//HISTDEPVOL,ONS,OFFS,DEPVOL,PROBSTOP,DEPDELAY,ARRDELAY,DELAY,PVAL,AVAL,Include,CUM_TIME2	
//WkTmOns,WkTmOffs,External_,WalkCost,RideCost,OperCost,TCost,External,Eliminated,StopLoc,XCoord,YCOORD
    
	tstop (long id1,long tripId1,string StopLbl1,long Stopidp1, long Stopids1,int StOrdr1,  
			string StopName1,double CumDist1, double CRdTm1, double undCRdTm1, double CRdTmC1, 
			short blnHist1, short blnInbd1,short blnIncl1, short blnExtr1,short blnElim1, double HistOns1, 
			 double HistOffs1, double HistDepVol1, double depDelay1, double arrDelay1, 
			 double probStoph1,double probStop1, double Ons1, double Offs1, double DepVol1,  
			 double dwlDelay1, double rideDelay1, double PVal1,double AVal1,double CRdTmE1, 
			 double hWkTmOns1, double hWkTmOffs1,double WkTmOns1, double WkTmOffs1, 
			 double WalkCost1,double RideCost1, double OperCost1,double TCost1,
			 double xc1,double yc1,double zc1,double nearDist1);

	tstop(long id1, string StopLbl1, long Stopidp1,long Stopids1, int StOrdr1,  string StopName1, double CumDist1, double CRdTm1, 
		double undCRdTm1, double CRdTmC1, short blnHist1,short blnInbd1,double blnIncl1,short blnExtr1, double HistOns1, 
		double HistOffs1, double HistDepVol1,long tripId1 =-1, double Ons1=0, double Offs1=0, double DepVol1=0, double probStop1=0, 
		double depDelay1=0,	double arrDelay1=0, double dwlDelay1=0, double rideDelay1=0, double PVal1=0,double AVal1=0,
		double CRdTmE1=0, double WkTmOns1=0, double WkTmOffs1=0, double WalkCost1=0,double RideCost1=0,double OperCost1=0,
		double TCost1=0,short blnElm1=0,long Edgeid1=-1, double xc1=0, double yc1=0 , double zc1=0,double nearDist1=-1 );

	tstop(long id1, string StopLbl1, long Stopidp1,long Stopids1, int StOrdr1, string StopName1, double CumDist1,  
		double CRdTm1,double undCRdTm1, double CRdTmC1,long tripId1 =-1, short blnHist1=0,short blnInbd1=0,short blnIncl1=0,short blnExtr1=-1,
		 short blnElim1=0,double HistOns1=0, double HistOffs1=0, double HistDepVol1=0,double xc1=0, double yc1=0, double zc1=0 );
	tstop() : id(0), tripid(-1), StopLbl(""), Stopidp(-1),Stopids(-1), StOrdr(-1),  StopName(""), 
		CumDist(0), CRdTm(0),undCRdTm(0), CRdTmC(0),blnHiStop(1),blnInbd(1),blnIncl(1),
		blnElim(0),blnExtr(0),HistOns(0),HistOffs(0),HistDepVol(0),Ons(0),Offs(0),DepVol(0),
		PVal(0),AVal(0),probStoph(0),probStop(0),depDelay(0), arrDelay(0), dwlDelay(0),
		rideDelay(0),CRdTmE(0),hWkTmOns(0),hWkTmOffs(0), WkTmOns(0), WkTmOffs(0), WalkCost(0), 
		RideCost(0),OperCost(0),TCost(0),Edgeid(-1),palong(0),lbl(0),xc(0),yc(0),zc(0),nearDist(0) {};
	tstop(const tstop& ts )  {
	tstop::id = ts.id;
	tstop::tripid = ts.tripid;
	tstop::StopLbl= ts.StopLbl;
	tstop::Stopidp = ts.Stopidp;
	tstop::Stopids = ts.Stopids;
	tstop::Edgeid = ts.Edgeid;
	tstop::StopName = ts.StopName;
	tstop::StOrdr = ts.StOrdr;
	tstop::blnHiStop= ts.blnHiStop; 
	tstop::blnInbd= ts.blnInbd; 
	tstop::blnExtr= ts.blnExtr; 
	tstop::blnIncl= ts.blnIncl;
	tstop::blnElim= ts.blnElim;
	tstop::palong= ts.palong;
	tstop::lbl= ts.lbl;
	tstop::CumDist= ts.CumDist;
	tstop::CRdTm= ts.CRdTm; 
	tstop::undCRdTm= ts.undCRdTm; 
	tstop::CRdTmC= ts.CRdTmC; 
	tstop::HistOns= ts.HistOns; 
	tstop::HistOffs= ts.HistOffs; 
	tstop::HistDepVol= ts.HistDepVol; 
	tstop::probStoph= ts.probStoph; 
	tstop::Ons= ts.Ons;
	tstop::Offs = ts.Offs;
	tstop::DepVol = ts.DepVol;
	tstop::probStop= ts.probStop; 
	tstop::depDelay= ts.depDelay; 
	tstop::arrDelay= ts.arrDelay; 
	tstop::dwlDelay= ts.dwlDelay; 
	tstop::rideDelay= ts.rideDelay;
	tstop::PVal= ts.PVal; 
	tstop::AVal= ts.AVal; 
	tstop::CRdTmE= ts.CRdTmE; 
	tstop::hWkTmOns= ts.hWkTmOns;
	tstop::hWkTmOffs= ts.hWkTmOffs; 
	tstop::WkTmOns= ts.WkTmOns;
	tstop::WkTmOffs= ts.WkTmOffs; 
	tstop::WalkCost= ts.WalkCost; 
	tstop::RideCost= ts.RideCost; 
	tstop::OperCost= ts.OperCost; 
	tstop::TCost= ts.TCost; 
	tstop::xc = ts.xc; 
	tstop::yc = ts.yc; 
	tstop::zc = ts.zc; 
	tstop::nearDist= ts.nearDist; 
	//  cout << "Ts[" << id << "]" << endl;
	//    ++copycons;
  }
  tstop& operator=(const tstop& ts) {
    //cout << "(" << id << ")=[" << ts.id << "]" << endl;
	tstop::id = ts.id;
	tstop::tripid = ts.tripid;
	tstop::StopLbl= ts.StopLbl;
	tstop::Stopidp = ts.Stopidp;
	tstop::Stopids = ts.Stopids;
	tstop::Edgeid = ts.Edgeid;
	tstop::StopName = ts.StopName;
	tstop::StOrdr = ts.StOrdr;
	tstop::blnHiStop= ts.blnHiStop; 
	tstop::blnInbd= ts.blnInbd; 
	tstop::blnExtr= ts.blnExtr; 
	tstop::blnIncl= ts.blnIncl;
	tstop::blnElim= ts.blnElim;
	tstop::palong= ts.palong;
	tstop::lbl= ts.lbl;
	tstop::CumDist= ts.CumDist;
	tstop::CRdTm= ts.CRdTm; 
	tstop::undCRdTm= ts.undCRdTm; 
	tstop::CRdTmC= ts.CRdTmC; 
	tstop::HistOns= ts.HistOns; 
	tstop::HistOffs= ts.HistOffs; 
	tstop::HistDepVol= ts.HistDepVol; 
	tstop::probStoph= ts.probStoph; 
	tstop::Ons= ts.Ons;
	tstop::Offs = ts.Offs;
	tstop::DepVol = ts.DepVol;
	tstop::probStop= ts.probStop; 
	tstop::depDelay= ts.depDelay; 
	tstop::arrDelay= ts.arrDelay; 
	tstop::dwlDelay= ts.dwlDelay; 
	tstop::rideDelay= ts.rideDelay;
	tstop::PVal= ts.PVal; 
	tstop::AVal= ts.AVal; 
	tstop::CRdTmE= ts.CRdTmE; 
	tstop::hWkTmOns= ts.hWkTmOns;
	tstop::hWkTmOffs= ts.hWkTmOffs; 
	tstop::WkTmOns= ts.WkTmOns;
	tstop::WkTmOffs= ts.WkTmOffs; 
	tstop::WalkCost= ts.WalkCost; 
	tstop::RideCost= ts.RideCost; 
	tstop::OperCost= ts.OperCost; 
	tstop::TCost= ts.TCost; 
	tstop::xc = ts.xc; 
	tstop::yc = ts.yc; 
	tstop::zc = ts.zc; 
	tstop::nearDist= ts.nearDist; 
    return *this;
  }

	bool operator==(tstop ts);

  // Create a stop from a stop pointer:
  tstop(const tstop* ts, const int& ix) 
    : id(ts->id+ix) {
    cout << "Copied stop " << *this << " from "
         << *ts << endl;
  }

	void tsr(double undCRdTm1=0, double CRdTmC1=0, double Ons1=0, double Offs1=0, double DepVol1=0, double probStop1=0, double depDelay1=0, 
		double arrDelay1=0, double dwlDelay1=0,double rideDelay1=0, double PVal1=0,double AVal1=0,double CRdTmE1=0, double WkTmOns1=0, double WkTmOffs1=0,
		double WalkCost1=0,double RideCost1=0,double OperCost1=0,double TCost1=0);

   void cdwlDelay(double unitOnTime, double unitOffTime, double hdwy ); 
   void cdwlDelayh(double unitOnTime, double unitOffTime, double hdwy ); 
   void cprobstop( double hdwy1); 
   void cprobstoph( double hdwy1); 
   void cDepVol(double depVol1);
   void cDepVolh(double depVol1);
   double tstop::cstopDelay();
   double tstop::thruVol();
   double tstop::thruVolh();
   double segDelay(double probjm1, double depDjm1, double stoProb); 
 
   void set_id (long id1) {id=id1;}

   long get_id () {return id;}

	void settripId (long tripid1) {tripid=tripid1;}

	long tripId () {return tripid;}  // trip Id


   void set_Edgeid (long Edgeid1) {Edgeid=Edgeid1;}

	long get_Edgeid () {return Edgeid;}

	void set_posalong (double palong1) {palong=palong1;}

	double get_posalong () {return palong;}

	void set_nearDist (double neardist1) {nearDist=neardist1;}

	double get_nearDist () {return nearDist;}

	void set_Stopidp (long Stopidp1) {Stopidp=Stopidp1;}

	long get_Stopidp () {return Stopidp;}  // predecessor

	void set_Stopids (long Stopids1) {Stopids=Stopids1;}

	long get_Stopids () {return Stopids;}  // predecessor

	void set_StOrdr (int StOrdr1) {StOrdr=StOrdr1;}

	int get_StOrdr () {return StOrdr;}  // predecessor

	void set_StopLbl (string StopLbl1) {StopLbl=StopLbl1;}

	string get_StopLbl () {return StopLbl;}

	void set_StopName (string StopName1) {StopName=StopName1;}

	string get_StopName () {return StopName;}

	void set_lbl (short lbl1) {lbl=lbl1;}

	short get_lbl () {return lbl;}

	void set_blnHist (short blnHist1) {blnHiStop=blnHist1;}

	short get_blnHist () {return blnHiStop;}

	void set_blnInbd (short blnInbd1) {blnInbd=blnInbd1;}

	short get_blnInbd () {return blnInbd;}

	void set_blnIncl (short blnIncl1) {blnIncl=blnIncl1;}

	short get_blnIncl () {return blnIncl;}

	void set_blnExtr (short blnExtr1) {blnExtr=blnExtr1;}

	short get_blnExtr () {return blnExtr;}

	void set_blnElim (short blnElim1) {blnElim=blnElim1;}

	short get_blnElim () {return blnElim;}

	void set_CumDist (double CumDist1) {CumDist=CumDist1;}

	double get_CumDist () {return CumDist;}

	void set_CRdTm (double CRdTm1) {CRdTm=CRdTm1;}

	double get_CRdTm () {return CRdTm;}

	void set_HistOns (double HistOns1) {HistOns=HistOns1;}

	double get_HistOns () {return HistOns;}

	void set_HistOffs (double HistOffs1) {HistOffs=HistOffs1;}

	double get_HistOffs () {return HistOffs;}

	void set_HistDepVol (double HistDepVol1) {HistDepVol=HistDepVol1;}

	double get_HistDepVol () {return HistDepVol;}

	void set_Ons (double Ons1) {Ons=Ons1;}

	double get_Ons () {return Ons;}

	void set_Offs (double Offs1) {Offs=Offs1;}

	double get_Offs () {return Offs;}

	void set_DepVol (double DepVol1) {DepVol=DepVol1;}

	double get_DepVol () {return DepVol;}

	void set_undCRdTm (double undCRdTm1) {undCRdTm=undCRdTm1;}

	double get_undCRdTm () {return undCRdTm;}

	void set_CRdTmC (double CRdTmC1) {CRdTmC=CRdTmC1;}

	double get_CRdTmC () {return CRdTmC;}

	void set_probStop (double probStop1) {probStop=probStop1;}

	double get_probStop () {return probStop;}

	void set_probStoph (double probStop1) {probStoph=probStop1;}

	double get_probStoph () {return probStoph;}

	void set_arrDelay (double arrDelay1) {arrDelay=arrDelay1;}

	double get_arrDelay () {return arrDelay;}

	void set_depDelay (double depDelay1) {depDelay=depDelay1;}

	double get_depDelay () {return depDelay;}

	void set_dwellDelay (double dwlDelay1) {dwlDelay=dwlDelay1;}

	double get_dwellDelay () {return dwlDelay;}

	void set_rideDelay (double rideDelay1) {rideDelay=rideDelay1;}

	double get_rideDelay () {return rideDelay;}

	void set_PVal (double PVal1) {PVal=PVal1;}

	double get_PVal () {return PVal;}

	void set_AVal (double AVal1) {AVal=AVal1;}

	double get_AVal () {return AVal;}

	void set_CRdTmE (double CRdTmE1) {CRdTmE=CRdTmE1;}

	double get_CRdTmE () {return CRdTmE;}

	void calCRdTmE (double maxCRdTm1) {CRdTmE=maxCRdTm1-CRdTm;}

	void set_hWkTmOns (double hWkTmOns1) {hWkTmOns=hWkTmOns1;}

	double get_hWkTmOns () {return hWkTmOns;}

	void set_hWkTmOffs (double hWkTmOffs1) {hWkTmOffs=hWkTmOffs1;}

	double get_hWkTmOffs () {return hWkTmOffs;}

	void set_WkTmOns (double WkTmOns1) {WkTmOns=WkTmOns1;}

	double get_WkTmOns () {return WkTmOns;}

	void set_WkTmOffs (double WkTmOffs1) {WkTmOffs=WkTmOffs1;}

	double get_WkTmOffs () {return WkTmOffs;}

	void set_WalkCost (double WalkCost1) {WalkCost=WalkCost1;}

	double get_WalkCost () {return WalkCost;}

	void set_RideCost (double RideCost1) {RideCost=RideCost1;}

	double get_RideCost () {return RideCost;}

	void set_OperCost (double OperCost1) {OperCost=OperCost1;}

	double get_OperCost () {return OperCost;}

	void set_TCost (double TCost1) {TCost=TCost1;}

	double get_TCost () {return TCost;}

	void set_xc (double xc1) {xc=xc1;}

	double get_xc () const {return xc;}

	void set_yc (double yc1) {yc=yc1;}

	double get_yc () const {return yc;}

	void set_zc (double zc1) {zc=zc1;}

	double get_zc () const {return zc;}

  friend std::ostream& operator<<(
    std::ostream& os, const tstop& ts) {
    return os << ts.id << "\t" << ts.tripid << "\t"<< ts.StopName<<"\t"<< ts.Edgeid <<"\t"
      << ts.CRdTm <<"\t" << ts.CRdTmC <<"\t" << ts.CRdTmE<<"\t"
	  << ts.HistOns<<"\t" << ts.HistOffs <<"\t" << ts.HistDepVol <<"\t"<< ts.Ons<<"\t"
	  << ts.Offs <<"\t"<< ts.DepVol <<"\t"<< ts.probStop <<"\t"<< ts.PVal<<"\t"
	  << ts.AVal <<"\t"<< ts.depDelay <<"\t" << ts.arrDelay<<"\t" << ts.rideDelay<<"\t" 
	  << ts.blnHiStop<<"\t"<< ts.blnInbd <<"\t"<< ts.blnIncl <<"\t"<< ts.blnExtr<<"\t"<< ts.blnElim<<"\t"<<ts.xc<<"\t"
	  <<ts.yc<<"\t"<<ts.zc<<"\t"<<endl ;
  }

    void ts_hdr(ostream& out)
    { 
     out << "id"<< "\t" << "TripId"<< "\t" << "StopName"<<"\t"<< "Edgeid" <<"\t"
      << "CRdTm (min)"<<"\t" << "CRdTmC (min)"<<"\t" << "CRdTmE (min)"<<"\t"
	  << "HistOns (Pax/Hr)"<<"\t" << "HistOffs (Pax/Hr)"<<"\t" << "HistDepVol (Pax/Hr)"<<"\t"<< "Ons (Pax/Hr)"<<"\t"
	  << "Offs (Pax/Hr)"<<"\t"<< "DepVol (Pax/Hr)"<<"\t"<< "probStop"<<"\t"<< "PVal"<<"\t"
	  << "AVal"<<"\t"<< "depDelay (secs)"<<"\t" << "arrDelay (secs)"<<"\t" << "rideDelay (Pers-hrs/hr)"<<"\t" 
	  << "blnHiStop"<<"\t"<< "blnInbd"<<"\t"<< "blnIncl"<<"\t"<< "blnExtr" <<"\t"<< "blnElim"<<"\t" <<"XC"<<"\t"<<"YC"<<"\t"<<"ZC"<<"\t"<<endl ;
  }
  
    void show_tsf(ostream& out)
    { 
		out << get_id() << "\t" << tripId() << "\t" << get_StopName() << "\t" << get_Edgeid() << "\t" 
		  << get_blnHist() << "\t" << get_CumDist() << "\t" << get_CRdTm()<< "\t"<< 
		  get_HistOns()<<"\t"<< get_HistOffs()<< "\t"<< get_HistDepVol() << "\t"<< 
		  get_Ons()<<"\t"<< get_Offs() << "\t"<<get_DepVol()<< "\t"<< 
		  get_PVal()<< "\t"<< get_AVal()<< get_rideDelay()<< "\t"<< endl;
	};

    void show_tsfhdr(ostream& out)
    { 
		out  << "Stopid" << "\t" << "TripId" << "\t" << "StopName" << "\t" << "Edgeid" << "\t" 
		 << "Hist" << "\t" << "CumDist" << "\t" << "CRdTm (min)"<<"\t"<<"HOns (Pax/hr)" << "\t"
		 <<"HOffs (Pax/hr)"<< "\t"<<"HDepVol (Pax/hr)" <<endl;
	};

	void show_tsdh(void)
    { 
		cout << "Stop Id: " << get_id() << endl;
		cout << "Trip Id: " << tripId() << endl;
		cout << "Stop Name: " << get_StopName() << endl;
		cout << "Edge Id: " << get_Edgeid() << endl;
		cout << "Hist: " << get_blnHist() << endl;
		cout << "Inbound: " << get_blnInbd() << endl;
		cout << "Cum. Dist: " << get_CumDist() << endl;
		cout << "Historic Cum. Time: " << get_CRdTm() << endl;
		cout << "Hist Ons: " << get_HistOns() << endl;
		cout << "Hist Offs: " << get_HistOffs() << endl;
		cout << "Hist Dep Vol: " << get_HistDepVol() << endl;
	};

	void show_tsd(void)
    { 
		cout << "Stop Id: " << get_id() << endl;
		cout << "Trip Id: " << tripId() << endl;
		cout << "Stop Name: " << get_StopName() << endl;
		cout << "Edge Id: " << get_Edgeid() << endl;
		cout << "Included: " << get_blnIncl() << endl;
		cout << "Inbound: " << get_blnInbd() << endl;
		cout << "Cum. Dist: " << get_CumDist() << endl;
		cout << "Computed Cum. Time: " << get_CRdTmC() << endl;
		cout << " Ons: " << get_Ons() << endl;
		cout << " Offs: " << get_Offs() << endl;
		cout << " Dep Vol: " << get_DepVol() << endl;
	};

		void tstop::serialize(ofstream& pts);
		void tstop::serializetext(ofstream& pts);
		void tstop::serializetexthdr(ofstream& pts);
		void tstop::deserialize(ifstream& pts);
		void tstop::deserializetext(ifstream& pts);

	~tstop(){
			//cout << "Deleting stop: " << *this << endl;
	}

protected:
			long id; // Stop id j
			long tripid; // trip Id incase of trip level run 
			long Stopidp; // predecessor Stop id i 
			long Stopids; // successor Stop id k
			long Edgeid;
			double palong;
			double nearDist;  // nearest distnace to the route
			int StOrdr; // Stop order
			string StopLbl; // Stop label
			string StopName; // Stop Name
	        short lbl; // labeled
			short blnHiStop; // true <>0 if historic
			short blnInbd; // true <>0 if inbound
			short blnExtr; // true <>0 if external
			short blnIncl; // true <>0 if included
			short blnElim; // true <>0 if eliminated
			double CumDist; // supplied Cumulative ride dist
			double CRdTm; // supplied Cumulative ride time
			double undCRdTm; // undelayed Cumulative ride time
			double CRdTmC;  // calculated Cumulative ride time
			double HistOns; 
			double HistOffs; 
			double HistDepVol; 
			double Ons; 
			double Offs; 
			double DepVol; 
			double probStoph; 
			double probStop; 
			double depDelay; 
			double arrDelay; 
			double dwlDelay; 
			double rideDelay; 
			double PVal;
			double AVal;
			double CRdTmE; 
			double hWkTmOns; 
			double hWkTmOffs;
			double WkTmOns; 
			double WkTmOffs;
			double WalkCost;
			double RideCost;
			double OperCost;
			double TCost;
         	double xc; // X Coordinate
         	double yc; // Y Coordinate
         	double zc; // Y Coordinate
};

void tstop::serialize(ofstream& pts)
{
 pts.write(reinterpret_cast<char *>(&id), sizeof(id));
 pts.write(reinterpret_cast<char *>(&tripid), sizeof(tripid));
 streamsize sizet=StopLbl.size();// store stop label's length
 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pts.write(StopLbl.c_str(), sizet+1); // write string label and final '\0'
 pts.write(reinterpret_cast<char *>(&Stopidp), sizeof(Stopidp));
 pts.write(reinterpret_cast<char *>(&Stopids), sizeof(Stopids));
 pts.write(reinterpret_cast<char *>(&Edgeid), sizeof(Edgeid));
 pts.write(reinterpret_cast<char *>(&palong), sizeof(palong));
 pts.write(reinterpret_cast<char *>(&nearDist), sizeof(nearDist));
 pts.write(reinterpret_cast<char *>(&StOrdr), sizeof(StOrdr));
 sizet=StopName.size();// store StopName's length
 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pts.write(StopName.c_str(), sizet+1); // write stopname including '\0' too
 pts.write(reinterpret_cast<char *>(&lbl),sizeof(lbl));
 pts.write(reinterpret_cast<char *>(&blnHiStop),sizeof(blnHiStop));
 pts.write(reinterpret_cast<char *>(&blnInbd),sizeof(blnInbd));
 pts.write(reinterpret_cast<char *>(&blnExtr),sizeof(blnExtr));
 pts.write(reinterpret_cast<char *>(&blnIncl),sizeof(blnIncl));
 pts.write(reinterpret_cast<char *>(&blnElim),sizeof(blnElim));
 pts.write(reinterpret_cast<char *>(&CumDist),sizeof(CumDist));
 pts.write(reinterpret_cast<char *>(&CRdTm),sizeof(CRdTm));
 pts.write(reinterpret_cast<char *>(&undCRdTm),sizeof(undCRdTm));
 pts.write(reinterpret_cast<char *>(&CRdTmC),sizeof(CRdTmC));
 pts.write(reinterpret_cast<char *>(&HistOns),sizeof(HistOns));
 pts.write(reinterpret_cast<char *>(&HistOffs),sizeof(HistOffs));
 pts.write(reinterpret_cast<char *>(&HistDepVol),sizeof(HistDepVol));
 pts.write(reinterpret_cast<char *>(&Ons),sizeof(Ons));
 pts.write(reinterpret_cast<char *>(&Offs),sizeof(Offs));
 pts.write(reinterpret_cast<char *>(&DepVol),sizeof(DepVol));
 pts.write(reinterpret_cast<char *>(&probStoph),sizeof(probStoph));
 pts.write(reinterpret_cast<char *>(&probStop),sizeof(probStop));
 pts.write(reinterpret_cast<char *>(&depDelay),sizeof(depDelay));
 pts.write(reinterpret_cast<char *>(&arrDelay),sizeof(arrDelay));
 pts.write(reinterpret_cast<char *>(&dwlDelay),sizeof(dwlDelay));
 pts.write(reinterpret_cast<char *>(&rideDelay),sizeof(rideDelay));
 pts.write(reinterpret_cast<char *>(&PVal),sizeof(PVal));
 pts.write(reinterpret_cast<char *>(&AVal),sizeof(AVal));
 pts.write(reinterpret_cast<char *>(&CRdTmE),sizeof(CRdTmE));
 pts.write(reinterpret_cast<char *>(&hWkTmOns),sizeof(hWkTmOns));
 pts.write(reinterpret_cast<char *>(&hWkTmOffs),sizeof(hWkTmOffs));
 pts.write(reinterpret_cast<char *>(&WkTmOns),sizeof(WkTmOns));
 pts.write(reinterpret_cast<char *>(&WkTmOffs),sizeof(WkTmOffs));
 pts.write(reinterpret_cast<char *>(&WalkCost),sizeof(WalkCost));
 pts.write(reinterpret_cast<char *>(&RideCost),sizeof(RideCost));
 pts.write(reinterpret_cast<char *>(&OperCost),sizeof(OperCost));
 pts.write(reinterpret_cast<char *>(&TCost),sizeof(TCost));
 pts.write(reinterpret_cast<char *>(&xc),sizeof(xc));
 pts.write(reinterpret_cast<char *>(&yc),sizeof(yc));
 pts.write(reinterpret_cast<char *>(&zc),sizeof(zc));

}


void tstop::deserialize(ifstream& pts)
{
 int len=0;
 char *p=0;
 pts.read(reinterpret_cast<char *>(&id), sizeof(id));
 pts.read(reinterpret_cast<char *>(&tripid), sizeof(tripid));
 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for stop Label
 pts.read(p, len+1); // copy string to temp, including '\0'
 StopLbl=p; // copy temp to data member
 delete[] p;
 pts.read(reinterpret_cast<char *>(&Stopidp), sizeof(Stopidp));
 pts.read(reinterpret_cast<char *>(&Stopids), sizeof(Stopids));
 pts.read(reinterpret_cast<char *>(&Edgeid), sizeof(Edgeid));
 pts.read(reinterpret_cast<char *>(&palong), sizeof(palong));
 pts.read(reinterpret_cast<char *>(&nearDist), sizeof(nearDist));
 pts.read(reinterpret_cast<char *>(&StOrdr), sizeof(StOrdr));


 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 pts.read(p, len+1); // copy string to temp, including '\0'
 StopName=p; // copy temp to data member
 delete[] p;

 pts.read(reinterpret_cast<char *>(&lbl),sizeof(lbl));
 pts.read(reinterpret_cast<char *>(&blnHiStop),sizeof(blnHiStop));
 pts.read(reinterpret_cast<char *>(&blnInbd),sizeof(blnInbd));
 pts.read(reinterpret_cast<char *>(&blnExtr),sizeof(blnExtr));
 pts.read(reinterpret_cast<char *>(&blnIncl),sizeof(blnIncl));
 pts.read(reinterpret_cast<char *>(&blnElim),sizeof(blnElim));
 pts.read(reinterpret_cast<char *>(&CumDist),sizeof(CumDist));
 pts.read(reinterpret_cast<char *>(&CRdTm),sizeof(CRdTm));
 pts.read(reinterpret_cast<char *>(&undCRdTm),sizeof(undCRdTm));
 pts.read(reinterpret_cast<char *>(&CRdTmC),sizeof(CRdTmC));
 pts.read(reinterpret_cast<char *>(&HistOns),sizeof(HistOns));
 pts.read(reinterpret_cast<char *>(&HistOffs),sizeof(HistOffs));
 pts.read(reinterpret_cast<char *>(&HistDepVol),sizeof(HistDepVol));
 pts.read(reinterpret_cast<char *>(&Ons),sizeof(Ons));
 pts.read(reinterpret_cast<char *>(&Offs),sizeof(Offs));
 pts.read(reinterpret_cast<char *>(&DepVol),sizeof(DepVol));
 pts.read(reinterpret_cast<char *>(&probStoph),sizeof(probStoph));
 pts.read(reinterpret_cast<char *>(&probStop),sizeof(probStop));
 pts.read(reinterpret_cast<char *>(&depDelay),sizeof(depDelay));
 pts.read(reinterpret_cast<char *>(&arrDelay),sizeof(arrDelay));
 pts.read(reinterpret_cast<char *>(&dwlDelay),sizeof(dwlDelay));
 pts.read(reinterpret_cast<char *>(&rideDelay),sizeof(rideDelay));
 pts.read(reinterpret_cast<char *>(&PVal),sizeof(PVal));
 pts.read(reinterpret_cast<char *>(&AVal),sizeof(AVal));
 pts.read(reinterpret_cast<char *>(&CRdTmE),sizeof(CRdTmE));
 pts.read(reinterpret_cast<char *>(&hWkTmOns),sizeof(hWkTmOns));
 pts.read(reinterpret_cast<char *>(&hWkTmOffs),sizeof(hWkTmOffs));
 pts.read(reinterpret_cast<char *>(&WkTmOns),sizeof(WkTmOns));
 pts.read(reinterpret_cast<char *>(&WkTmOffs),sizeof(WkTmOffs));
 pts.read(reinterpret_cast<char *>(&WalkCost),sizeof(WalkCost));
 pts.read(reinterpret_cast<char *>(&RideCost),sizeof(RideCost));
 pts.read(reinterpret_cast<char *>(&OperCost),sizeof(OperCost));
 pts.read(reinterpret_cast<char *>(&TCost),sizeof(TCost));
 pts.read(reinterpret_cast<char *>(&xc),sizeof(xc));
 pts.read(reinterpret_cast<char *>(&yc),sizeof(yc));
 pts.read(reinterpret_cast<char *>(&zc),sizeof(zc));
}


tstop::tstop(long id1,long tripId1=-1,string StopLbl1="", long Stopidp1=-1,long Stopids1=-1, int StOrdr1=-1,  string StopName1="", 
			 double CumDist1=0, double CRdTm1=0, double undCRdTm1=0, double CRdTmC1=0, short blnHist1=1, 
			 short blnInbd1=1,short blnIncl1=1, short blnExtr1=0, short blnElim1=0,double HistOns1=-1, 
			 double HistOffs1=-1, double HistDepVol1=0, double depDelay1=3.5, double arrDelay1=5, 
			 double probStoph1=0,double probStop1=0, double Ons1=-1, double Offs1=-1, double DepVol1=0,  
			 double dwlDelay1=0, double rideDelay1=0, double PVal1=0,double AVal1=0,double CRdTmE1=0, 
			 double hWkTmOns1=0, double hWkTmOffs1=0,double WkTmOns1=0, double WkTmOffs1=0, 
			 double WalkCost1=0,double RideCost1=0, double OperCost1=0,double TCost1=0,
			 double xc1=-1,double yc1=-1,double zc1=-1,double nearDist1=-1) 
{
	tstop::id = id1;
	tstop::tripid = tripId1;
	tstop::Stopidp = Stopidp1;
	tstop::Stopids = Stopids1;
	tstop::StopLbl= StopLbl1;
	tstop::StopName = StopName1;
	tstop::StOrdr = StOrdr1;
	tstop::blnHiStop=blnHist1; 
	tstop::blnInbd=blnInbd1; 
	tstop::blnExtr=blnExtr1; 
	tstop::blnIncl=blnIncl1;
	tstop::blnElim=blnElim1;
	tstop::CumDist=CumDist1;
	tstop::CRdTm=CRdTm1; 
	tstop::undCRdTm=undCRdTm1; 
	tstop::CRdTmC=CRdTmC1; 
	tstop::HistOns=HistOns1; 
	tstop::HistOffs=HistOffs1; 
	tstop::HistDepVol=HistDepVol1; 
	tstop::probStoph=probStoph1; 
	tstop::Ons= Ons1;
	tstop::Offs = Offs1;
	tstop::DepVol = DepVol1;
	tstop::probStop=probStop1; 
	tstop::depDelay=depDelay1; 
	tstop::arrDelay=arrDelay1; 
	tstop::dwlDelay=dwlDelay1; 
	tstop::rideDelay=rideDelay1;
	tstop::PVal=PVal1; 
	tstop::AVal=AVal1; 
	tstop::CRdTmE=CRdTmE1; 
	tstop::hWkTmOns=hWkTmOns1;
	tstop::hWkTmOffs=hWkTmOffs1; 
	tstop::WkTmOns=WkTmOns1;
	tstop::WkTmOffs=WkTmOffs1; 
	tstop::WalkCost=WalkCost1; 
	tstop::RideCost=RideCost1; 
	tstop::OperCost=OperCost1; 
	tstop::TCost=TCost1; 
}

void tstop::tsr(double undCRdTm1, double CRdTmC1, double Ons1, double Offs1, double DepVol1, double probStop1, double depDelay1, 
		double arrDelay1, double dwlDelay1,double rideDelay1, double PVal1,double AVal1,double CRdTmE1, double WkTmOns1, double WkTmOffs1,
		double WalkCost1,double RideCost1,double OperCost1,double TCost1) 
{
	tstop::undCRdTm = undCRdTm1;
	tstop::CRdTmC = CRdTmC1;
	tstop::Ons= Ons1;
	tstop::Offs = Offs1;
	tstop::DepVol = DepVol1;
	tstop::probStop=probStop1; 
	tstop::depDelay=depDelay1; 
	tstop::arrDelay=arrDelay1; 
	tstop::dwlDelay=dwlDelay1; 
	tstop::rideDelay=rideDelay1;
	tstop::PVal=PVal1; 
	tstop::AVal=AVal1; 
	tstop::CRdTmE=CRdTmE1; 
	tstop::WkTmOns=WkTmOns1;
	tstop::WkTmOffs=WkTmOffs1; 
	tstop::WalkCost=WalkCost1; 
	tstop::RideCost=RideCost1; 
	tstop::OperCost=OperCost1; 
	tstop::TCost=TCost1; 
}

void tstop::cprobstop(double hdwy1) 
{
	probStop = 1 - exp(-(hdwy1/60*(Ons+Offs))); // hdwy in minutes ons+offs in pax/hr
}

void tstop::cprobstoph(double hdwy1) 
{
	probStoph = 1 - exp(-(hdwy1/60*(HistOns+HistOffs))); // hdwy in mins, ons+offs in pax/hr
}

void tstop::cDepVol(double depVol1) 
{
	DepVol = depVol1 - Offs + Ons;
}

void tstop::cDepVolh(double depVol1) 
{
	HistDepVol = depVol1 - HistOffs + HistOns;
}
void tstop::cdwlDelay(double unitOnTime, double unitOffTime, double hdwy ) 
{
	if (unitOnTime > 0 && unitOnTime > 0) {
		dwlDelay = (unitOnTime * Ons  + unitOffTime * Offs ) * hdwy/60 ;
	}
}

double tstop::cstopDelay() 
{
		return ((arrDelay+depDelay)*probStop/60) ;
}

void tstop::cdwlDelayh(double unitOnTime, double unitOffTime, double hdwy ) 
{
	if (unitOnTime > 0 && unitOnTime > 0) {
		dwlDelay = (unitOnTime * HistOns  + unitOffTime * HistOffs) * hdwy/60 ;
	}
}

double tstop::thruVol() 
{
	return (DepVol  - Offs);
}

double tstop::thruVolh() 
{
	return (HistDepVol  - HistOffs);
}

double tstop::segDelay(double probjm1, double depDjm1, double stoProb) 
{
  return (probjm1 * depDjm1 + dwlDelay + stoProb * arrDelay );
}


void tstop::serializetext(ofstream& pts)
{
 pts <<id<<"\t"<<tripid<<"\t"<<Stopidp<<"\t"<<Stopids<<"\t"<<Edgeid<<"\t"<<palong<<"\t"<<StOrdr<<"\t"<<StopLbl<<"\t" 
	 <<StopName<<"\t"<<xc<<"\t"<<yc<<"\t"<<zc<<"\t"<<lbl<<"\t"<<blnHiStop<<"\t"<<blnInbd<<"\t"<<blnExtr<<"\t"<<blnIncl<<"\t"<< blnElim<<"\t"
	 <<CumDist<<"\t"<<CRdTm<<"\t"<<undCRdTm<<"\t"<<CRdTmC<<"\t"<<HistOns<<"\t"<< HistOffs<<"\t"
	 <<HistDepVol<<"\t"<<Ons<<"\t"<<Offs<<"\t"<<DepVol<<"\t"<<probStoph<<"\t"<< probStop<<"\t"<<depDelay<< "\t"
	<<arrDelay<<"\t"<<dwlDelay<<"\t"<<rideDelay<<"\t"<< PVal<<"\t"<<AVal<<"\t"<<CRdTmE<<"\t"<<hWkTmOns<<"\t"
	<<hWkTmOffs<<"\t"<<WkTmOns<<"\t"<< WkTmOffs<<"\t"<<WalkCost<<"\t"<<RideCost<<"\t"<<OperCost<<"\t"<<TCost<<endl;
}
void tstop::serializetexthdr(ofstream& pts)
{
 pts <<"id"<<"\t"<<"TripId"<<"\t"<<"Stopidp"<<"\t"<<"Stopids"<<"\t"<<"Edgeid"<<"\t"<<"palong"<<"\t"<<"StOrdr"<<"\t"<<"StopLbl"<<"\t" 
	 <<"StopName"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"z"<<"\t"<<"lbl"<<"\t"<<"blnHist"<<"\t"<<"blnInbd"<<"\t"<<"blnExtr"<<"\t"<<"blnIncl"<<"\t"<<"blnElim"<<"\t"
	 <<"CumDist"<<"\t"<<"CRdTm (min)"<<"\t"<<"undCRdTm (min)"<<"\t"<<"CRdTmC (min)"<<"\t"<<"H.Ons (Pax/Hr)"<<"\t"<< "H.Offs (Pax/Hr)"<<"\t"
	 <<"H.DepVol (Pax/Hr)"<<"\t"<<"Ons (Pax/Hr)"<<"\t"<<"Offs (Pax/Hr)"<<"\t"<<"DepVol (Pax/Hr)"<<"\t"<<"probStoph"<<"\t"<< "probStop"<<"\t"<<"depDelay"<< "\t"
	<<"arrDelay (Secs)"<<"\t"<<"dwlDelay (Hrs)"<<"\t"<<"RideDelay (Pers-Hrs/Hr)"<<"\t"<<"PVal"<<"\t"<<"AVal"<<"\t"<<"CRdTmE (Mins)"<<"\t"<<"hWkTmOns (Pax-Min/Hr)"<<"\t"
	<<"hWkTmOffs (Pax-Min/Hr)"<<"\t"<<"WkTmOns (Pax-Min/Hr)"<<"\t"<< "WkTmOffs (Pax-Min/Hr)"<<"\t"<<"WalkCost ($/Pd)"<<"\t"<<"RideCost ($/Pd)"<<"\t"
	<<"OperCost ($/Pd)"<<"\t"<<"TCost ($/Pd)"<<endl;
}


void tstop::deserializetext(ifstream& pts)
{
 pts >>id>>tripid>>Stopidp>>Stopids>>Edgeid>>palong>>StOrdr>>StopLbl 
	 >>StopName>>xc>>yc>>zc>>lbl>>blnHiStop>>blnInbd>>blnExtr>>blnIncl>> blnElim
	 >>CumDist>>CRdTm>>undCRdTm>>CRdTmC>>HistOns>> HistOffs
	 >>HistDepVol>>Ons>>Offs>>DepVol>>probStoph>> probStop>>depDelay
	>>arrDelay>>dwlDelay>>rideDelay>> PVal>>AVal>>CRdTmE>>hWkTmOns
	>>hWkTmOffs>>WkTmOns>> WkTmOffs>>WalkCost>>RideCost>>OperCost>>TCost;
}


// dp class derived from tstop  

class dptStop 
{
private:
	tstop tsem;
	long i,j,k,l,m,q;
	string dpkey;
public:
	dptStop() : i(-1),j(-1),k(-1),l(-1),m(-1) {};
	dptStop(const dptStop& dpts )  {
		dptStop::i = dpts.i;
		dptStop::j = dpts.j;
		dptStop::k = dpts.k;
		dptStop::l = dpts.l;
		dptStop::m = dpts.m;
		dptStop::q = dpts.q;
		dptStop::dpkey = dpts.dpkey;
		dptStop::tsem = dpts.tsem;
	}
/*  dptStop(tstop tst,long i1=-1, long j1=-1, long l1=-1,long m1=-1,long q1 = -1,
	  string dpkey1="-1")
    : tsem(tst),dpkey(makey(i1,tst.get_Stopidp(),tst.get_id(),tst.get_Stopids(),m1)),i(i1), 
	j(tst.get_Stopidp()), k(tst.get_id()), l(tst.get_Stopids()), m(m1),q(q1)  {
      //cout << "Creating DPStop Object from a stop object: " << *this << endl;
  } */
  dptStop(tstop tst,long i1=-1, long j1=-1, long k1=-1,long l1=-1,long m1=-1,long q1 = -1)
    : tsem(tst),dpkey(makey(i1,j1,tst.get_id(),l1,m1)),i(i1), 
	j(j1), k(tst.get_id()), l(l1), m(m1),q(q1)  {
      //cout << "Creating DPStop Object from a stop object: " << *this << endl;
  }


 /* dptStop(long i1=-1, long j1=-1, long k1=-1,long l1=-1,long m1=-1)
    : i(i1), j(j1), k(k1), l(l1), m(m1) {
		//cout << "Creating DPStop Object with out a stop object: " << *this << endl;
  }
  */
  tstop get_tstop() const { return tsem; }

  void set_tstop(tstop& tsp) { tsem=tsp; }

 dptStop& operator=(dptStop& dpts) {
    dptStop::i = dpts.get_i();
    dptStop::j = dpts.get_j();
    dptStop::k = dpts.get_k();
    dptStop::l = dpts.get_l();
    dptStop::m = dpts.get_m();
    dptStop::q = dpts.get_q();
    dptStop::dpkey = dpts.get_dpkey();
	dptStop::tsem = dpts.get_tstop();

    return *this;
  }
  friend ostream&
  operator<<(ostream& os, const dptStop& ts) {
    return os <<ts.q<<"\t"<< ts.i << "\t" << ts.j << "\t" << ts.k<< "\t" << ts.l
	<<"\t"<<ts.m<<"\t"<< ts.dpkey<< "\t"<< ts.tsem<<endl;
  }


 void dptStop::serializetext(ofstream& pts)
{
 pts <<q<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<l<<"\t"<<m<<"\t"<<dpkey<<"\t";
 tsem.serializetext(pts);
}

void dptStop::serializetext3d(ofstream& pts)
{
 pts <<q<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<dpkey<<"\t";
 tsem.serializetext(pts);
}

void dptStop::serializetext2(ofstream& pts)
{
 pts <<q<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<l<<"\t"<<m<<"\t"<<dpkey<<"\t";
}

void dptStop::serializetext2d3(ofstream& pts)
{
 pts <<q<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<dpkey<<"\t";
}

void dptStop::serializetexthdr(ofstream& pts)
{
 pts<<"ScId"<<"\t"<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"l"<<"\t"<<"m"<<"\t"<<"dpkey"<<"\t";
 tsem.serializetexthdr(pts);
}

void dptStop::serializetext3dhdr(ofstream& pts)
{
 pts<<"ScId"<<"\t"<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"dpkey"<<"\t";
 tsem.serializetexthdr(pts);
}

void dptStop::deserializetext(ifstream& pts)
{
 pts >>q>>i>>j>>k>>l>>m>>dpkey;
 tsem.deserializetext(pts);
}

void dptStop::deserializetext3d(ifstream& pts)
{
 pts >>q>>i>>j>>k>>dpkey;
 tsem.deserializetext(pts);
}

void dptStop::serialize(ofstream& pts)
{
 pts.write(reinterpret_cast<char *>(&q), sizeof(q));
 pts.write(reinterpret_cast<char *>(&i), sizeof(i));
 pts.write(reinterpret_cast<char *>(&j), sizeof(j));
 pts.write(reinterpret_cast<char *>(&k), sizeof(k));
 pts.write(reinterpret_cast<char *>(&l), sizeof(l));
 pts.write(reinterpret_cast<char *>(&m), sizeof(m));

 streamsize sizet=dpkey.size();// store dp key's length
 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pts.write(dpkey.c_str(), sizet+1); // write string dpkey and final '\0'
 tsem.serialize(pts);
}


void dptStop::deserialize(ifstream& pts)
{
 int len=0;
 char *p=0;

 pts.read(reinterpret_cast<char *>(&q), sizeof(q));
 pts.read(reinterpret_cast<char *>(&i), sizeof(i));
 pts.read(reinterpret_cast<char *>(&j), sizeof(j));
 pts.read(reinterpret_cast<char *>(&k), sizeof(k));
 pts.read(reinterpret_cast<char *>(&l), sizeof(l));
 pts.read(reinterpret_cast<char *>(&m), sizeof(m));

 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for stop Label
 pts.read(p, len+1); // copy string to temp, including '\0'
 dpkey=p; // copy temp to data member
 delete[] p;
 tsem.deserialize(pts);
}


virtual ~dptStop() { 
//	  delete tsem;
    //cout << "Deleting dpstop: " << *this << endl;
  }

	void set_i (long i1) {i=i1;}
	long get_i () {return i;}
	void set_j (long j1) {j=j1;}
	long get_j () {return j;}
	void set_k (long k1) {k=k1;}
	long get_k () {return k;}
	void set_l (long l1) {l=l1;}
	long get_l () {return l;}
	void set_m (long m1) {m=m1;}
	long get_m () {return m;}
	void set_q (long q1) {q=q1;}
	long get_q () {return q;}
	void set_dpkey(string dpkey1) {dpkey=dpkey1;}
	string get_dpkey () {return dpkey;}

	string makey(long i1,long j1,long k1,long l1,long m1,char *sp="|")
	{
		return (to_string(i1) + *sp + to_string(j1) + *sp + to_string(k1) + *sp + to_string(l1) + *sp +  to_string(m1));
	}
	string makey(dptStop& dpts,char *sp="|")
	{
		return (to_string(dpts.get_i()) + *sp + to_string(dpts.get_j()) + *sp + to_string(dpts.get_k()) + *sp + to_string(dpts.get_l()) + *sp +  to_string(dpts.get_m()));
	}
	string makey3d(dptStop& dpts,char *sp="|")
	{
		return (to_string(dpts.get_i()) + *sp + to_string(dpts.get_j()) + *sp + to_string(dpts.get_k()));
	}
	string makey5d(dptStop& dpts,char *sp="|")
	{
		return (to_string(dpts.get_i()) + *sp + to_string(dpts.get_j()) + *sp + to_string(dpts.get_k()) + *sp + to_string(dpts.get_l()) + *sp +  to_string(dpts.get_m()));
	}
	string makey7d(dptStop& dpts,char *sp="|")
	{
		return (to_string(dpts.get_i()) + *sp + to_string(dpts.get_j()) + *sp + to_string(dpts.get_k()) + *sp + to_string(dpts.get_l()) + *sp +  to_string(dpts.get_m()));
	}

};



// TripDP class derived from dptStop  

class tripDPs 
{
private:
	dptStop dptSChld;
	long tripId;
	double walk;
	double ride;
	string tDPkey;
public:
	tripDPs() : tripId(-1), walk(10),ride(5) {};
	tripDPs(const tripDPs& tdpts )  {
		tripDPs::tripId = tdpts.tripId;
		tripDPs::walk = tdpts.walk;
		tripDPs::ride = tdpts.ride;
		tripDPs::tDPkey = tdpts.tDPkey;
		tripDPs::dptSChld = tdpts.get_dptStop();
	}
/*  dptStop(tstop tst,long i1=-1, long j1=-1, long l1=-1,long m1=-1,long q1 = -1,
	  string dpkey1="-1")
    : dptSChld(tst),dpkey(makey(i1,tst.get_Stopidp(),tst.get_id(),tst.get_Stopids(),m1)),i(i1), 
	j(tst.get_Stopidp()), k(tst.get_id()), l(tst.get_Stopids()), m(m1),q(q1)  {
      //cout << "Creating DPStop Object from a stop object: " << *this << endl;
  } */
  tripDPs(dptStop trpDP,long tripid1=-1, double walk = -1, double ride=-1)
    : dptSChld(trpDP),tDPkey(makey(tripId))  {
      //cout << "Creating DPStop Object from anothe DPStop object: " << *this << endl;
  }


 /* dptStop(long i1=-1, long j1=-1, long k1=-1,long l1=-1,long m1=-1)
    : i(i1), j(j1), k(k1), l(l1), m(m1) {
		//cout << "Creating DPStop Object with out a stop object: " << *this << endl;
  }
  */
  dptStop get_dptStop() const { return dptSChld; }

  void set_dptStop(dptStop& ts) { dptSChld=ts; }

 tripDPs& operator=( tripDPs& trdp) {
    tripDPs::tripId = trdp.get_tripId();
    tripDPs::walk = trdp.get_walk();
    tripDPs::ride = trdp.get_ride();
    tripDPs::tDPkey = trdp.get_tDPkey();
	tripDPs::dptSChld = trdp.get_dptStop();
    return *this;
  }
  friend ostream&
  operator<<(ostream& os, const tripDPs& trdp) {
    return os <<trdp.tripId<<"\t"<<trdp.walk<<"\t"<<trdp.ride<<"\t"<< trdp.tDPkey<< "\t"<< trdp.dptSChld<<endl;
  }


 void tripDPs::serializetext(ofstream& pts)
{
 pts <<tripId<<"\t"<<walk<<"\t"<<ride<<"\t"<<tDPkey<<"\t";
 dptSChld.serializetext(pts);
}

void tripDPs::serializetext3d(ofstream& pts)
{
 pts <<tripId<<"\t"<<walk<<"\t"<<ride<<"\t"<<tDPkey<<"\t";
 dptSChld.serializetext(pts);
}

void tripDPs::serializetext2(ofstream& pts)
{
 pts <<tripId<<"\t"<<walk<<"\t"<<ride<<"\t"<<tDPkey<<"\t";
}

void tripDPs::serializetext2d3(ofstream& pts)
{
 pts <<tripId<<"\t"<<walk<<"\t"<<ride<<tDPkey<<"\t";
}

void tripDPs::serializetexthdr(ofstream& pts)
{
 pts<<"TripId"<<"\t"<<"WalkPara"<<"\t"<<"RidePara"<<"\t"<<"dpkey"<<"\t";
 dptSChld.serializetexthdr(pts);
}

void tripDPs::serializetext3dhdr(ofstream& pts)
{
 pts<<"TripId"<<"\t"<<"WalkPara"<<"\t"<<"RidePara"<<"\t"<<"dpkey"<<"\t";
 dptSChld.serializetexthdr(pts);
}

void tripDPs::deserializetext(ifstream& pts)
{
 pts >>tripId>>walk>>ride>>tDPkey;
 dptSChld.deserializetext(pts);
}

void tripDPs::deserializetext3d(ifstream& pts)
{
 pts>>tripId>>walk>>ride>>tDPkey;
 dptSChld.deserializetext(pts);
}

void tripDPs::serialize(ofstream& pts)
{
 pts.write(reinterpret_cast<char *>(&tripId), sizeof(tripId));
 pts.write(reinterpret_cast<char *>(&walk), sizeof(walk));
 pts.write(reinterpret_cast<char *>(&ride), sizeof(ride));

 streamsize sizet=tDPkey.size();// store dp key's length
 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pts.write(tDPkey.c_str(), sizet+1); // write string dpkey and final '\0'
 dptSChld.serialize(pts);
}


void tripDPs::deserialize(ifstream& pts)
{
 int len=0;
 char *p=0;

 pts.read(reinterpret_cast<char *>(&tripId), sizeof(tripId));
 pts.read(reinterpret_cast<char *>(&walk), sizeof(walk));
 pts.read(reinterpret_cast<char *>(&ride), sizeof(ride));

 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for stop Label
 pts.read(p, len+1); // copy string to temp, including '\0'
 tDPkey=p; // copy temp to data member
 delete[] p;
 dptSChld.deserialize(pts);
}


virtual ~tripDPs() { 
//	  delete dptSChld;
    //cout << "Deleting dpstop: " << *this << endl;
  }

	void set_tripId (long tripId1) {tripId=tripId1;}
	long get_tripId () {return tripId;}
	void set_walk (double walk1) {walk=walk1;}
	double get_walk () {return walk;}
	void set_ride (long ride1) {ride=ride1;}
	double get_ride () {return ride;}
	void set_tDPkey(string tDPkey1) {tDPkey=tDPkey1;}
	string get_tDPkey () {return tDPkey;}

	string makey(long tripid1, string sp="|")
	{
		return (to_string(tripid1) + sp + dptSChld.makey5d(dptSChld));
	}
	string makey3d(tripDPs& dpts,char *sp="|")
	{
		return (to_string(dpts.get_tripId()) + *sp + dptSChld.makey3d(dptSChld) );
	}
	string makey5d(tripDPs& dpts,char *sp="|")
	{
		return (to_string(dpts.get_tripId()) + *sp + dptSChld.makey5d(dptSChld) );
	}
	string makey7d(tripDPs& dpts,char *sp="|")
	{
		return (to_string(dpts.get_tripId()) + *sp + dptSChld.makey7d(dptSChld) );
	}

};




class parcel {
public:

    friend istream& operator>>( istream&, parcel& );
	
	enum EnumLbl { fixed=-1, free=0, labeled = 1 };

	void parcelUp(long pid1,long eoid1,long evp1,bool blnHist,short lbl1, long frid1, long toid1, 
		long origon1,long origoff1, double palong1, short on1,double lndsf1, double grarea1, 
		double lvarea1, double valtot1, double valand1,double valbld1,
        double cfact1,double pval1,double aval1,string ptype1) ; // constructor

		void parcel::serialize(ofstream& parc);
		void parcel::deserialize(ifstream& parc);
		void parcel::texthdr(ofstream& par);
		void parcel::serializetext(ofstream& par);
		void parcel::serializetexthdr(ofstream& par);

		parcel():  pid(0), pacid("-1"),eoid(0), evp(0),blnHist(1),  lbl(0), 
		hwkoncost(0),hwkoffcost(0), hwkxoncost(0), hwkxoffcost(0),
		wkoncost(0),wkoffcost(0),  wkxoncost(0),  wkxoffcost(0),  frid(0),  toid(0),ons(0.0),offs(0.0),
		origon(-1),origoff(-1),palong(0),on(-1), lndsf(0), grarea(0), lvarea(0), valtot(0),valand(0),
		 valbld(0),  cfact(1), pval(0), aval(0),ptype("-1"),lucd("-X"),xc(0),yc(0),zc(0),pnote("") {};

    void parc(long pid1=0,long eoid1=0,long evp1=0,bool blnHist1=1, short lbl1=0, double cost1=-1.0, double ons1=0.0,
		double offs1=0.0,double wkoncost1=0.0, double wkoffcost1=0.0, double wkxoncost1=0.0, double wkxoffcost1=0.0, 
		long frid1=-1, long toid1=-1, long origon1=-1, long origoff1=-1, double palong1=0, short on1=-1, double lndsf1=0, 
		double grarea1=0, double lvarea1=0, double valtot1=0, double valand1=0, double valbld1=0, double cfact1=1,
		double pval1=0,double aval1=0,string ptype1="-1", string lucd1="-X", double xc1=0, double yc1=0,double zc1=0) 
	{
   	pid = pid1; // parcel id
	eoid = eoid1; // edge id
	ons = ons1; // origin id
	offs = offs1; // origin id
	origon = origon1; // origin id
	origoff = origoff1; // origin id
	evp = evp1;  // parcel pointer 
    blnHist = blnHist1;
	lbl = lbl1;  // parcel label
	hwkoncost= wkoncost1;  // parcel cost
	hwkoffcost= wkoffcost1;  // parcel cost
	hwkxoncost= wkoncost1;  // parcel cost
	hwkxoffcost= wkoffcost1;  // parcel cost
	wkxoncost= wkoncost1;  // parcel cost
	wkxoffcost= wkoffcost1;  // parcel cost
	frid = frid1;  // parcel from vertex
	toid = toid1;  // parcel to vertex
	on = on1;  // going from or to parcel
	lndsf = lndsf1; // land sq ft
    grarea = grarea1;  // gross floor area
	lvarea = lvarea1; // living area
	valtot = valtot1; // total valuation 
	valand = valand1; // land valuation 
	valbld = valbld1; // Bldg valuation 
	cfact = cfact1; // competition factor 
	pval = pval1;  // parcel production strength factor
	aval = aval1; // parcel attraction strength factor
	xc = xc1; // x-coordinate
	yc = yc1; // y-coordinate
	zc = zc1; // y-coordinate
	};

	void set_id (long pid1) {pid=pid1;}

	long get_id () const {return pid;}

	void set_origon (long origon1) {origon=origon1;}

	long get_origon () const {return origon;}

	void set_origoff (long origoff1) {origoff=origoff1;}

	long get_origoff () const {return origoff;}

	void set_eoid (long eoid1) {eoid=eoid1;}

	long get_eoid () {return eoid;}

	void set_evp (long evp1) {evp=evp1;}

	long get_evp () {return evp;}

//	void set_par1 (parcel* par) {par1=par;}

	void set_hist (bool hist1) {blnHist=hist1;}

	bool get_hist () {return blnHist;}

	void set_lbl (short lbl1) {lbl=lbl1;}

	short get_lbl () {return lbl;}

	void set_hwkoncost (double hwkoncost1) {hwkoncost=hwkoncost1;}

	double get_hwkoncost () const {return hwkoncost;}

	void set_hwkoffcost (double hwkoffcost1) {hwkoffcost=hwkoffcost1;}

	double get_hwkoffcost () const {return hwkoffcost;}

	void set_hwkxoncost (double hwkxoncost1) {hwkxoncost=hwkxoncost1;}

	double get_hwkxoncost () const {return hwkxoncost;}

	void set_hwkxoffcost (double hwkxoffcost1) {hwkxoffcost=hwkxoffcost1;}

	double get_hwkxoffcost () const {return hwkxoffcost;}

	void set_wkoncost (double wkoncost1) {wkoncost=wkoncost1;}

	double get_wkoncost () const {return wkoncost;}

	void set_wkoffcost (double wkoffcost1) {wkoffcost=wkoffcost1;}

	double get_wkoffcost () const {return wkoffcost;}

	void set_wkxoncost (double wkxoncost1) {wkxoncost=wkxoncost1;}

	double get_wkxoncost () const {return wkxoncost;}

	void set_wkxoffcost (double wkxoffcost1) {wkxoffcost=wkxoffcost1;}

	double get_wkxoffcost () const {return wkxoffcost;}

	void set_frid (long frid1) {frid=frid1;}

	long get_frid () {return frid;}

	void set_toid (long toid1) {toid=toid1;}

	long get_toid ()  {return toid;}

	void set_ons (double ons1) {ons=ons1;}

	double get_ons () const {return ons;}

	void set_offs (double offs1) {offs=offs1;}

	double get_offs () const {return offs;}

	void set_on (short on1) {on=on1;}		// on = 1 , off = 0

	short get_on () const {return on;}

	void set_palong (double palong1) {palong=palong1;}

	double get_palong () {return palong;}

	void set_lndsf (double lndsf1) {lndsf=lndsf1;}

	double get_lndsf () {return lndsf;}

	void set_grarea (double grarea1) {grarea=grarea1;}

	double get_grarea () {return grarea;}

	void set_lvarea (double lvarea1) {lvarea=lvarea1;}

	double get_lvarea () {return lvarea;}

	void set_valtot (double valtot1) {valtot=valtot1;}

	double get_valtot () {return valtot;}

	void set_valand (double valand1) {valand=valand1;}

	double get_valand () {return valand;}

	void set_valbld (double valbld1) {valbld=valbld1;}

	double get_valbld () {return valbld;}

	void set_cfact (double cfact1) {cfact=cfact1;}

	double get_cfact () {return cfact;}

	void set_pval (double pval1) {pval=pval1;}

	double get_pval () const {return pval;}

	void set_aval (double aval1) {aval=aval1;}

	double get_aval () const {return aval;}

	void set_xc (double xc1) {xc=xc1;}

	double get_xc () const {return xc;}

	void set_yc (double yc1) {yc=yc1;}

	double get_yc () const {return yc;}

	void set_zc (double zc1) {zc=zc1;}

	double get_zc () const {return zc;}

	void set_pacid (string pacid1) {pacid=pacid1;}

	string get_pacid () {return pacid;}

	void set_ptype (string ptype1) {ptype=ptype1;}

	string get_ptype () {return ptype;}

	void set_lucd (string lucd1) {lucd=lucd1;}

	string get_lucd () {return lucd;}

	void set_pnote (string enote1) {pnote=enote1;}

	string get_pnote () {return pnote;}

	void cHOns (double sHons,double sPval) {
		if (sPval >0) {
		ons = pval/sPval * sHons;
		} else {
			ons = 0;
		}
		hwkxoncost = ons * hwkoncost;
	}  // parcel on strength
	void cHOffs (double sHoffs,double sAval) {
		if (sAval >0) {
			offs = aval/sAval * sHoffs;
		} else {
			offs = 0;
		}
		hwkxoffcost = offs * hwkoffcost;
	}  //parcel off strength

	void cOns () {
//		ons = pval/sPval * sOns;
		wkxoncost = ons * wkoncost;
	}  // parcel on strength
	void cOffs () {
//		offs = aval/sAval * sOffs;
		wkxoffcost = offs * wkoffcost;
	}  //parcel off strength

	void cPval (double sPval) {pval = pval / sPval;}
	void cAval (double sAval) {aval = aval / sAval;}

	void show_parcelhdr(ostream& out)
    { 
		out << "id" << "\t" << "pacid" << "\t" << "ptype" << "\t" << "eoid" << "\t" 
		<< "origon"<< "\t" << "origoff" << "\t" << "grarea" << "\t" << "lvarea" << "\t" 
		<< "pval" << "\t"<<"aval" << "\t" << "valbld" << "\t" << "valand"<< "\t" <<
		"valtot" << endl;
	};
 	void show_parcel(ostream& out)
    { 
		out << get_id() << "\t" << get_pacid() << "\t" << get_ptype() << "\t" << 
			get_eoid()<< "\t" << get_origon() << "\t"<< get_origoff() << "\t" << 
			get_grarea() << "\t" << get_lvarea() << "\t" << get_pval() << "\t" <<
			get_aval() << "\t" <<get_valbld() << "\t" << get_valand()<< "\t" << 
			get_valtot() << endl;
	};

	void read_parcel(istream& in)
    { 
		in >>pid>> evp>> frid>>toid>> palong>>origon>>origoff >>hwkoncost>>
			hwkoffcost>>lbl>>pnote>>eoid;
	};


	virtual ~parcel()
	{
	//	if(!ev1)
	 //{
	//	cout << "parcel Object is empty! "<<endl;}
	  //  else {cout << "parcel Object is deleted! "<<endl;
	           //delete [] ev1;
	 //  }
	}

    bool operator==(parcel par1);

// parcelOId	PID_LONG	PTYPE	LAND_SF	GROSS_AREA	LIVING_ARE	FY2003_TOT	FY2003_LAN	
// FY2003_BLD	CompFact	PVal	AVal	InOns	InOffs	OutOns	OutOffs	EdgeId

protected:
	unsigned long pid;  // parcel object id long
           	 string pacid; // Parcel id string
	         short lbl; // label
	         bool blnHist; // if this is Historic it will be 0, alternative it will be 1
	         short on; // if this is on (walking from parcel to stop) it will be 0 if off (walking from stop) it will be 1
         	 double lndsf; // land sq ft 
         	 double grarea; // gross area sq ft 
         	 double lvarea; // living area sq ft 
         	 double valtot; // total valuation
         	 double valand; // land valuation
         	 double valbld; // building valuation
         	 double cfact; // competition factor
         	 double pval; // production strength factor
         	 double aval; // attraction strength factor
			 double ons;  // historical in/outbound boardings (ons) strength stopons*pval/sum(pval)
			 double offs;  // historical in/outbound alightings (offs) strength stopoffs*aval/sum(aval)
         	 double hwkoncost; // historic offs path cost
         	 double hwkoffcost; // historic offs path cost
         	 double hwkxoncost; // historic ons path cost x ons
         	 double hwkxoffcost; // historic offs path cost x offs
			 double wkoncost; // alternative ons path cost
         	 double wkoffcost; // alternative offs path cost
         	 double wkxoncost; // alternative ons path cost x ons
         	 double wkxoffcost; // alternative offs path cost x offs
         	 double palong; // total path cost
         	 double xc; // X Coordinate
         	 double yc; // Y Coordinate
         	 double zc; // Y Coordinate
			 long eoid;  // edge id
			 long evp;  // edge pointer id
			 long frid;  // alighting stop id
			 long toid;  // boarding stop id
			 long origon;  // parcel assigned stop id
			 long origoff;  // parcel assigned stop id
        	string ptype; //Land use type
        	string lucd; //Land use code
        	string pnote; //parcel Note text
};


void parcel::parcelUp(long pid1,long eoid1,long evp1,bool blnHist, short lbl1,long frid1, long toid1,
			 long origon1=-1,long origoff1=-1, double palong1=0,short on1=-1,double lndsf1=0, 
			 double grarea1=0, double lvarea1=0, double valtot1=0,double valand1=0,
			 double valbld1=0,double cfact1=0, double pval1=0,double aval1=0,
			 string ptype1="102")
{
	
	parcel::pid = pid1;
	parcel::pacid = "-";
	parcel::eoid = eoid1;
	parcel::origon = origon1;
	parcel::origoff = origoff1;
	parcel::evp = evp1;
	parcel::lbl = lbl1;
	parcel::hwkoncost= 0;
	parcel::hwkoffcost= 0;
	parcel::frid = frid1;
	parcel::toid = toid1;
	parcel::palong = palong1;
	parcel::on = on1;
	parcel::lndsf = lndsf1; // land sq ft
    parcel::grarea = grarea1;  // gross floor area
	parcel::lvarea = lvarea1; // living area
	parcel::valtot = valtot1; // total valuation 
	parcel::valand = valand1; // land valuation 
	parcel::valbld = valbld1; // Bldg valuation 
	parcel::cfact = cfact1; // competition factor 
	parcel::pval = pval1;  // parcel production strength factor
	parcel::aval = aval1; // parcel attraction strength factor
	parcel::ons=0;
	parcel::offs=0;
	parcel::wkoncost=0;
	parcel::wkoffcost=0;
	parcel::hwkxoncost=0;
	parcel::hwkxoffcost=0;
	parcel::wkxoncost=0;
	parcel::wkxoffcost=0;
	parcel::ptype = ptype1; // parcel land use type
	parcel::lucd = "-X"; // parcel land use type
//	return parcel;
}
/*
parcel::parcel()
{
	parcel::pid = 0;
	pacid = "-";
	pfc = "-";
	parcel::eoid = 0;
	parcel::origon = 0;
	parcel::origoff = 0;
	parcel::evp = 0;
	parcel::lbl = 0;
	parcel::hwkoncost= 0;
	parcel::hwkoffcost= 0;
	parcel::hwkxoncost= 0;
	parcel::hwkxoffcost= 0;
	parcel::frid = 0;
	parcel::toid = 0;
	parcel::palong = 0;
	parcel::on = 0;
	parcel::lndsf = 0; // land sq ft
    parcel::grarea = 0;  // gross floor area
	parcel::lvarea = 0; // living area
	parcel::valtot = 0; // total valuation 
	parcel::valand = 0; // land valuation 
	parcel::valbld = 0; // Bldg valuation 
	parcel::cfact = 0; // competition factor 
	parcel::pval = 0;  // parcel production strength factor
	parcel::aval = 0; // parcel attraction strength factor
	parcel::ons=0;
	parcel::offs=0;
	wkoncost=0;
	wkoffcost=0;
    frid =-1;
	toid =-1;
    stopoff ="-1";
	stopon ="-1";
	parcel::ptype = "-1"; // parcel land use type
	parcel::lucd = "-X"; // parcel land use type
}
*/

//	parcel::ons = ons1; // parcel boardings 
//	parcel::offs = offs1; // parcel alightings

bool parcel::operator==(parcel par1)
 {
   if(pid!=par1.pid)
      return false;
   if(pacid!=par1.pacid)
      return false;
   if(frid!=par1.frid)
      return false;
   if(toid!=par1.toid)
      return false;
   return true;
 }


void parcel::serialize(ofstream& parc)
{
 parc.write(reinterpret_cast<char *>(&pid),sizeof(pid));
 size_t sizet=pacid.size();// store pacid's length
 parc.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 parc.write(pacid.c_str(), sizet+1); // write final '\0' too
 parc.write(reinterpret_cast<char *>(&lbl),sizeof(lbl));
 parc.write(reinterpret_cast<char *>(&on),sizeof(on));
 parc.write(reinterpret_cast<char *>(&lndsf),sizeof(lndsf));
 parc.write(reinterpret_cast<char *>(&grarea),sizeof(grarea));
 parc.write(reinterpret_cast<char *>(&lvarea),sizeof(lvarea));
 parc.write(reinterpret_cast<char *>(&valtot),sizeof(valtot));
 parc.write(reinterpret_cast<char *>(&valand),sizeof(valand));
 parc.write(reinterpret_cast<char *>(&valbld),sizeof(valbld));
 parc.write(reinterpret_cast<char *>(&cfact),sizeof(cfact));
 parc.write(reinterpret_cast<char *>(&pval),sizeof(pval));
 parc.write(reinterpret_cast<char *>(&aval),sizeof(aval));
 parc.write(reinterpret_cast<char *>(&ons),sizeof(ons));
 parc.write(reinterpret_cast<char *>(&offs),sizeof(offs));
 parc.write(reinterpret_cast<char *>(&hwkoncost),sizeof(hwkoncost));
 parc.write(reinterpret_cast<char *>(&hwkoffcost),sizeof(hwkoffcost));
 parc.write(reinterpret_cast<char *>(&hwkxoncost),sizeof(hwkxoncost));
 parc.write(reinterpret_cast<char *>(&hwkxoffcost),sizeof(hwkxoffcost));
 parc.write(reinterpret_cast<char *>(&wkoncost),sizeof(wkoncost));
 parc.write(reinterpret_cast<char *>(&wkoffcost),sizeof(wkoffcost));
 parc.write(reinterpret_cast<char *>(&wkxoncost),sizeof(wkxoncost));
 parc.write(reinterpret_cast<char *>(&wkxoffcost),sizeof(wkxoffcost));
 parc.write(reinterpret_cast<char *>(&palong),sizeof(palong));
 parc.write(reinterpret_cast<char *>(&xc),sizeof(xc));
 parc.write(reinterpret_cast<char *>(&yc),sizeof(yc));
 parc.write(reinterpret_cast<char *>(&zc),sizeof(zc));
 parc.write(reinterpret_cast<char *>(&eoid),sizeof(eoid));
 parc.write(reinterpret_cast<char *>(&evp),sizeof(evp));
 parc.write(reinterpret_cast<char *>(&frid),sizeof(frid));
 parc.write(reinterpret_cast<char *>(&toid),sizeof(toid));
 parc.write(reinterpret_cast<char *>(&origon),sizeof(origon));
 parc.write(reinterpret_cast<char *>(&origoff),sizeof(origoff));
 sizet=ptype.size();// store ptype's length
 parc.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 parc.write(ptype.c_str(), sizet+1); // write final '\0' too
 sizet=lucd.size();// store lucd's length
 parc.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 parc.write(lucd.c_str(), sizet+1); // write final '\0' too
 sizet=pnote.size();// store pacid's length
 parc.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 parc.write(pnote.c_str(), sizet+1); // write final '\0' too

}

void parcel::deserialize(ifstream& parc)
{
 int len=0;
 char *p=0;
 parc.read(reinterpret_cast<char *>(&pid), sizeof(pid));
 
 parc.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 parc.read(p, len+1); // copy string to temp, including '\0'
 pacid=p; // copy temp to data member
 delete[] p;


 parc.read(reinterpret_cast<char *>(&lbl),sizeof(lbl));
 parc.read(reinterpret_cast<char *>(&on),sizeof(on));
 parc.read(reinterpret_cast<char *>(&lndsf),sizeof(lndsf));
 parc.read(reinterpret_cast<char *>(&grarea),sizeof(grarea));
 parc.read(reinterpret_cast<char *>(&lvarea),sizeof(lvarea));
 parc.read(reinterpret_cast<char *>(&valtot),sizeof(valtot));
 parc.read(reinterpret_cast<char *>(&valand),sizeof(valand));
 parc.read(reinterpret_cast<char *>(&valbld),sizeof(valbld));
 parc.read(reinterpret_cast<char *>(&cfact),sizeof(cfact));
 parc.read(reinterpret_cast<char *>(&pval),sizeof(pval));
 parc.read(reinterpret_cast<char *>(&aval),sizeof(aval));
 parc.read(reinterpret_cast<char *>(&ons),sizeof(ons));
 parc.read(reinterpret_cast<char *>(&offs),sizeof(offs));
 parc.read(reinterpret_cast<char *>(&hwkoncost),sizeof(hwkoncost));
 parc.read(reinterpret_cast<char *>(&hwkoffcost),sizeof(hwkoffcost));
 parc.read(reinterpret_cast<char *>(&hwkxoncost),sizeof(hwkxoncost));
 parc.read(reinterpret_cast<char *>(&hwkxoffcost),sizeof(hwkxoffcost));
 parc.read(reinterpret_cast<char *>(&wkoncost),sizeof(wkoncost));
 parc.read(reinterpret_cast<char *>(&wkoffcost),sizeof(wkoffcost));
 parc.read(reinterpret_cast<char *>(&wkxoncost),sizeof(wkxoncost));
 parc.read(reinterpret_cast<char *>(&wkxoffcost),sizeof(wkxoffcost));
 parc.read(reinterpret_cast<char *>(&palong),sizeof(palong));
 parc.read(reinterpret_cast<char *>(&xc),sizeof(xc));
 parc.read(reinterpret_cast<char *>(&yc),sizeof(yc));
 parc.read(reinterpret_cast<char *>(&zc),sizeof(zc));
 parc.read(reinterpret_cast<char *>(&eoid),sizeof(eoid));
 parc.read(reinterpret_cast<char *>(&evp),sizeof(evp));
 parc.read(reinterpret_cast<char *>(&frid),sizeof(frid));
 parc.read(reinterpret_cast<char *>(&toid),sizeof(toid));
 parc.read(reinterpret_cast<char *>(&origon),sizeof(origon));
 parc.read(reinterpret_cast<char *>(&origoff),sizeof(origoff));

 parc.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 parc.read(p, len+1); // copy string to temp, including '\0'
 ptype=p; // copy temp to data member
 delete[] p;
 parc.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 parc.read(p, len+1); // copy string to temp, including '\0'
 lucd=p; // copy temp to data member
 delete[] p;
 parc.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for name
 parc.read(p, len+1); // copy string to temp, including '\0'
 pnote=p; // copy temp to data member
 delete[] p;
}

void parcel::serializetext(ofstream& par)
{
 par <<pid<<"\t"<<pacid<<"\t"<<lbl<<"\t"<<on<<"\t"<<lndsf<<"\t"<<grarea<<"\t"<<lvarea<<"\t"<<valtot
	 <<"\t"<<valand<<"\t"<<valbld<<"\t"<<cfact<<"\t"<<pval<<"\t"<<aval<<"\t"<<ons<<"\t"<<offs<<"\t"<<hwkoncost
	 <<"\t"<<hwkoffcost<<"\t"<<wkoncost<<"\t"<<wkoffcost<<"\t"<<palong<<"\t"<<xc<<"\t"<<yc<<"\t"<<zc<<"\t"<<eoid<<"\t"<<evp
	 <<"\t"<<frid<<"\t"<<toid<<"\t"<<origon<<"\t"<<origoff<<"\t"<<ptype<<"\t"<<lucd<<"\t"<<xc<<"\t"<<yc<<"\t"<<zc<<"\t"<<pnote<<endl; 

}

void parcel::serializetexthdr(ofstream& par)
{
 par <<"pid"<<"\t"<<"pacid"<<"\t"<<"lbl"<<"\t"<<"on"<<"\t"<<"lndsf"<<"\t"
	 <<"grarea"<<"\t"<<"lvarea"<<"\t"<<"valtot"<<"\t"<<"valand"<<"\t"<<"valbld"<<"\t"
	 <<"cfact"<<"\t"<<"pval"<<"\t"<<"aval"<<"\t"<<"ons"<<"\t"<<"offs"<<"\t"<<"hwkoncost"
	 <<"\t"<<"hwkoffcost"<<"\t"<<"wkoncost"<<"\t"<<"wkoffcost"<<"\t"<<"palong"<<"\t"<<"xc"<<"\t"
	 <<"yc"<<"\t"<<"zc"<<"\t"<<"eoid"<<"\t"<<"evp"<<"\t"<<"frid"<<"\t"<<"toid"<<"\t"<<"origon"<<"\t"<<"origoff"
	 <<"\t"<<"ptype"<<"\t"<<"lucd"<<"\t"<<"XC"<<"\t"<<"YC"<<"\t"<<"ZC"<<"\t"<<"pnote"<<endl; 

}

void parcel::texthdr(ofstream& par)
{
 par <<"pid\tpacid\tpfc\tlbl\ton\tlndsf\tgrarea\tlvarea\tvaltot\t";
 par << "valand\tvalbld\tcfact\tpval\taval\tons\toffs\t";
 par << "hwkoncost\thwkoffcost\twkoncost\twkoffcost\tpalong\txc\tyc\tzc\teoid\t";
 par << "evp\tfrid\ttoid\torigon\torigoff\tptype\tlucd\txc\tyc\tzc\tstopoff\tpnote"<< endl; 

}


// To calculate Parcel aval / pval/ons/offs totals:
class ParcelOnAccum {
	long id;
	double pval;
	double aval;
	double ons; 
	double offs; 
	double hwkons; 
	double hwkoffs; 
	double hwktmons; 
	double hwktmoffs; 
	double wkons; 
	double wkoffs; 
	double wktmons; 
	double wktmoffs; 
	double wkxoncost;
	double wkxoffcost;
	long ParCnt;
public:
  ParcelOnAccum() : id(0), pval(0), aval(0), ons(0), offs(0), 
	  hwktmons(0),  hwktmoffs(0), hwkons(0), hwkoffs(0), wktmons(0), 
	  wktmoffs(0),wkxoncost(0),wkxoffcost(0), wkons(0), wkoffs(0),ParCnt(0) {}
  void operator()(const parcel& parc) {
    id = parc.get_origon();
    pval += parc.get_pval();
    aval += parc.get_aval();
    ons += parc.get_ons();
    hwkons += parc.get_hwkoncost();
    wkons += parc.get_wkoncost();
    hwktmons += parc.get_ons()*parc.get_hwkoncost();
    wktmons += parc.get_ons()*parc.get_wkoncost();
    wkxoncost += parc.get_wkxoncost();
	ParCnt += 1;
  }

    void  AccumHeader(ofstream& os) {
    os<< "Summary of parcel attributes "<<endl;
    os<< "Stop" << "\t" << " PVal " << "\t" << " Aval "<< "\t" << "h.ons "<< "\t"<<
	"h.offs " << "\t" << " ons " << "\t" << " offs " << "\t" << "Cnt"<<"\t"<<
	"hwkoffs"<<"\t"<<"hwkons"<<"\t"<<"wkoffs"<<"\t"<<" wkons"<<"\t"<< 
	"hwktmoffs"<<"\t"<<"hwktmons"<<"\t"<<"wktmoffs"<<"\t"<<"wktmons"<<"\t"<<
	"wkxoffs"<<"\t"<<"wkxons"<< endl;
  }

  friend ostream&  operator<<(ostream& os, const ParcelOnAccum& pa) {
    /*<< "Summary of parcel attributes " <<endl
    << "Stop" << "\t" << " PVal " << "\t" << " Aval "<< "\t" << "h. ons "<< "\t"<<
	"h.In offs " << "\t" << " In ons " << "\t" << " In offs " << "\t" << "ParCnt"<< endl */
	return os <<pa.id<<"\t"<<pa.pval<<"\t"<<pa.aval<<"\t"<<pa.ons<<"\t"<<pa.offs<< "\t" 
		<<pa.ons<<"\t"<<pa.offs<< "\t"<<pa.ParCnt<<"\t"<<pa.hwkoffs<<"\t"<< pa.hwkons<<"\t"  
		<< pa.wkoffs<< "\t" << pa.wkons<< "\t" << pa.hwktmoffs<< "\t" <<pa.hwktmons<<"\t"
		<<pa.wktmoffs<<"\t"<<pa.wktmons<<"\t"<<pa.wkxoffcost<<"\t"<<pa.wkxoncost<<endl;
  }
	long get_id () {return id;}
	double get_pval() {return pval;}
	double get_aval() {return aval;}
	double get_ons() {return ons;} 
	double get_offs() {return offs;} 
	double get_hwkons() {return hwkons;} 
	double get_hwkoffs() {return hwkoffs;} 
	double get_hwktmons() {return hwktmons;} 
	double get_hwktmoffs() {return hwktmoffs;} 
	double get_wkons() {return wkons;} 
	double get_wkoffs() {return wkoffs;} 
	double get_wktmons() {return wktmons;} 
	double get_wktmoffs() {return wktmoffs;} 
	double get_wkxoncost() {return wkxoncost;} 
	double get_wkxoffcost() {return wkxoffcost;} 
	long get_ParCnt() {return ParCnt;}

};

class ParcelOffAccum {
	long id;
	double pval;
	double aval;
	double ons; 
	double offs; 
	double hwkons; 
	double hwkoffs; 
	double hwktmons; 
	double hwktmoffs; 
	double wkons; 
	double wkoffs; 
	double wktmons; 
	double wktmoffs; 
	double wkxoncost;
	double wkxoffcost; 
	long ParCnt;
public:
  ParcelOffAccum() : id(0), pval(0), aval(0),  ons(0), offs(0), 
	  hwktmons(0),  hwktmoffs(0), hwkons(0), hwkoffs(0), wktmons(0), 
	  wktmoffs(0),wkxoncost(0),wkxoffcost(0), wkons(0), wkoffs(0),ParCnt(0) {}
  void operator()(const parcel& parc) {
    id = parc.get_origoff();
    pval += parc.get_pval();
    aval += parc.get_aval();
    offs += parc.get_offs();
    hwkoffs += parc.get_hwkoffcost();
    wkoffs += parc.get_wkoffcost();
    hwktmoffs += parc.get_offs()*parc.get_hwkoffcost();
    wktmoffs += parc.get_offs()*parc.get_wkoffcost();
    wkxoffcost += parc.get_wkxoffcost();
	ParCnt += 1;
  }

  void  AccumHeader(ofstream& os) {
    os<< "Summary of parcel attributes "<<endl;
    os<< "Stop" << "\t" << " PVal " << "\t" << " Aval "<< "\t" << "h.ons "<< "\t"<<
	"h.offs " << "\t" << " ons " << "\t" << " offs " << "\t" << "Cnt"<<"\t"<<
	"hwkoffs"<<"\t"<<"hwkons"<<"\t"<<"wkoffs"<<"\t"<<" wkons"<<"\t"<< 
	"hwktmoffs"<<"\t"<<"hwktmons"<<"\t"<<"wktmoffs"<<"\t"<<"wktmons"<<"\t"<<
	"wkxoffs"<<"\t"<<"wkxons"<< endl;
  }
  friend ostream&  operator<<(ostream& os, const ParcelOffAccum& pa) {
    return os <<pa.id<<"\t"<<pa.pval<<"\t"<<pa.aval<<"\t"<<pa.ons<<"\t"<<pa.offs<< "\t" 
		<<pa.ons<<"\t"<<pa.offs<< "\t"<<pa.ParCnt<<"\t"<<pa.hwkoffs<<"\t"<< pa.hwkons<<"\t"  
		<< pa.wkoffs<< "\t" << pa.wkons<< "\t" << pa.hwktmoffs<< "\t" <<pa.hwktmons<<"\t"
		<<pa.wktmoffs<<"\t"<<pa.wktmons<<"\t"<<pa.wkxoffcost<<"\t"<<pa.wkxoncost<<endl;
  }
	long get_id () {return id;}
	double get_pval() {return pval;}
	double get_aval() {return aval;}
	double get_ons() {return ons;} 
	double get_offs() {return offs;} 
	double get_hwkons() {return hwkons;} 
	double get_hwkoffs() {return hwkoffs;} 
	double get_hwktmons() {return hwktmons;} 
	double get_hwktmoffs() {return hwktmoffs;} 
	double get_wkons() {return wkons;} 
	double get_wkoffs() {return wkoffs;} 
	double get_wktmons() {return wktmons;} 
	double get_wktmoffs() {return wktmoffs;} 
	double get_wkxoncost() {return wkxoncost;} 
	double get_wkxoffcost() {return wkxoffcost;} 
	long get_ParCnt() {return ParCnt;}

};


class EdgeAccum {
			 long id;  // edge id
	         long eoid; // edge predecessor
	         long lbl; // label
	         long dirn; // forward = 1, reverse=-1, boundary=-2,unscanned=0;
			 long	toStop; // 1 for cost of travel to a stop (ons) 0 otherwise (offs) 
			 double ecost; // edge cost
         	 double scost; // start path cost
         	 double tcost; // total path cost
         	 double slen; // length of shape
         	 double palong; // total path cost
			 long orig;  // edge to id
			 long tway;  // two way or one way indicator
        	string enote; //Edge Note pointer
			long Egcnt; // count 
public:
  EdgeAccum() : id(0), orig(0), ecost(0),scost(0),tcost(0), lbl(0),
	  dirn(1),Egcnt(0),palong(0),slen(0),eoid(-1),toStop(0),tway(1),enote("") {}
  void operator()(const edgev& Eg) {
	id = Eg.get_id();
	orig = Eg.get_orig();
	dirn += Eg.get_dirn();
	lbl += Eg.get_lbl();
    palong += Eg.get_palong();
    ecost += Eg.get_cost();
    scost += Eg.get_scost();
    tcost += Eg.get_tcost();
    toStop += Eg.get_toStop();
    slen += Eg.get_slen();
    tway += Eg.get_tway();
	Egcnt += 1;
  }

  friend ostream&  operator<<(ostream& os, const EdgeAccum& EgAc) {
	/*  << "Summary of Vertex attributes by id " <<endl
    << "Id" << "\t"<< "Stop" << "\t" << " id " << "\t" << " cost "<< "\t" << "TotCost "
	<< "\t"<< "ParCnt"<< "\t"<<"index "<< "\t"<<"LowLink "<< endl */
    return os << EgAc.id<< "\t" << EgAc.orig << "\t" << EgAc.ecost << "\t" << EgAc.scost 
		<< "\t" << EgAc.tcost <<  "\t" << EgAc.Egcnt<< "\t" << EgAc.dirn <<"\t" << EgAc.lbl 
		<<  endl;
  }

	long get_id () {return id;}

	long get_orig () {return orig;}

	long get_eoid () {return eoid;}

	long get_lbl () {return lbl;}

	long get_dirn () {return dirn;}

	double get_cost () {return ecost;}

	long get_tway () {return tway;}

	double get_tcost () {return tcost;}

	double get_scost () {return scost;}

	double get_slen () {return slen;}

	string get_enote () {return enote;}

	long get_toStop () {return toStop;}

	void show_edgehdr(ostream& out)
    { 
		out << "eid" << "\t" << "orig" << "\t" << "ecost"<< "\t" << "scost"<< 
		"\t" << "tcost" << "\t" << "dirn" << "\t" << "lbl" << endl;
	};
 	void show_edge(ostream& out)
    { 
		out << get_id() << "\t" << get_orig() << "\t" << get_cost() << "\t" 
		<< get_scost() << "\t" <<get_tcost() << "\t" << get_dirn() << "\t" 
		<< get_lbl() << endl;
	};
	long get_Egcnt() {return Egcnt;} 

};




class pointdp {

 public:
	pointdp(long pnode1=9999999, double cost1=100000.0);

	void set_pn (unsigned long p1) {pnode=p1;}
	unsigned long get_pn() 	{ return pnode;	}

	void set_cost (double cost1) {cost=cost1;}
    double get_cost() { return cost; }

	virtual ~pointdp(){ 
		//	cout << "Point DP Object is deleted! "<<endl;
			//delete[] pd;
	}

protected:
	unsigned long pnode;
	double cost;
	pointdp* pd;

};

pointdp::pointdp(long pnode1, double cost1)
{
pnode = pnode1;
cost=cost1;
}

class nodedp {
public:
		nodedp(long nodeid1, short lbl1);

	void set_ndid (unsigned long nodeid1) {nodeid=nodeid1;}
	unsigned long get_pn() 	{return nodeid;}

	void set_lbl (short lbl1) {lbl=lbl1;}
    short get_lbl() { return lbl; }

	virtual	~nodedp(){delete[] np;
	//		cout << "Point Node DP Object is deleted! "<<endl;
			}
protected:
	unsigned long nodeid;
	         short lbl;
    nodedp* np;
};

nodedp::nodedp(long nodeid1=0, short lbl1=0)
{
nodeid = nodeid1;
lbl=lbl1;
}

template <class Tg>
class genar {
public:

	Tg& operator[]( int );
	genar( int );
    
	~genar() {delete[] a;}

	int get_size() { return size; }

private:
	Tg* a;
	int size;
	genar();
};

template <class Tg>
genar<Tg>::genar (int s) {
		a = new Tg[ size = s ];
	}

template <class Tg>
Tg& genar<Tg>::operator[] (int i) {
		if (i < 0 || i >= size ) {
			cerr << "index " << i << "out of bounds : ";
			return a[0];
		}
		return a[i];
	}


 
template <class Tg>
ostream& operator<< (ostream& out, genar< Tg >& a) {
		for ( int i = 0; i < a.get_size(); i++) 
			out << a [i] << endl; 
			return out;
	}



template <class T2d>
class dynamic_2d_array
{
public:
  dynamic_2d_array(int row, int col) : m_row(row),
                                       m_col(col),
                                       m_data((row != 0 && col != 0) ? new T2d[row * col] : NULL){}

  dynamic_2d_array(const dynamic_2d_array& src) : m_row(src.m_row),
                                                  m_col(src.m_col),
                                                  m_data((src.m_row != 0 && src.m_col != 0) ? new T2d[src.m_row * src.m_col] : NULL)
  {
    for(int r = 0;r < m_row; ++r)
      for(int c = 0; c < m_col; ++c)
        (*this)[r][c] = src[r][c];
  }

  virtual ~dynamic_2d_array()
  {
    if(m_data)
      delete []m_data;
  }
  
  inline T2d* operator[](int i) { return (m_data + (m_col * i)); }

  inline T2d const*const operator[](int i) const {return (m_data + (m_col * i)); }


private:
  dynamic_2d_array& operator=(const dynamic_2d_array&);
  const int m_row;
  const int m_col;
  T2d* m_data; 
};


template<typename T> const T& minT(const T& a, const T& b) {
  return (a < b) ? a : b;
}

template<typename T> const T& maxT(const T& a, const T& b) {
  return (a >= b) ? a : b;
}

template<typename Tx> const Tx& minTx(const Tx& a, const Tx& b) {
  return (a.get_cost() < b.get_cost()) ? a : b;
}


template<class T> class objCount {
  static int count;
public:
  objCount() { ++count; }
  objCount(const objCount<T>&) { ++count; }
  ~objCount() { --count; }
  static int getCount() { return count; }
};
 
template<class T> int objCount<T>::count = 0;


//ifstream& openfile(ifstream& infile);
// function prototypes
ifstream& openfile(ifstream& infile1);
const unsigned MaxStrLen = 200;

class inputfilename {
public:
	inputfilename(); //constructor
    inputfilename(const char iname[]); //constructor
	~inputfilename();
	void set(const char string[]);
	const char* get();
private:
	char iname[MaxStrLen + 1]; // +1 for \0
};
// define set method for inputfilename
void inputfilename::set(const char string[]) {
	::strncpy_s (iname, string, MaxStrLen);
	iname[MaxStrLen] = '\0';
}
// define get method for inputfilename
const char* inputfilename::get() {
	return iname;
}



//: C05:Sortable.h
// Template specialization.
#ifndef SORTABLE_H
#define SORTABLE_H
 
template<class T1, class T2>
class mapsort : public std::map<T1,T2> 
{
public:
  void sort();
};
 
/*
// Partial specialization for pointers:
template<class T1*, class T2*>
class mapsort<T1 , T2> : public std::map<T1 , T2> 
{
public:
  void sort();
};
*/ 
template<class T1, class T2>
void mapsort<T1, T2>::sort() 
{
  std::map<T1, T2>::iterator iter;
  temp = (*iter).begin;
  for (iter = .begin(); iter != .end(); iter++)
  {
      if ((*iter)->get_cost >  temp->get_cost)
      tempend->next = p;
      end = p;
  }
}

#endif // SORTABLE_H ///:~
// map container sorter
template <typename container>
void sort(const container& c)
{
  double x,y;
  container::const_iterator it, temp;
  temp = c.begin;
  for(container::const_iterator it=c.begin();it!=c.end();++it)
  {
    x = temp->last->get_cost();
    y = it->last->get_cost();
	if (x>y) {
     c.push_down(temp);
	  cout<<it->first<<""<<it->last<<endl;
	}
  }
}

template <typename m, typename n, typename o, typename p >
n& sortEdge(m& m1, n& n1, o& o1, p& p1)
{
  double x=0,y=0;
  long frid=0, toid=0,pred=0,idx=0;
  m::iterator it;
  m::iterator temp;
  typedef pair <p,o> obj_pair;
  for(it=m1.begin();it!=m1.end();++it)
  {
    frid = it->second.get_frid();
    toid = it->second.get_toid();
		temp = m1.find(toid);
		if (temp!=m1.end())
		{
			idx = idx+1;
			temp->second.set_evp(it->second.get_frid());
			n1.insert(obj_pair(idx,it->second));
		} else {
			p1 = inf;
			n1.insert(obj_pair(p1,it->second));
		}

  }
  return n1;
}

template <typename m,  typename o,typename n, typename p >
n& sortEdgeObj( m& m1,o& o1,n& n1,p& p1)
{
	double x=0,y=0;
	static long idx=0;
	long frid=0, toid=0,pred=0;
	m::iterator it;
	m::iterator iu;
	o* pO=&o1;
	o o2;
	m::iterator temp;
	typedef pair <p,o> obj_pair;
	n1.insert(obj_pair(idx,o1));
    frid = o1.get_frid();
    toid = o1.get_toid();
	temp = m1.find(toid);
		if (temp!=m1.end())
		{
			idx = idx+1;
			o2 = temp->second;
			m1.erase(temp);
			// recusively update the order of edges for polyines (multi-lines)
			sortEdgeObj(m1,o2,n1,p1);
		}

  return n1;
}


// Delete pointers in an STL sequence container.
#ifndef PURGE_H
#define PURGE_H
#include <algorithm>
 
template<class Seq> void purge(Seq& c) {
  typename Seq::iterator i;
  for(i = c.begin(); i != c.end(); ++i) {
    delete *i;
    *i = 0;
  }
}
 
// Iterator version:
template<class InpIt> void purge(InpIt begin, InpIt end) {
  while(begin != end) {
    delete *begin;
    *begin = 0;
    ++begin;
  }
}
#endif // PURGE_H ///:~

// update the stop attributes using accumulator deque
#ifndef UPDATESTOP_H
#define UPDATESTOP_H
#include <algorithm>

template <typename v , typename s,  typename o>
s& updStopDmndWlkTm(v& vpa, s& s1, o& o1,bool blnHist=true)
{
  v::iterator vpait, vtemp;
  s::iterator si1, stemp;
  long ip;
  tstop* pstop;
  for(vpait=vpa.begin();vpait!=vpa.end();vpait++)
  {
	  
    ip = (vpait->get_id());
	if (ip>0) {
		si1= s1.find(ip);
		 pstop = &(si1->second);
	 if (o1) {
		 if (blnHist) {
			pstop->set_Ons(pstop->get_HistOns());
			pstop->set_hWkTmOns(vpait->get_hwktmons());
			pstop->set_WkTmOns(vpait->get_hwktmons());
			pstop->set_PVal(vpait->get_pval());
		 } else 
		 {
			pstop->set_Ons(vpait->get_ons());
			pstop->set_WkTmOns(vpait->get_wktmons());
		 }
	 }
	 else {
		 if (blnHist) {
			pstop->set_Offs(pstop->get_HistOffs());
			pstop->set_hWkTmOffs(vpait->get_hwktmoffs());
			pstop->set_WkTmOffs(vpait->get_hwktmoffs());
			pstop->set_AVal(vpait->get_aval());
		 } else
		 {
			pstop->set_Offs(vpait->get_offs());
			pstop->set_WkTmOffs(vpait->get_wktmoffs());
			}
	 }
	} // if a valid stop id (>0) is present 

   }
  return  s1;
}


#endif // UPDATESTOP_H ///:~

// update the stop attributes using accumulator deque
#ifndef CALCPARCELDMND_H
#define CALCPARCELDMND_H
#include <algorithm>

template < typename s, typename p,  typename o>
p& calcParcelDemand0(s& s1,p& mmaparStop, o& o1,bool blnHist=true)
{

s::iterator si1;
typedef p::iterator parit;
parit parit1;
long ip;
tstop* pstop;

  for(si1=s1.begin();si1!=s1.end();si1++)
  {
	pstop = &(si1->second);
	ip = pstop->get_id();
// get the range of parcel objects in the parcel-stop map with current stop id
	 	pair<parit, parit> pastrange = mmaparStop.equal_range(ip);

		calcParcelDemand1(pastrange.first, pastrange.second, pstop, o1, blnHist);
  }
return mmaparStop;
}

template <typename s, typename p, typename o>
p& calcParcelDemand1(s& s1,p& mmaparStop, o& o1,bool blnHist=true)
{
s::iterator si1, stemp;
typedef p::iterator parit, ptemp;
parit parit1;
long ip;
   parcel* parp;
   tstop* pstop;

  for(si1=s1.begin();si1!=s1.end();si1++)
  {
	pstop = &(si1->second);
	ip = pstop->get_id();
// get the range of parcel objects in the parcel-stop map with current stop id
	 	pair<parit, parit> pastrange = mmaparStop.equal_range(ip);
		for (parit1 = pastrange.first; parit1!=pastrange.second;parit1++)
	    {
	   	   parp = &(parit1->second);
		   // calculate  the fraction of demand (boarding/aligting) values for the parcel

			if (o1) {
				if (blnHist) {
					parp->cHOns(pstop->get_HistOns(),pstop->get_PVal()); //parp->set_hinons(parp->get_pval()* pstop->get_HistOns());
				} else
				{
					parp->cOns(); //parp->set_hinoffs(parp->get_aval()* pstop->get_HistOffs());
				}
			} else
			{
				if (blnHist) {
					parp->cHOffs(pstop->get_HistOffs(),pstop->get_AVal()); //parp->set_hinoffs(parp->get_aval()* pstop->get_HistOffs());
				} else
				{
					parp->cOffs(); //parp->set_hinoffs(parp->get_aval()* pstop->get_HistOffs());
				}
				
			}
		}
  }

  return  mmaparStop;
}

// Iterator version:
template<class parIt,typename s, typename o> 
void calcParcelDemand1(parIt begin, parIt end,s& pstop,o& o1, bool blnHist) 
{
  parIt parIt1;
  parIt1 = begin;
  parcel* parp;

  while(parIt1 != end) {
// get the range of parcel objects in the parcel-stop map with current stop id
	   	   parp = &(parIt1->second);
		   // calculate  the fraction of demand (boarding/aligting) values for the parcel

			if (o1) {
				if (blnHist) {
					parp->cHOns(pstop->get_HistOns(),pstop->get_PVal()); //parp->set_hinons(parp->get_pval()* pstop->get_HistOns());
				} else {
					parp->cOns(); //parp->set_hinons(parp->get_pval()* pstop->get_HistOns());
				}
			} else {
				if (blnHist) {
					parp->cHOffs(pstop->get_HistOffs(),pstop->get_AVal()); //parp->set_hinoffs(parp->get_aval()* pstop->get_HistOffs());
				} else {
					parp->cOffs(); //parp->set_hinoffs(parp->get_aval()* pstop->get_HistOffs());
				}
			}
			parIt1++;
  }
}

// Iterator version to find the edge object that the parcel is located at
// p - parcel object
// parIt - range of edge and eoid mapping for an edge
// e1 collection of edges ordered by eoid
// o - edge object
// v1 collection of vertices
// x - vertex object
template<typename p,class parIt,typename e, typename o,typename v,typename x> 
o& findEdgePart(p& parc, parIt begin, parIt end,e& e1, o& o1, v& v1, x& x1) 
{
	parIt ePartIt;
	ePartIt = begin;
	long id=0,eid=0,vid1=0,vid2=0;
	double vcost=0,xcost=0,posalong=0,xdist=0,ecost=0,epalong=0;
	bool bnEdge=false;
	v::iterator vit;
	v::iterator vit2;
	typedef e::iterator eit;
	eit eit1;
	eit eit2;
	o* pO=&o1;
	o o2;
	typedef pair<long,o> ePair;
	e e2;
	e e3;
	posalong = parc->get_palong();
		while(ePartIt != end) {
			// sum the cost of the edge pieces to determine where the object is located 
			// and find out the two vertices that it belongs. If j = 1 then there is no need to do this

	// get the edge objects in the range with current edge id
			o1 = (ePartIt->second);
				vid1 = o1.get_toid();
				e2.insert(ePair(vid1,o1));
			ePartIt++;
		}
		for (eit1=e2.begin();eit1!=e2.end();++eit1)
		{// look for the relation between to and from for successive objects
			eit2 = e2.find(eit1->second.get_frid());
			if ( eit2== e2.end())
			{// this is the begining of the object collection
				o2=eit1->second;
			}
		}
		for (eit1=e2.begin();eit1!=e2.end();++eit1)
		{
			o1 = eit1->second;
			e3.insert(ePair(o1.get_frid(),o1));
		}
		e2.clear();
		e2 = sortEdgeObj(e3,o2,e2,id);
		xcost = 0;
		for (eit1=e2.begin();eit1!=e2.end();++eit1)
		{
			pO= &(eit1->second);
			xcost = xcost + pO->get_cost();
		}
		eit1 = e2.begin();
		while (eit1!=e2.end()) 
		{
			pO= &(eit1->second);
			epalong = epalong + pO->get_cost()/xcost;
			pO->set_palong( epalong );
			if (posalong<epalong) { // this edge is the requested link
				epalong = (pO->get_cost()-(epalong-posalong)*xcost)/pO->get_cost();
				pO->set_palong( epalong );
				bnEdge=true;
				break;
			}
			eit1++;
		}
		pO->set_toStop(bnEdge);
return *pO;
}

#endif // CALCPARCELDMND_H ///:~

template <typename m , typename o, typename f>
m& VertexReport(m& m1, o& o1, f& f1)
{
m::iterator mit;
long vid;
int il=0,j=0;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     vid = mit->first;
     o1 = &(mit->second);
// get the vertex id & vertex object
//		k =(*mapverit).first;
//		vxp =  &(mapverit->second);
		if (o1->get_lbl()==0) { il++;}
		if (o1->get_orig()>0) { 
			if (o1->get_tcost()==inf && o1->get_cost()<inf )
			{o1->set_tcost(o1->get_cost());}
		j++;}
  }
  f1 << " End Vertex Export. "<<j<<" of " <<m1.size()<< " veritices assigned "<<il<<" Were not assigned to any stop."<<endl;
  return m1;
}


#ifndef ASSIGNPARCSTOP_H
#define ASSIGNPARCSTOP_H
#include <algorithm>

template <typename v , typename e, typename p, typename h>
p& assignParcStop(v& vx1, v& vx2, e& evp, p& parp,h& h1)
{
  double c1=0,v1=vx1->get_cost(), v2=vx2->get_cost(),x1=evp->get_cost(),posa=parp->get_palong();
  bool forwd=true;
  // find the indifference point and see which side the parcel is located
  if (v1<v2) {
	forwd = true ;
  } else {
	forwd = false; 
  }
	c1 = (abs(v2-v1)+x1)/(2*x1); // from vx1
//	if (forwd) { //  Cost of v1 < v2
		if (posa<c1) // this is closer to vx1
		{ // posalong < indifference point v1 is governing (forwd)
			if (parp->get_on()==0) 
			{ // offs
				parp->set_origoff(vx1->get_orig());
				if (h1) {
					parp->set_hwkoffcost(v1+x1*c1);
				} else {
					parp->set_wkoffcost(v1+x1*c1);
				}
			} else 	{ // ons 
				if (h1) {
					parp->set_hwkoncost(v1+x1*(c1));
				} else {
					parp->set_wkoncost(v1+x1*(c1));
				}
				parp->set_origon(vx1->get_orig());
			}
		} else { // posalong > indifference point, hence closer to vx2 (forwd)
			if (parp->get_on()==0) 
			{ // offs
				parp->set_origoff(vx2->get_orig());
				if (h1) {
					parp->set_hwkoffcost(v2+x1*(1-posa));
				} else {
					parp->set_wkoffcost(v2+x1*(1-posa));
				}
			} else 	{ // ons 
				if (h1) {
					parp->set_hwkoncost(v1+x1*(1-posa));
				} else {
					parp->set_wkoncost(v1+x1*(1-posa));
				}
				parp->set_origon(vx2->get_orig());
			}
		}
//	} 
		/* else { // backwd
		if (posa<c1) // this is closer to vx2
		{ // posalong < indifference point v2 is governing (backwd)
			if (parp->get_on()==0) 
			{ // offs
				parp->set_origoff(vx2->get_orig());
				if (h1) {
					parp->set_hwkoffcost(v2+x1*c1);
				} else {
					parp->set_wkoffcost(v2+x1*c1);
				}
			} else 	{ // ons 
				if (h1) {
					parp->set_hwkoncost(v2+x1*(c1));
				} else {
					parp->set_wkoncost(v2+x1*(c1));
				}
				parp->set_origon(vx2->get_orig());
			}
		} else { // posalong > indifference point, hence closer to vx1 (backwd)
			if (parp->get_on()==0) 
			{ // offs
				parp->set_origoff(vx1->get_orig());
				if (h1) {
					parp->set_hwkoffcost(v1+x1*(1-posa));
				} else {
					parp->set_wkoffcost(v1+x1*(1-posa));
				}
			} else 	{ // ons 
				if (h1) {
					parp->set_hwkoncost(v1+x1*(1-posa));
				} else {
					parp->set_wkoncost(v1+x1*(1-posa));
				}
				parp->set_origon(vx1->get_orig());
			}
		}
	} */
  return  parp;
}

template <typename v , typename e, typename p, typename h>
p& assignParcStopx(v& vx1, v& vx2, e& evp, p& parp,h& h1)
{
  double c1=0,x1=0,x2=0; // 
	long k=0;
	bool forwd = (vx2->getcost() - vx1->getcost())>=0;
	if (forwd) {
		x1=vx1->get_cost();
		x2=vx2->get_cost();
		if (c1 < parp->get_palong()) // this is closer to vx2
		{
			k = vx2->get_orig();
		} else {
			k = vx1->get_orig();
		}
	} else {
		x1=vx2->get_cost();
		x2=vx1->get_cost();
		if (c1 < parp->get_palong()) // this is closer to vx1
		{
			k = vx2->get_orig();
		} else {
			k = vx1->get_orig();
		}
	}
	
	if (c1 > parp->get_palong()) // this is closer to vx1
	{
		if (parp->get_on()==0) 
		{
			parp->set_origoff(vx1->get_orig());
			if (h1) {
				parp->set_hwkoffcost(vx1->get_cost()+evp->get_cost()*c1);
			} else {
				parp->set_wkoffcost(vx1->get_cost()+evp->get_cost()*c1);
			}

		}
		else
		{
			if (h1) {
				parp->set_hwkoncost(vx1->get_cost()+evp->get_cost()*c1);
			} else {
				parp->set_wkoncost(vx1->get_cost()+evp->get_cost()*c1);
			}
			parp->set_origon(vx1->get_orig());
		}
	}
	else  // closer to vx2
	{
		if (parp->get_on()==0) {
			if (h1) {
				parp->set_hwkoffcost(vx2->get_cost()+evp->get_cost()*c1);
			} else {
				parp->set_wkoffcost(vx2->get_cost()+evp->get_cost()*c1);
			}
			parp->set_origoff(vx2->get_orig());
		}
		else
		{
			if (h1) {
				parp->set_hwkoncost(vx2->get_cost()+evp->get_cost()*c1);
			} else {
				parp->set_wkoncost(vx2->get_cost()+evp->get_cost()*c1);
			}
			parp->set_origon(vx2->get_orig());
		}
	}

}
#endif // ASSIGNPARCSTOP_H ///:~


#ifndef STOPASSIGN_H
#define STOPASSIGN_H

class StopAssign {
	long id;
	long cnt;
	long serNo;
public:
  StopAssign() : id(0), cnt(0), serNo(0) {}
  void operator()(const tstop& ts) {
    serNo +=1;
  }
  void set_id 	(long id1) {id=id1;}
  void set_cnt 	(long cnt1) {cnt=cnt1;}
  long get_id() { return id; }
  long get_cnt() { return cnt; }
};

class objAssign {
	long id;
	long cnt;
	long serNo;
public:
  objAssign() : id(0), cnt(0), serNo(0) {}
  template <typename o>
  void operator()(const o& obj) {
    serNo +=1;
  }
  void set_id 	(long id1) {id=id1;}
  void set_cnt 	(long cnt1) {cnt=cnt1;}
  long get_id() { return id; }
  long get_cnt() { return cnt; }
};


template <typename Seq1, typename Seq2, typename Seq3>  
  Seq3& objsCntWithId(const Seq1& c1,Seq2& c2, Seq3& c3) 
  {
  typename Seq1::iterator sit1;
  typename Seq2::iterator sit2;
  StopAssign sa;
	for(sit1 = c1.begin(); sit1 != c1.end(); ++sit1) {
      pair<sit2, sit2> objrange = c2.equal_range((sit1->second).get_id());
	  size_t j = distance(objrange.first,objrange.second);
	
	  sa.set_id((sit1->second).get_id());
	  sa.set_cnt(j);
	  c3->push_back(sa);
	}
  }


  template <typename Seq1, typename Seq2, typename Seq3>  
  Seq3& objsWithId(const Seq1& c1,Seq2& c2, Seq3& c3) 
  {
  typename Seq1::iterator sit1;
  typename Seq2::iterator sit2;
  StopAssign sa;
	for(sit1 = c1.begin(); sit1 != c1.end(); ++sit1) {
      pair<sit2, sit2> objrange = c2.equal_range((sit1->second).get_id());
	  size_t j = distance(objrange.first,objrange.second);
	
	  sa.set_id((sit1->second).get_id());
	  sa.set_cnt(j);
	  c3->push_back(sa);
	}
  }

  template <typename Seq1, typename Seq2, typename Seq3>  
  Seq3& objsWithOrig(const Seq1& c1,Seq2& c2, Seq3& c3) 
  {
  typename Seq1::iterator sit1;
  typename Seq2::iterator sit2;
  StopAssign sa;
	for(sit1 = c1.begin(); sit1 != c1.end(); ++sit1) {
      pair<sit2, sit2> objrange = c2.equal_range((sit1->second).get_id());
	  size_t j = distance(objrange.first,objrange.second);
	  sa.set_id((sit1->second).get_id());
	  sa.set_cnt(j);
	  c3->push_back(sa);
	}
  }


#endif // STOPASSIGN_H ///:~

// Calculate undelayed time in an STL sequence container.
#ifndef UNDELTM_H
#define UNDELTM_H
#include <algorithm>
 
template<class Seq> void undelTm(Seq& c) {
  typename Seq::iterator psit,psit2;
  tstop* pstop;
  psit2 = c.begin();
  for(psit = c.begin(); psit != c.end(); ++psit) {
		pstop = &(psit->second);
	    pstop->set_undCRdTm(pstop->get_CRdTm());
	    pstop->set_CRdTmC(pstop->get_CRdTm());
  }
}


// Iterator version:
template<class InpIt> void undelTm(InpIt begin, InpIt end,double unOn, double unOff, double hdWay) {
  tstop* pstop;
  tstop* pstop0;
  InpIt stIt;
  int ordr=0;
  double t2=0,t3=0,t4=0,t5=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
  stIt = begin;
  while(stIt != end) {
		pstop = &(stIt->second);
		pstop->cprobstoph(hdWay);
		if (begin==stIt) {
			pstop0=pstop;
			pstop->set_dwellDelay(0);
			t2 =  pstop->segDelay(1,0.0,0.0);
			pstop->set_undCRdTm(pstop->get_CRdTm());
			pstop->set_CRdTmC(pstop->get_CRdTm());
			pstop->set_probStoph(1.0);
			pstop->set_probStop(1.0);
		}
		else
		{
			pstop->cdwlDelayh(unOn,unOff,hdWay);
			t2 = pstop->segDelay(pstop0->get_probStoph(),pstop0->get_depDelay(),pstop->get_probStoph());
			pstop->set_undCRdTm(pstop0->get_undCRdTm()  + (pstop->get_CRdTm() - pstop0->get_CRdTm() - t2/60 ));
			if (pstop->get_undCRdTm()<0) {
			cout<<" Negative Undelayed ride time = "<<pstop->get_undCRdTm()<<endl<<
				"  Please review ride time data between stops "<<pstop0->get_StopName()<<
				" and "<<pstop->get_StopName()<<endl;
			}
		}
		if (pstop->get_posalong()<0) {
			pstop->set_posalong(0);
		}
		if (pstop->get_StOrdr()<=0) {
			pstop->set_StOrdr(pstop->get_id());
		}
		if (pstop->get_Stopidp()<=0) {
			pstop->set_Stopidp(pstop0->get_id());
		}
			pstop0=&(stIt->second);
		stIt++;
  }
}


// Iterator version:
template<typename m, class InpIt> 
m& undelTMap(m& m1,double unOn, double unOff, double hdWay, std::ofstream &logfile) {
  tstop* pstop;
  tstop* pstop0;
  InpIt stIt;
  InpIt begin = m1.begin();;
  InpIt end= m1.end();

  int ordr=0;
  double t2=0,t3=0,t4=0,t5=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
  stIt = begin;
  while(stIt != end) {
		pstop = &(stIt->second);
		pstop->cprobstoph(hdWay);
		if (begin==stIt) {
			pstop0=pstop;
			pstop->set_dwellDelay(0);
			t2 =  pstop->segDelay(1,0.0,0.0);
			pstop->set_undCRdTm(pstop->get_CRdTm());
			pstop->set_CRdTmC(pstop->get_CRdTm());
			pstop->set_probStoph(1.0);
			pstop->set_probStop(1.0);
		}
		else
		{
			pstop->cdwlDelayh(unOn,unOff,hdWay);
			t2 = pstop->segDelay(pstop0->get_probStoph(),pstop0->get_depDelay(),pstop->get_probStoph());
			pstop->set_undCRdTm(pstop0->get_undCRdTm()  + (pstop->get_CRdTm() - pstop0->get_CRdTm() - t2/60 ));
			if (pstop->get_undCRdTm()<0) {
			cout<<" Negative Undelayed ride time = "<<pstop->get_undCRdTm()<<endl<<
				"  Please review ride time data between stops "<<pstop0->get_StopName()<<
				" and "<<pstop->get_StopName()<<endl;
			logfile<<" Negative Undelayed ride time = "<<pstop->get_undCRdTm()<<endl<<
				"  Please review ride time data between stops "<<pstop0->get_StopName()<<
				" and "<<pstop->get_StopName()<<endl;
			
			}
		}
		if (pstop->get_posalong()<0) {
			pstop->set_posalong(0);
		}
		if (pstop->get_StOrdr()<=0) {
			pstop->set_StOrdr(pstop->get_id());
		}
		if (pstop->get_Stopidp()<=0) {
			pstop->set_Stopidp(pstop0->get_id());
		}
			pstop0=&(stIt->second);
		stIt++;
  }
return m1;
}


// Iterator version for individual trips :
template<typename m, class InpIt> 
m& undelTripMap(m& m1,double unOn, double unOff, double hdWay, std::ofstream &logfile) {
  tstop* pstop;
  tstop* pstop0;
  InpIt stIt;
  InpIt begin = m1.begin();;
  InpIt end= m1.end();

  int ordr=0;
  double t2=0,t3=0,t4=0,t5=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
  stIt = begin;
  while(stIt != end) {
		pstop = &(stIt->second);
		pstop->cprobstoph(hdWay);
		pstop->cprobstop(hdWay);
		ordr++;  // increment the stop count
		if (begin==stIt) {
			pstop0=pstop;
			pstop->set_dwellDelay(0);
			t2 =  pstop->segDelay(1,0.0,0.0);
			pstop->set_undCRdTm(pstop->get_CRdTm());
			pstop->set_CRdTmC(pstop->get_CRdTm());
			pstop->set_probStoph(1.0);
			pstop->set_probStop(1.0);
		}
		else
		{
			pstop->cdwlDelayh(unOn,unOff,hdWay);
			t2 = pstop->segDelay(pstop0->get_probStoph(),pstop0->get_depDelay(),pstop->get_probStoph());
			if (ordr==1) { // ignore the dwell time at the first stop
				pstop->set_undCRdTm(pstop0->get_undCRdTm()  + (pstop->get_CRdTm() - pstop0->get_CRdTm()  ));
			} else { 
				pstop->set_undCRdTm(pstop0->get_undCRdTm()  + (pstop->get_CRdTm() - pstop0->get_CRdTm() - t2/60 ) ); // pstop0->get_dwellDelay()/60 - pstop0->get_depDelay()/60*pstop0->get_probStop() - pstop->get_arrDelay()/60*pstop->get_probStop()));
			}
			if (pstop->get_undCRdTm()<0) {
				cout<<" Negative Undelayed ride time = \t"<<pstop->get_undCRdTm()<<endl<<
					"  Please review ride time data between stops "<<pstop0->get_StopName()<<
					" and "<<pstop->get_StopName()<<"\t dwell time (mins) : "<<pstop0->get_dwellDelay()<<
					"\t Arrival Delay (Secs) : "<<pstop->get_arrDelay()<<"\t Dep Delay (Secs) : "<<pstop0->get_depDelay()
					<<"\t and Stop Number : "<<ordr<<endl;
				logfile<<" Negative Undelayed ride time = \t"<<pstop->get_undCRdTm()<<endl<<
					"  Please review ride time data between stops : \t"<<pstop0->get_StopName()<<
					" and \t"<<pstop->get_StopName()<<"\t dwell time :"<<pstop0->get_dwellDelay()
					<<"\tStop Number"<<ordr<<endl;
			}
		}
		if (pstop->get_posalong()<0) {
			pstop->set_posalong(0);
		}
		if (pstop->get_StOrdr()<=0) {
			pstop->set_StOrdr(pstop->get_id());
		}
		if (pstop->get_Stopidp()<=0) {
			pstop->set_Stopidp(pstop0->get_id());
		}
			pstop0=&(stIt->second);
		stIt++;
  }
return m1;
}


#endif // UNDELTM_H ///:~

#ifndef IMPACTCALC_H //Walk, Ride, Operating cost
#define IMPACTCALC_H
#include <algorithm>

// Iterator version:
template<class InpIt> 
void ImpactCalc(InpIt begin, InpIt end,globcost& gc, pdhdway& pdh,bool blnHist) {
	//initialize the stop variables needed for the program
    tstop stop0;
	tstop stop1;
	tstop stop2;
    tstop* pstop0;
	tstop* pstop;
	tstop* pstop2;
    pstop0 = &stop0;
    pstop = &stop1;
    pstop2 = &stop2;
	InpIt stIt;
	double thruVol=0, depVol=0,dblsegDVij=0,dblT1=0,dblT2=0,dblT3=0,dblT4=0,dblT5=0,dblT6=0,dblT7=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
	double dblD1=0,dblD2=0,dblD3=0,dblD4=0,dblD5=0,dblarrDelay=0,dbldepDelay=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
	double dblregCumTmi=0,dblregCumTmj=0,dblregCumTmk=0,dblundCumTmi=0,dblundCumTmj=0,dblundCumTmk=0;
	double dblpredCumTm=0;
	double dblsegRegij=0,dblsegRegjk=0,dblsegUndij=0,dblsegUndjk=0, dblWkToTime=0;
	double dblsegRdTime=0, dblsegRunTime=0,dblWkCost=0, dblRdCost=0,dblOpCost=0;
	long stopi=0,stopj=0,stopk=0;
	double dblwkToTime=0;
	double dblHdwy = pdh.get_hdway(); // headway in minutes
	double dblBegTm = pdh.get_begtm();
	double dblPdLen = pdh.get_pdlen();
	double dblWkSpd = gc.get_walkspd()*3600/1000;
	double dblUnitWkCost = gc.get_walkcost();
	double dblUnitRdCost = gc.get_ridecost();
	double dblUnitOpCost = pdh.get_opercost();
	double dblunitOnTm=gc.get_unitontm();
	double dblunitOffTm=gc.get_unitofftm();
	double dblOffs=0, dblOns=0; 
	stIt = begin;

	while (stIt != end) {
		pstop = &(stIt->second);
		if (!pstop->get_blnExtr()) {
			pstop->set_probStoph(1);
			pstop->set_probStop(1);
		} else {
			if (blnHist) {
				pstop->cprobstoph(dblHdwy);
			} 
			pstop->cprobstop(dblHdwy);
		}
		if (stIt==begin) {
			// calculate preliminary departing volume
			if (pstop->get_blnHist() && pstop->get_blnExtr()) 
			{
				depVol = pstop->get_HistDepVol();

			} else
			{
				depVol = pstop->get_HistDepVol() - pstop->get_HistOns() +
							pstop->get_HistOffs() + pstop->get_Ons() - pstop->get_Offs();
			}
			pstop->set_DepVol(depVol);
			pstop->set_undCRdTm(pstop->get_CRdTm());
			pstop->set_CRdTmC(pstop->get_CRdTm());
			pstop->set_depDelay(0);
			pstop->set_arrDelay(0);
			stIt++;
		}
		else if (++stIt == end)
		{
			if (!pstop->get_blnExtr() && !pstop->get_blnHist()) {
				pstop->set_DepVol(pstop->get_HistDepVol());
				}
			else {
				pstop->set_DepVol(pstop0->get_DepVol()+pstop->get_Ons()-pstop->get_Offs());
			}
		}
		else 
		{
			pstop->set_DepVol(pstop0->get_DepVol()+pstop->get_Ons()-pstop->get_Offs());
			stIt++;
		}
			stop0 = *pstop;
	}
	stIt = begin;
	while(stIt != end) {
		pstop = &(stIt->second);
		if (begin==stIt) {
			// calculate preliminary departing volume
			depVol = pstop->get_DepVol();
            dblsegDVij = 0;
		// copy the current stop object into the previous stop object
			stop0 = *pstop;
		// get the next stop to form an ijk pair (i=j, k = j+1)
			stop2 = (++stIt)->second;
			--stIt;
		}
		else if ((++stIt) != end)
		{
			stop2 = stIt->second;
			(--stIt);
		}
		else
		{
			(--stIt);
			stop2 = stIt->second;
		}

		// get the next stop to form an ijk pair (i=j, k = j+1)
		dblOffs = pstop->get_HistOffs();
		dblOns = pstop->get_HistOns();
		dblarrDelay = pstop->get_arrDelay();
		dbldepDelay = pstop->get_depDelay();

	    dblregCumTmi=pstop0->get_CRdTm();
	    dblregCumTmj=pstop->get_CRdTm();
	    dblregCumTmk=pstop2->get_CRdTm();
	    dblundCumTmi=pstop0->get_undCRdTm();
	    dblundCumTmj=pstop->get_undCRdTm();
	    dblundCumTmk=pstop2->get_undCRdTm();
	    stopi=pstop0->get_id();
	    stopj=pstop->get_id();
	    stopk=pstop2->get_id();
		dblsegRegij = dblregCumTmj-dblregCumTmi;
		dblsegRegjk = dblregCumTmk-dblregCumTmj;
		dblsegUndij = dblundCumTmj-dblundCumTmi;
		dblsegUndjk = dblundCumTmk-dblundCumTmj;
		thruVol = depVol - pstop->get_HistOns();
		dblwkToTime = pstop->get_hWkTmOffs()+pstop->get_hWkTmOns();


// Alighting Time 
		dblT2 = dblunitOffTm * dblOffs * (dblHdwy/60); // Delay per stop Secs per hour
// Boarding Time 
        dblT3 = dblunitOnTm * dblOns * (dblHdwy/60); // Delay per stop Secs per hour
// Decelerating time  
		dblT4 = (dblarrDelay * pstop->get_probStop());  // Delay per stop Secs
// Accelerating Time 
        dblT5 = (pstop0->get_depDelay() * pstop0->get_probStop());  // Delay per stop Secs
// Decelerating delay experienced by Alighting passengers
        dblT6 = (dblarrDelay * dblOffs); // Delay per stop Secs per hour
// Accelerating delay experienced by Boarding passengers
        dblT7 = (pstop0->get_depDelay() * dblOns ); // Delay per stop Secs per hour
// Alighting delay experienced by through passengers
        dblD2 = dblT2*(dblHdwy/60) * (thruVol + (0.5 * dblOffs)); // Delay per stop Secs per hour 
// Boarding delay experienced by through passengers
        dblD3 = (dblT3*(dblHdwy/60) * thruVol); // Delay per stop Secs per hour
// Decelerating delay experienced by through passengers 
        dblD4 = (dblT4 * thruVol); // Delay per stop Secs secs per hour
// Accelerating delay experienced by through passengers
        dblD5 = (dblT5 * thruVol); // Delay per stop Secs per hour
//  
		dblT1 = (dblT2 + dblT3 + dblT4 + dblT5)/3600;  // stopping delay per stop (hrs) 
		dblD1 = (dblD2 + dblD3 + dblD4 + dblD5 + dblT6 + dblT7)/3600; // Riding Delay (Person-hrs/hr)     
		pstop->set_rideDelay(dblD1);      
//      wkCost = Cw1/Vwalk * ((wkToTime)/1000*60) * PdLen  
        dblWkCost = ((dblUnitWkCost * (dblwkToTime/60)) * dblPdLen);
       //Vr = 25 * 5280/3280.83333 

     
    // Segment Riding Cost
             dblsegRdTime = (dblD1* 60)  + (((dblsegUndjk*depVol)+(dblsegUndij*dblsegDVij))/(2));  //pers-min/hr 
    // Segment Operating Time
             dblsegRunTime = ((dblsegUndjk+dblsegUndij)/2) + (dblT1 * 60); // segment ride time min
 
//      if (dpW) then
//         rdCost = (SegRdTime/60 * PdLen)
    // Incremental Operating Cost
//         OpCost = ((SegRunTime/60)  * (60 / PHdway) * PdLen)
//      else
//        dblRdCost = dblUnitRdCost * dblPdLen * (dblD1  + (dblsegUndjk * depVol +dblsegUndij * dblsegDVij)/(2*60)) ;
   // Incremental Operating Cost
        dblOpCost = (dblUnitOpCost * (dblsegRunTime/60)  * (60 / dblHdwy) * dblPdLen);
        dblRdCost = dblUnitRdCost * dblsegRdTime / 60 * dblPdLen;
/*
		if ((pstop->get_blnExtr()) && (pstop->get_blnHist())) {
	         dblWkCost = 0;
	         dblRdCost = 0;
	         dblOpCost = 0;
	         dblwkToTime = 0;
	         dblsegRdTime = 0;
	         dblsegRunTime = 0;
		}
 */
		pstop->set_WalkCost(dblwkToTime);
         pstop->set_RideCost(dblsegRdTime); 
         pstop->set_OperCost(dblsegRunTime);
         double dbltotCost =  (dblWkCost+dblRdCost+dblOpCost); 
         pstop->set_TCost(dbltotCost);
// update the predicted cumulative time table
		
		 if (stIt!=begin) {
            dblpredCumTm = dblpredCumTm + dblsegUndij +
			(pstop0->get_probStop()*pstop0->get_depDelay()+dblT2+dblT3+pstop->get_probStop()*dblarrDelay)/60;   
		 } else {
            dblpredCumTm = dblregCumTmj;          
		 }
         pstop->set_CRdTmC(dblpredCumTm);
// Pstopi - probability of stopping for the prior stop (set =1 for the first stop)
           dblsegUndij = dblsegUndjk;
           dblsegDVij = depVol;
           dblsegRegij = dblsegRegjk;

      // loop over stops
		stop0 = *pstop;
			(++stIt);
	}
}
#endif // Walk, Ride, Operating cost IMPACTCALC_H ///:~

#ifndef IMPACTCALCMAP_H //Walk, Ride, Operating cost
#define IMPACTCALCMAP_H
#include <algorithm>

// Map Container version:

template<typename m, typename g, typename h,typename q, class InpIt> 
m& ImpactCalcMap(m& m1, g& gc, h& pdh,q& q1) {
	//initialize the stop variables needed for the program
    tstop stop0;
	tstop stop1;
	tstop stop2;
    tstop* pstop0= &stop0;
	tstop* pstop= &stop1;
	tstop* pstop2= &stop2;
	InpIt stIt;
	InpIt begin= m1.begin();
	InpIt end = m1.end();
	double thruVol=0, depVol=0,dblsegDVij=0,dblT1=0,dblT2=0,dblT3=0,dblT4=0,dblT5=0,dblT6=0,dblT7=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
	double dblD1=0,dblD2=0,dblD3=0,dblD4=0,dblD5=0,dblarrDelay=0,dbldepDelay=0; // delay components t2 - alighting delay, t3 - boarding delay, t4 - arrival , t5 - departure delay
	double dblregCumTmi=0,dblregCumTmj=0,dblregCumTmk=0,dblundCumTmi=0,dblundCumTmj=0,dblundCumTmk=0;
	double dblpredCumTm=0;
	double dblsegRegij=0,dblsegRegjk=0,dblsegUndij=0,dblsegUndjk=0, dblWkToTime=0;
	double dblsegRdTime=0, dblsegRunTime=0,dblWkCost=0, dblRdCost=0,dblOpCost=0;
	long stopi=0,stopj=0,stopk=0;
	double dblwkToTime=0;
	double dblHdwy = pdh.get_hdway(); // headway in minutes
	double dblBegTm = pdh.get_begtm();
	double dblPdLen = pdh.get_pdlen();
	double dblWkSpd = gc.get_walkspd()*(gc.get_unitconv()); // in ft/sec
	double dblUnitWkCost = gc.get_walkcost();
	double dblUnitRdCost = gc.get_ridecost();
	double dblUnitOpCost = pdh.get_opercost();
	double dblunitOnTm=gc.get_unitontm();
	double dblunitOffTm=gc.get_unitofftm();
	double dblOffs=0, dblOns=0; 
	stIt = m1.begin();

	while (stIt != m1.end()) {
		pstop = &(stIt->second);
		pstop->cprobstoph(dblHdwy);
		pstop->cprobstop(dblHdwy);
		if (stIt==begin) {
			// calculate preliminary departing volume
			depVol = pstop->get_HistDepVol() - pstop->get_HistOns() +
							pstop->get_HistOffs() + pstop->get_Ons() - pstop->get_Offs();
//			}
			pstop->set_DepVol(depVol);
			pstop->set_undCRdTm(pstop->get_CRdTm());
			pstop->set_CRdTmC(pstop->get_CRdTm());
//			pstop->set_depDelay(0);
			pstop->set_arrDelay(0);
			if (!pstop->get_blnExtr()) {
				pstop->set_probStop(1);
			}
		}
		else if (stIt != end)
		{
			depVol = pstop0->get_DepVol()+pstop->get_Ons()-pstop->get_Offs();
			pstop->set_DepVol(depVol);
		} 
			stop0 = *pstop;
			stIt++;
	} // while compute Dep. Vol. 
		if (stIt-- != end)
		{stop1 = (stIt->second);}
		if (stIt-- != end)
		{stop0 = (stIt->second);}

	if (!pstop->get_blnExtr() && !pstop->get_blnHist()) {
				pstop->set_DepVol(pstop->get_HistDepVol());
		}
		else if (pstop0->get_id()==pstop->get_id()) {
			pstop->set_DepVol(pstop->get_HistDepVol());
		} else {
			pstop->set_DepVol(pstop0->get_DepVol()+pstop->get_Ons()-pstop->get_Offs());
		}


	stIt = m1.begin();
	while(stIt != m1.end()) {
		pstop = &(stIt->second);
			depVol = pstop->get_DepVol();
		if (stIt==m1.begin()) {
			// calculate preliminary departing volume
            dblsegDVij = 0;
		// copy the current stop object into the previous stop object
			stop0 = *pstop;
		// get the next stop to form an ijk pair (i=j, k = j+1)
			stop2 = (++stIt)->second;
			--stIt;
		} else if ((++stIt) != m1.end()) // k is not the last stop
		{
			stop2 = stIt->second;
			(--stIt);
		} else { // k is the last stop in the set
			(--stIt);
			stop2 = stIt->second;
			stop2.set_probStop(1.0);
			stop2.set_depDelay(0.0);
		}

		// get the next stop to form an ijk pair (i=j, k = j+1)
		dblOffs = pstop->get_Offs();
		dblOns = pstop->get_Ons();
		dblarrDelay = pstop->get_arrDelay();
		dbldepDelay = pstop->get_depDelay();

	    dblregCumTmi=pstop0->get_CRdTm();
	    dblregCumTmj=pstop->get_CRdTm();
	    dblregCumTmk=pstop2->get_CRdTm();
	    dblundCumTmi=pstop0->get_undCRdTm();
	    dblundCumTmj=pstop->get_undCRdTm();
	    dblundCumTmk=pstop2->get_undCRdTm();
	    stopi=pstop0->get_id();
	    stopj=pstop->get_id();
	    stopk=pstop2->get_id();
		pstop->set_Stopidp(stopi);
		pstop->set_Stopids(stopk);

		dblsegRegij = dblregCumTmj-dblregCumTmi;
		dblsegRegjk = dblregCumTmk-dblregCumTmj;
		dblsegUndij = dblundCumTmj-dblundCumTmi;
		dblsegUndjk = dblundCumTmk-dblundCumTmj;
//		if ( pstop0->get_DepVol() <=0 && q1 >0 && q1 <9999)  // removed becuase ride delays became very high for edge segments and bias costs
//		{
//			thruVol = inf;
//		} else {
			thruVol = pstop0->get_DepVol() - pstop->get_Ons();
//		}
		dblwkToTime = pstop->get_WkTmOffs()+pstop->get_WkTmOns();


// Alighting Time 
		dblT2 = dblunitOffTm * dblOffs * (dblHdwy/60); // Delay per stop Secs per hour per bus
// Boarding Time 
        dblT3 = dblunitOnTm * dblOns * (dblHdwy/60); // Delay per stop Secs per hour per bus
// Decelerating time  
		dblT4 = (dblarrDelay * pstop->get_probStop());  // Delay per stop Secs
// Accelerating Time 
        dblT5 = (pstop0->get_depDelay() * pstop0->get_probStop());  // Delay per stop Secs
// Decelerating delay experienced by Alighting passengers
        dblT6 = (dblarrDelay * dblOffs); // Delay per pax-stop Secs per hour
// Accelerating delay experienced by Boarding passengers
        dblT7 = (pstop0->get_depDelay() * dblOns ); // Delay per pax-sec per hour
// Alighting delay experienced by through passengers
        dblD2 = dblT2 * (thruVol + (0.5 * dblOffs)); // Delay per Pax-sec per hour 
// Boarding delay experienced by through passengers
        dblD3 = (dblT3 * thruVol); // Delay per stop Pax-Secs per hour
// Decelerating delay experienced by through passengers Pax-Secs per hour
        dblD4 = (dblT4 * thruVol); // Delay per stop Pax-Secs per hour
// Accelerating delay experienced by through passengers
        dblD5 = (dblT5 * thruVol); // Delay per stop Pax-Secs per hour
//  
		dblT1 = (dblT2 + dblT3 + dblT4 + dblT5)/3600;  // stopping delay per stop (hrs) 
		dblD1 = (dblD2 + dblD3 + dblD4 + dblD5 + dblT6 + dblT7)/3600; // Riding Delay (Pax-hrs/hr)     
		pstop->set_rideDelay(dblD1);      
//      wkCost = Cw1/Vwalk * ((wkToTime)/1000*60) * PdLen  
        dblWkCost = ((dblUnitWkCost * (dblwkToTime/60)) * dblPdLen);
		
       //Vr = 25 * 5280/3280.83333 
     
    // Segment Riding Cost
             dblsegRdTime = (dblD1* 60)  + (((dblsegUndjk*depVol)+(dblsegUndij*dblsegDVij))/(2));  //pers-min/hr 
    // Segment Operating Time
             dblsegRunTime = ((dblsegUndjk+dblsegUndij)/2) + (dblT1 * 60); // segment ride time min
 
//      if (dpW) then
//         rdCost = (SegRdTime/60 * PdLen)
    // Incremental Operating Cost
//         OpCost = ((SegRunTime/60)  * (60 / PHdway) * PdLen)
//      else
//        dblRdCost = dblUnitRdCost * dblPdLen * (dblD1  + (dblsegUndjk * depVol +dblsegUndij * dblsegDVij)/(2*60)) ;
   // Ride Cost
        dblRdCost = dblUnitRdCost * dblsegRdTime / 60 * dblPdLen;
   // Incremental Operating Cost
        dblOpCost = (dblUnitOpCost * (dblsegRunTime/60)  * (60 / dblHdwy) * dblPdLen);

		pstop->set_WalkCost(dblWkCost);
         pstop->set_RideCost(dblRdCost); 
         pstop->set_OperCost(dblOpCost);
         double dbltotCost =  (dblWkCost+dblRdCost+dblOpCost); 
         pstop->set_TCost(dbltotCost);
// update the predicted cumulative time table
		
		 if (stIt!=m1.begin()) {
            dblpredCumTm = dblpredCumTm + dblsegUndij +
			(pstop0->get_probStop()*pstop0->get_depDelay()+dblT2+dblT3+pstop->get_probStop()*dblarrDelay)/60;   
		 } else {
            dblpredCumTm = dblregCumTmj;          
		 }
         pstop->set_CRdTmC(dblpredCumTm);
// Pstopi - probability of stopping for the prior stop (set =1 for the first stop)
           dblsegUndij = dblsegUndjk;
           dblsegDVij = depVol;
           dblsegRegij = dblsegRegjk;

      // loop over stops
		stop0 = *pstop;
			(++stIt);
	}
	return m1;
}

#endif // Walk, Ride, Operating cost IMPACTCALC_H ///:~


// Calculate undelayed time in an STL sequence container.
#ifndef MAXMINRTM_H
#define MAXMINRTM_H
#include <algorithm>
// find the minimum run time from the stop list
template<class Seq> double minRTm(Seq& c) {
  typename Seq::iterator psit,psit2;
  tstop* pstop;
  psit2 = c.begin();
  pstop = &(psit2->second);
  double rtm1 = pstop->get_CRdTm();
  for(psit = c.begin(); psit != c.end(); ++psit) {
		pstop = &(psit->second);
		if (rtm1 > pstop->get_CRdTm()) {
			rtm1 = pstop->get_CRdTm();
		}
  }
  return rtm1;
}
// find the maximum run time from the stop list
template<class Seq> double maxRTm(Seq& c) {
  typename Seq::iterator psit,psit2;
  tstop* pstop;
  psit2 = c.begin();
  pstop = &(psit2->second);
  double rtm1 = pstop->get_CRdTm();
  for(psit = c.begin(); psit != c.end(); ++psit) {
		pstop = &(psit->second);
		if (rtm1 < pstop->get_CRdTm()) {
			rtm1 = pstop->get_CRdTm();
		}
  }
  return rtm1;
}
// find the run time from the boarding stop to the end
template<class Seq> void RTmE(Seq& c, double RTm1) {
  typename Seq::iterator psit;
  tstop* pstop;
  for(psit = c.begin(); psit != c.end(); ++psit) {
		pstop = &(psit->second);
		pstop->calCRdTmE(RTm1);
  }
}

// parcel list processing
template <typename s>
s& rideTime2End(s& tsrtmap)
{
	// mmaparStop - Stop Assigned -> Parcel Object multi-map 
    long ip;
	double dblRt=0; 
	tstop* pstop;
 tsrtmit = tsrtmap.begin();
 while (tsrtmit !=tsrtmap.end()) {
	 pstop = &(tsrtmit->second);
	 ip = pstop->get_id();
	 if (dblRt < pstop->get_CRdTm()) {
		 dblRt=pstop->get_CRdTm();
	 }
	tsrtmit++;
 }
 tsrtmit = tsrtmap.begin();
 while (tsrtmit !=tsrtmap.end()) {
	 pstop = &(tsrtmit->second);
	 ip = pstop->get_id();
	 pstop->set_CRdTmE(dblRt-pstop->get_CRdTm()); 
	tsrtmit++;
 }
return tsrtmap;
}

#endif // MAXMINRTM_H ///:~

#ifndef SINGLETONREPORT_H
#define SINGLETONREPORT_H
// A singleton. Will automatically report the
// statistics as the program terminates:
class gcostReport {
  static gcostReport gr;
  gcostReport() {} // Private constructor
public:
  ~gcostReport() {
    std::cout << "\n-------------------\n"
 /*     << "globcost creations: " << globcost::create
      << "\nCopy-Constructions: " 
      << globcost::copycons
      << "\nAssignments: " << globcost::assign
      << "\nDestructions: " << globcost::destroy */
      << std::endl;
  }
};

#endif  // SINGLETONREPORT_H

#ifndef SYNCHSTOPVERTEXMAP_H
#define SYNCHSTOPVERTEXMAP_H

// map container sorter
template <typename v , typename s, typename d>
d& dpSynchStopVertexMap(v& v1, s& s1,d& d1)
{  
double x,y;
v::iterator vit;
s::iterator si1;
long sid,vid;
  for(vit=v1.begin();vit!=v1.end();vit++)
  {
	  vid = vit->second.get_id();
	  sid = vit->second.get_orig();
	  if (sid>0) {
	    x = vit->second.get_cost();
		si1= s1.find(sid);
		y = si1->second.get_CRdTm();
		if (x>=inf) {
			vit->second.set_cost(y);
			vit->second.set_tcost(y);
			d1.insert(dblng_Pair(y,vit->second.get_id()));
		}
	  }
  }
  return d1;
}

// map Stop - Vertex container cleanup
template <typename u, typename v, typename w, typename x, typename y>
u& dpSynchLngVertexMap(u& u1, v& v1, w& w1, x& x1, y& y1)
{  
u::iterator uit;
v::iterator vit;
typedef pair <x,y> obj_pair;

	for(uit=u1.begin();uit!=u1.end();uit++)
	{
		x1 = uit->first;
		y1 = uit->second;
		vit= v1.find(y1);
		if (vit!=v1.end()) 
		{
				w1.insert(obj_pair(x1,y1));
		}
	}
  return u1;
}


// routine to remove stop related vertices from a multi-map using double values
// remove the current Vertex from the origin map before doing vertex voronoi routine 
// mmapStopVx - u - relation between stop and vertex ids'
// mmapVx0R - w - vertex Id's and related cost (reverse of mmapVx0)  
// if stop lng id is not found in alternative stop map 
// the mapping goes from lng-lng -> lng-dbl -> dbl-lng
template <typename u, typename v, typename w, typename d, typename h>
d& dpSynchDblVertexMap(u& u1, v& v1, w& w1, d& d1, h& h1)
{  
double x=0;
long sid=0,vid=0;
u::iterator uit;
v::iterator vit;
w::iterator wit;
d::iterator dit;

	for(uit=u1.begin();uit!=u1.end();uit++)
	{
		sid = uit->first;
		vit= v1.find(sid);
		// if the stop is not present (changed the external removal on  4/8/2009 then delete vertices associated with it
		if (vit==v1.end()) { // || (vit->second.get_blnExtr() && (!h1)) ) {
		    vid = uit->second;
			wit= w1.find(vid);
			if (wit!=w1.end())
			{
				x = wit->second;
				dit= d1.find(x);
				if (dit!=d1.end())
				{
					d1.erase(dit);
				} 
			} 
		} 
	}
  return d1;
}

#endif // SYNCHSTOPVERTEXMAP_H

// make a new map by changing the trip key to be the sorting key
template <typename m, typename n, typename o, typename p, typename t >
n& remapdHdwy2pdKey( m& m1, n& n1, o& o1,p& p1, t& triPd)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
		o1 = mit->second;
		p1 = o1.get_pdKey();
		n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}



#ifndef SORTPARCELONOFFID_H
#define SORTPARCELONOFFID_H

// map container update vertex 
template <typename v , typename s, typename o>
o& updateVert(v& v1, s& s1,o& o1)
{
  
double x,y;
v::iterator vit;
s::iterator si1;
long sid,vid;
  for(vit=v1.begin();vit!=v1.end();vit++)
  {
	  vid = vit->second.get_id();
	  sid = vit->second.get_orig();
	  if (sid>0) {
	    x = vit->second.get_cost();
		si1= s1.find(sid);
		y = si1->second.get_CRdTm();
		if (x>=inf) {
			vit->second.set_cost(y);
			vit->second.set_tcost(y);
			o1.insert(dblng_Pair(y,vit->second.get_id()));
		}
	  }
  }
  return o1;
}

// Sort the parcel data into a multi-map with either on or off id associated with parcel object
// p- parcel multimap in any order, sp1 parcel multimap in stop off -> onoff=0 or on -> onoff=1 order
template <typename p, typename sp, typename o >
sp& SortOnOffIdParcObject(p& p1, sp& sp1, o& onoff)
{
	p::iterator pit;
	long stid=0;
	parcel parc;
	for(pit=p1.begin();pit!=p1.end();pit++)
	{
		parc = pit->second;
		if (onoff) {
			stid = parc.get_origon();
		} else {
			stid = parc.get_origoff();
		}
		if (stid>0) {
			sp1.insert(PE_Pair(stid,parc));
		}
	} 
return sp1;
}

#endif // SORTPARCELONOFFID_H ///:~


#ifndef DPDATACREATION_H
#define DPDATACREATION_H
// map DP result with Walk and Ride Costs 
template <typename u, typename v, typename w, typename x, typename y, typename t>
v& dpTripWalkRideMap(u& u1, v& v1, w& w1, x& x1, y& y1, t& t1 )
{  
u::iterator uit;  // DP resutl Map data
v::iterator vit;  // trip dp, walk / ride map
typedef pair <w,y> obj_pair;  // x = tripID, y = tripDP
y* py; 
py=&y1;
	for(uit=u1.begin();uit!=u1.end();uit++)
	{
		x1 = uit->second;
		py->set_dptStop(x1);
		py->set_tripId( w1);
		py->set_walk((int) (t1->get_walkcost()*10));
		py->set_ride((int) (t1->get_ridecost()*10));
		v1.insert(obj_pair(w1,*py));
	}
  return v1;
}

#endif // DP_Result_map_creation ///:~


#ifndef READRECORDATA_H
#define  READRECORDATA_H

// Read input record into object data in their respective object


template <typename o, typename m >
o& readgcost1(string& rec1, o& o1, m& maphdrit, char *seps ) // *seps = "\t" 
{
   static int rno = 0; //record number

   int i=0; int ib=10;

//COSTWALK,COSTRIDE,UNITONTM,UNITOFFTM,MAXWLKDIST,PROPENSITY,FILESTEM,NOPERIODS,
//WALKSPD,FILEPATH
			 double walkcost;  // vertex id
	         short nopds; // no of periods
	         short maxskip; // Maximum no of skips allowed
         	 float unitontm;  // unit boarding time
             float unitofftm;  // unit alighting time
			 float maxwalkdist; // total cost
			 double ridecost;  // ride cost 
			 float propensity;  // propensity ratio
	         float walkspd; // walk speed 
	         string  filestem; // stem for file naming 
	         string  filepath; // path name for output files 

   string f1, fldval,fldhdr;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
 
   rno++;

  mapflds1 = recoread(rec1);
  mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;

	fldhdr = maphdrit->second;


     std::transform(fldhdr.begin(), fldhdr.end(), fldhdr.begin(), to_lower());
	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "COSTWALK")   // Walk Cost
	{
		walkcost = fromString<double>(fldval);//val1,&stop1);
	  if (walkcost>0) {
		  o1.set_walkcost(walkcost);
	  }
	  else
	  {
          o1.set_walkcost(0);	 
          o1.set_nopds(0);	 
		  return o1;
	  }
	}
	else if (f1 == "COSTRIDE") // Ride Cost
	{
      ridecost = fromString<double>(fldval); //strtod(val1,&stop1);
	  if (walkcost>0) {
		  o1.set_ridecost(ridecost);
	  }
	}
	else if (f1 == "UNITONTM") // Unit On Time
	{
		unitontm = fromString<float>(fldval);
		  o1.set_unitontm(unitontm);
	}
	else if (f1 == "UNITOFFTM") // Unit Off Time
	{
		unitofftm = fromString<float>(fldval); //strtod(val1,&stop1);
		  o1.set_unitofftm(unitofftm);
	}
	else if (f1 == "MAXWLKDIST") // Max Walk Distance
	{
			maxwalkdist = fromString<float>(fldval); //strtod(val1,&stop1);
			o1.set_maxwalkdist(maxwalkdist);
	}
	else if (f1 == "PROPENSITY") // propensity
	{
		propensity = fromString<float>(fldval); //strtod(val1,&stop1);
		o1.set_propensity(propensity); 
	}
	else if (f1 == "FILESTEM") // FILE STEM used for naming output files
	{
		    o1.set_filestem(filestem);
	}
	else if (f1 == "NOPERIODS" || f1 == "PERIOD") 
		// Period of this run used to extract the headway and length of period from the period headway file
	{
		nopds = fromString<short>(fldval);  //strtol(val1,&stop1,ib);
			o1.set_nopds(nopds); 
	}
	else if (f1 == "WALKSPD") // WALKing SPeeD
	{
			walkspd = fromString<float>(fldval); //strtod(val1,&stop1);
		    o1.set_walkspd(walkspd);
	}
	else if (f1 == "FILEPATH") // FILE path for saving output files
	{
		    o1.set_filepath(filepath);
	}
	else if (f1 == "MAXSKIP") // No of Periods included in the analysis
	{
		maxskip = fromString<short>(fldval); 
			o1.set_maxskip(maxskip); 
	}
	else if (f1 == "DPDIMENSION") // No of Periods included in the analysis
	{
		dpdimension = fromString<int>(fldval); 
			o1.set_dpdimension(dpdimension); 
	}
   // get the next fld
	mapfldsit++;
	maphdrit++;
}
   if (o1.get_walkcost()>0)
   {
	   o1.show_globcosthdr(cout);
	   o1.show_globcost(cout);
   }

   return o1;
}

template <typename o, typename m >
o& readpdhdway1(string& rec1, o& o1, m& maphdrit, char *seps ) // *seps = "\t" 
{
   static int rno = 0; //record number

	int i=0; int ib=10; int j=0;

//FLDPERD,FLDBEGT,FLDHDWAY,FLDPDLEN,COSTOPER
	         short pdId; // period Id
			 float begtm; // Beginning of Time Period 
			 float hdway;  // headway 
			 float pdlen;  // period length 
	         float opercost; // Operating Cost 

   const int bufsz = 1000; // Buffer size;

   string f1, fldval, fldhdr;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
   rno++;

  mapflds1 = recoread(rec1);
 
   mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;
	fldhdr = maphdrit->second;
	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "FLDPERD")  // period id
	{
		pdId = fromString<short>(fldval); //strtol(val1,&stop1,ib);
	  if (pdId>0) {
		  o1.set_pdId(pdId);
	  }
	  else
	  {
		  o1.set_pdId(0);
		  return o1;
	  }
	}
	if (f1 == "FLDBEGT") // Period Begining Time
	{
		begtm = fromString<float>(fldval); //strtod(val1,&stop1);
		  o1.set_begtm(begtm);
	  }
	if (f1 == "FLDHDWAY") // headway
	{
		hdway = fromString<float>(fldval); //strtod(val1,&stop1);
		  o1.set_hdway(hdway);
	  }
	if (f1 == "FLDPDLEN") // period length
	{
		pdlen = fromString<float>(fldval); //strtod(val1,&stop1);
		  o1.set_pdlen(pdlen);
	  }
	if (f1 == "COSTOPER") // operating cost
	{
		opercost = fromString<float>(fldval); //strtod(val1,&stop1);
		o1.set_opercost(opercost);
	} 
   // get the next fld
	mapfldsit++;
	maphdrit++;
}

   return o1;
}

template <typename r, typename o >
o& readInputFiles(r& rec1, o& o1, char *seps ) // *seps = "\t" 
{

// routename global cost table , period headway table, landuse table, stop table, edge table, vertex table , parcel table , scenario table
			 
   const int bufsz = 1000; // Buffer size;

   string fldval;
   typedef pair < string, string >  fld_pair;
	int pos;
	pos = rec1.find(seps);    // position of seps in str
	fldval = rec1.substr (0,pos);   // get from seps to the end
	fldval = (stringUpper<string>(fldval)); 
	if (fldval == "GCOSTTABLE" || fldval == "GLOBALCOSTTABLE") {
		fldval = rec1.substr (pos+1);   // the name of the Global Cost Table
		o1.set_tblgcost(fldval);
	} else if (fldval == "ROUTENAME") {
		fldval = rec1.substr (pos+1);   // get the name of the Transit Route to be run
		o1.set_rtename(fldval);
	} else if (fldval == "DIRNAME") {
		fldval = rec1.substr (pos+1);   // get the Route Direction
		o1.set_dirname(fldval);
	} else if (fldval == "PERIODTABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Period Table
		o1.set_tblperiod(fldval);
	} else if (fldval == "LANDUSETABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Land Use Table
//		o1.insert(fldpair("LANDUSETABLE",fldval));
		o1.set_tbllanduse(fldval);
	} else if (fldval == "TRIPTABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Trip Table
		o1.set_tbltrip(fldval);
	} else if (fldval == "STOPTABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Stop Table
		o1.set_tblstop(fldval);
	} else if (fldval == "EDGETABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Edge Table
		o1.set_tbledge(fldval);
	} else if (fldval == "VERTEXTABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Vertex Table
		o1.set_tblvertex(fldval);
	} else if (fldval == "PARCELTABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Parcel Table
		o1.set_tblparcel(fldval);
	} else if (fldval == "SCENETABLE") {
		fldval = rec1.substr (pos+1);   // get the name of the Scene Table
		o1.set_tblscene(fldval);
	}  else if (fldval == "SRID" ) {
		fldval = rec1.substr (pos+1);   // get the Level of Traffic Stress 
		o1.set_srid(from_string<int>(fldval));
	}  else if (fldval == "LTS" || fldval == "RTSTRESS") {
		fldval = rec1.substr (pos+1);   // get the Level of Traffic Stress 
		o1.set_lts(from_string<int>(fldval));
	} else if (fldval == "RUNHISTORIC" || fldval == "HISTORIC") {
		fldval = stringUpper<string>(rec1.substr (pos+1));   // get if Historic Run is to be run first
		if (fldval == "FALSE" || fldval == "F" || fldval == "0" )
			{
				o1.set_historic(false);
			} else {
				o1.set_historic(true);
			}
	} else if (fldval == "EUCLID" || fldval == "EUCLIDEAN") {
		fldval = stringUpper<string>(rec1.substr (pos+1));   // get if Euclidean Distance is to be used
		if (fldval == "TRUE" || fldval == "T" || fldval == "1" )
			{
				o1.set_euclid(true);
			} else {
				o1.set_euclid(false);
			}
	} else if (fldval == "PERIODRUN" || fldval == "PDRUN") {
		fldval = stringUpper<string>(rec1.substr (pos+1));   // get if Euclidean Distance is to be used
		if (fldval == "TRUE" || fldval == "T" || fldval == "1" )
			{
				o1.set_pdRun(true);
			} else {
				o1.set_pdRun(false);
			}
	} else if (fldval == "PARCEL"  || fldval == "PARCELEVEL") {
		fldval = stringUpper<string>(rec1.substr (pos+1));   // check if parcel data 
		if (fldval == "TRUE" || fldval == "T" || fldval == "1" )
			{
				o1.set_parcel(true);
			} else {
				o1.set_parcel(false);
			}
	} else if (fldval == "INPUTDBASE" || fldval == "INPUTDB") {
		fldval = rec1.substr (pos+1);   // get input database file 
		o1.set_inputDB(fldval);
	} else if (fldval == "RESULTDBASE" || fldval == "RESULTDB") {
		fldval = rec1.substr (pos+1);   // get input database file 
		o1.set_resultDB(fldval);
	} else if (fldval == "PERIODFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Period Table
		o1.set_fileperiod(fldval);
	} else if (fldval == "LANDUSEFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Land Use Table
//		o1.insert(fldpair("LANDUSETABLE",fldval));
		o1.set_filelanduse(fldval);
	} else if (fldval == "TRIPFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Trip Table
		o1.set_filetrip(fldval);
	} else if (fldval == "STOPFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Stop Table
		o1.set_filestop(fldval);
	} else if (fldval == "EDGEFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Edge Table
		o1.set_fileedge(fldval);
	} else if (fldval == "VERTEXFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Vertex Table
		o1.set_filevertex(fldval);
	} else if (fldval == "PARCELFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Parcel Table
		o1.set_fileparcel(fldval);
	} else if (fldval == "SCENEFILE") {
		fldval = rec1.substr (pos+1);   // get the name of the Scene Table
		o1.set_filescene(fldval);
	}
   return o1;
}

template <typename o, typename m >
o& readlucodes1(string& rec1, o& o1, m& maphdrit, char *seps ) // *seps = "\t" 
{
   static int rno = 0; //record number

	int i=0; int ib=10;

//PTYPE,DESC_,LU,KEYFLD,ONCOEF,OFFCOEF
			 string pType;  // Property type code (DOR Codes)
			 string Desc;  // Description of Code  
			 string LUC; // Text Land Use Code (general)
         	 string KeyProp;  // Key Field/property from parcel object 
             float OnCoeff;  // Property Code production per Key Field
			 float OffCoeff; // Property Code Attraction per Key Field

//   const int bufsz = 1000; // Buffer size;
//   const char *val1;

   string f1, fldval,fldhdr;
   int ci=0;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
   rno++;
  mapflds1 = recoread(rec1);

  mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;

	fldhdr = maphdrit->second;
	f1 = (stringUpper<string>(fldhdr));
     std::transform(fldhdr.begin(), fldhdr.end(), fldhdr.begin(), to_lower());

	if (f1 == "PTYPE" || f1 == "SICTYPE") // Property Type
	{  
	  //ci = fromString<long>(fldval); // strtol(val1,&stop1,ib);
	  //if (ci>0) {
		  o1.set_pType(fldval);
	  //}
	  //else
	  //{
		//  o1.set_pType("0");
		//  return o1;
	  //}
	}
	else if (f1 == "DESC") // Description
	{  Desc = fldval;
		  o1.set_Desc(Desc);
	}
	else if (f1 == "LU" || f1 == "LUCODE" || f1 == "LANDUSECODE" ) // Land Use Code
	{
		LUC = fldval;
		o1.set_LUC(LUC);
	  }
	else if (f1 == "KEYFLD" )   // Property Type Text Code
	{
		KeyProp = fldval;
		o1.set_KeyProp(KeyProp);
	}
	else if (f1 == "ONCOEF") // Unit On rate per key field
	{  
		OnCoeff = fromString<float>(fldval); //strtod(val1,&stop1);
		o1.set_OnCoeff(OnCoeff);
	}
	else if (f1 == "OFFCOEF") // Unit Off rate per key field
	{
		OffCoeff = fromString<float>(fldval); //strtod(val1,&stop1);
		o1.set_OffCoeff(OffCoeff);
	} 
   // get the next fld
	mapfldsit++;
	maphdrit++;
}
   if (fromString<long>(o1.get_pType())>0)
   {
	   //o1.show_lucodeshdr(cout);
	   //o1.show_lucodes(cout);
   }

   return o1;
}

template <typename o, typename m >
 o& readpts1(string& rec1, o& o1, m& maphdrit , char *seps ) // *seps = "\t" 
{
   static int rno = 0; //record number

	int i1 = 0, ib=10, j1=0;
//Key,StopNamei,Stopi,StopNamej,Stopj,StopNamek,Stopk,StopNamel,Stopl,StopNamem,Stopm,CUM_TIME,ONS,OFFS,WalkCost,RideCost,OperCost,TCost

	long i;
	long j;
	long k;
	long l;
	long m;
	string StopName; // Stop Name
	double CRdTm; // supplied Cumulative ride time
	double Ons; 
	double Offs; 
	double WalkCost;
	double RideCost;
	double OperCost;
	double TCost;


  // const char *val1;
   string  f1, fldval,fldhdr;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
   rno++;
  mapflds1 = recoread(rec1);

 mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i1 = mapfldsit->first;
	fldval = mapfldsit->second;
//	val1 = fldval.c_str();
	fldhdr = maphdrit->second;

	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "STOPI") // tstop id
	{  i = fromString<long>(fldval); 
	  if (i>0) {
		  o1.set_i(i);
	  }
	  else
	  {
		  return o1;
	  }
	}
	else if (f1 == "STOPJ") // J
	{
		  j = fromString<long>(fldval); //strtol(val1,&stop1,ib);
		  o1.set_j(j);
	}
	else if (f1 == "STOPK") // k
	{
		  k = fromString<long>(fldval); //strtol(val1,&stop1,ib);
		  o1.set_k(k);
	}
	else if (f1 == "STOPL") // l
	{
		  l = fromString<long>(fldval); //strtol(val1,&stop1,ib);
		  o1.set_l(l);
	}
	else if (f1 == "STOPM") // l
	{
		  m = fromString<long>(fldval); //strtol(val1,&stop1,ib);
		  o1.set_m(m);
	}
	else if (f1 == "KEY") // k
	{
		StopName = fldval;
	}
	else if (f1 == "CUM_TIME") // cum. ride time
	{ 
		CRdTm = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_CRdTm(CRdTm);
	}
	else if (f1 == "ONS") // Ons
	{
		Ons = fromString<double>(fldval); 
		o1.set_Ons(Ons);} 
	else if (f1 == "OFFS") // Offs
	{	
		Offs = fromString<double>(fldval);
		o1.set_Offs(Offs);
	} 
	else if (f1 == "WALKCOST") // Walk Cost
	{	
		WalkCost = fromString<double>(fldval); 
		o1.set_WalkCost(WalkCost);
	}
	else if (f1 == "RIDECOST") // Ride Cost
	{	
		RideCost = fromString<double>(fldval); 
		o1.set_RideCost(RideCost);
	}
	else if (f1 == "OPERCOST") // Operating Cost
	{	
		OperCost = fromString<double>(fldval); 
		o1.set_OperCost(OperCost);
	}
	else if (f1 == "TCOST") // Total Cost
	{	
		TCost = fromString<double>(fldval); 
		o1.set_TCost(TCost);
	}

   // get the next fld
	mapfldsit++;
	maphdrit++;
}
   if (o1.get_i()>0)
   {
	   o1.set_dpkey(o1.makey(o1.get_i(),o1.get_j(),o1.get_k(),o1.get_l(),o1.get_m()));
	cout<<"i "<<o1.get_i()<<" j "<<o1.get_j()<<" k "<<o1.get_k()<<" l "<<o1.get_l()<<" m "<<o1.get_m()<<", Key : "<<o1.get_dpkey()<<
		", Cum Rd Tm "<<o1.get_CRdTm()<<", Ons "<<o1.get_Ons()<<", Offs "<<o1.get_Offs()<<endl;
   }

   return o1;
}

template <typename o, typename m >
o& readtstop1(string& rec1, o& o1, m& maphdrit, char *seps ) // *seps = "\t" 
{
   static int rno = 0; //record number

	int i=0; int ib=10; int j=0;
//OBJECTID,ORDER,LABEL,STOPNAME,CUM_DIST,CUM_TIME,UNDCRTM,CUM_TIME_P,HISTORIC,INBOUND,HISTONS,HISTOFFS,	
//HISTDEPVOL,ONS,OFFS,DEPVOL,PROBSTOP,DEPDELAY,ARRDELAY,DELAY,PVAL,AVAL,INCLUDE,CUM_TIME2	
//WKTMONS,WKTMOFFS,EXTERNAL,WALKCOST,RIDECOST,OPERCOST,TCOST,EXTERNAL,ELIMINATED,STOPLOC
	unsigned long id; // Stop id
//			 long Stopidp; // predecessor Stop id 
	long Edgeid;  // edge at which the stop is located
	double posalong; // pos along the edge
	int StOrdr; // Stop order
	 string StopLbl; // Stop label
	 string StopName; // Stop Name
	short blnHist; // true <>0 if historic
	short blnInbd; // true <>0 if inbound
	short blnExtr; // true <>0 if external
	short blnIncl; // true <>0 if included
	short blnElim; // true <>0 if eliminated
	double CumDist; // supplied Cumulative ride dist
	double CRdTm; // supplied Cumulative ride time
	double undCRdTm; // undelayed Cumulative ride time
	double CRdTmC;  // calculated Cumulative ride time
	double HistOns; 
	double HistOffs; 
	double HistDepVol; 
	double Ons; 
	double Offs; 
	double DepVol; 
	double probStop; 
	double depDelay; 
	double arrDelay; 
	double rideDelay; 
	double PVal;
	double AVal;
	double CRdTmE; 
	double WkTmOns; 
	double WkTmOffs;
	double WalkCost;
	double RideCost;
	double OperCost;
	double TCost;
	double xc;
	double yc;
	string  f1, fldval,fldhdr;
	typedef map<int, string, less<int>> mapflds;
	mapflds :: iterator fld1map_Iter, fld2map_Iter;
	typedef pair < int, string >  fld_pair;
	mapflds mapflds1;
	mapflds::iterator mapfldsit;
	mapflds1 = recoread(rec1);

	mapfldsit = mapflds1.begin();
	while ( mapfldsit != mapflds1.end())
	{
	i = mapfldsit->first;
	fldval = mapfldsit->second;
	//	val1 = fldval.c_str();
	fldhdr = maphdrit->second;

	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "OBJECTID" || f1 == "ID" ) // o id
	{  
		if (fldval.length()>0) {
			id = fromString<long>(fldval); //strtol(val1,&stop1,ib);
			if (id>0) {
			  o1.set_id(id);
			} else {
				return o1;
			}
		} else {
		  rno++;
		  o1.set_id(rno);
		  id = o1.get_id();
		}
	}
	else if (f1 == "ORDER" || f1 == "ORDER_") // Order
	{
		  StOrdr = fromString<int>(fldval); //strtol(val1,&stop1,ib);
		  o1.set_StOrdr(StOrdr);
	}
	else if (f1 == "EDGEID" || f1 == "EID") // EdgeID
	{
		  Edgeid = fromString<int>(fldval); //strtol(val1,&stop1,ib);
		  o1.set_Edgeid(Edgeid);
	}
	else if (f1 == "PALONG" || f1 == "POSALONG") // EdgeID
	{
		  posalong = fromString<double>(fldval); //position along the edge for the stop
		  o1.set_posalong(posalong);
	}
	else if (f1 == "LABEL" || f1 == "STOPIDTEXT" || f1 == "STOP_ID" || f1 == "STOP") // stop label
	{
		StopLbl = fldval;
		o1.set_StopLbl(StopLbl);
	}
	else if (f1 == "STOPNAME" || f1 == "STOP_NAME" )  //stop name
	{
		StopName = fldval;
		o1.set_StopName(StopName);
	  }
	else if (f1 == "CUM_DIST") // cum. distance
	{
		CumDist = fromString<double>(fldval); // strtod(val1,&stop1);
		o1.set_CumDist(CumDist);
	} 
	else if (f1 == "CUM_TIME") // cum. ride time
	{ 
		CRdTm = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_CRdTm(CRdTm);
	}
	else if (f1 == "HISTORIC")
	{
		blnHist = fromString<short>(fldval); // strtol(val1,&stop1,ib);
		o1.set_blnHist(blnHist);
	} 
	else if (f1 == "INBOUND") // Inbound
	{
		blnInbd = fromString<short>(fldval); //strtol(val1,&stop1,ib);
		o1.set_blnInbd(blnInbd);
	}
	else if (f1 == "INCLUDE") // Include
	{
		blnIncl = fromString<short>(fldval); //strtol(val1,&stop1,ib);
		o1.set_blnIncl(blnIncl);
	}
	else if (f1 == "HISTONS") // HistOns
	{
		HistOns = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_HistOns(HistOns);} 
	else if (f1 == "HISTOFFS") // HistOffs
	{	HistOffs = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_HistOffs(HistOffs);
	} 
	else if (f1 == "HISTDEPVOL") // HistDepVol
	{	
		HistDepVol = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_HistDepVol(HistDepVol);
	} 
	else if (f1 == "ONS") // Ons
	{
		Ons = fromString<double>(fldval); 
		o1.set_Ons(Ons);} 
	else if (f1 == "OFFS") // Offs
	{	
		Offs = fromString<double>(fldval);
		o1.set_Offs(Offs);
	} 
	else if (f1 == "DEPVOL") // DepVol
	{	
		DepVol = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_DepVol(DepVol);
	} 
	else if (f1 == "PROBSTOP") // DepVol
	{	
		probStop = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_probStop(probStop);
	}
	else if (f1 == "DEPDELAY") // Departure Delay
	{	
		depDelay = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_depDelay(depDelay);
	}
	else if (f1 == "ARRDELAY") // Arrival Delay
	{	
		arrDelay = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_arrDelay(arrDelay);
	}
	else if (f1 == "DELAY") // total Delay
	{	
		rideDelay = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_rideDelay(rideDelay);
	}
	else if (f1 == "AVAL") // Attraction sum 
	{	
		AVal = fromString<double>(fldval); 
		o1.set_AVal(AVal);
	}
	else if (f1 == "PVAL") // Production sum 
	{	
		PVal = fromString<double>(fldval); 
		o1.set_PVal(PVal);
	}
	else if (f1 == "CUM_TIME2") // Ride time from stop to the end of the line (for ons) WkTmOns 
	{	
		CRdTmE = fromString<double>(fldval); 
		o1.set_CRdTmE(CRdTmE);
	}
	else if (f1 == "CUM_TIME_P") // cum. ride time
	{ 
		CRdTmC = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_CRdTmC(CRdTmC);
	}
	else if (f1 == "UNDCRTM") // Ride time from stop to the end of the line (for ons) WkTmOns 
	{	
		undCRdTm = fromString<double>(fldval); 
		o1.set_undCRdTm(undCRdTm);
	}
	else if (f1 == "WKTMONS") // Walk time summary for ons
	{	
		WkTmOns = fromString<double>(fldval); 
		o1.set_WkTmOns(WkTmOns);
	}
	else if (f1 == "WKTMOFFS") // Walk time summary for ons
	{	
		WkTmOffs = fromString<double>(fldval); 
		o1.set_WkTmOffs(WkTmOffs);
	}
	else if (f1 == "WALKCOST") // Walk Cost
	{	
		WalkCost = fromString<double>(fldval); 
		o1.set_WalkCost(WalkCost);
	}
	else if (f1 == "RIDECOST") // Ride Cost
	{	
		RideCost = fromString<double>(fldval); 
		o1.set_RideCost(RideCost);
	}
	else if (f1 == "OPERCOST") // Operating Cost
	{	
		OperCost = fromString<double>(fldval); 
		o1.set_OperCost(OperCost);
	}
	else if (f1 == "TCOST") // Total Cost
	{	
		TCost = fromString<double>(fldval); 
		o1.set_TCost(TCost);
	}
	else if (f1 == "EXTERNAL" || f1 == "EXTERNAL_") // External stop
	{	
		blnExtr = fromString<short>(fldval); 
		o1.set_blnExtr(blnExtr);
	} 
	else if (f1 == "ELIMINATE") // ELIMINATE stop
	{	
		blnElim = fromString<short>(fldval); 
		o1.set_blnElim(blnElim);
	} 
	else if (f1 == "XCOORD" || f1 == "XC" || f1 == "X") // X coordinate
	{	
		xc = fromString<double>(fldval);
		o1.set_xc(xc);
	} 
	else if (f1 == "YCOORD" || f1 == "YC" || f1 == "Y") // Y coordinate
	{	
		yc = fromString<double>(fldval);
		o1.set_yc(yc);
	} 

	// get the next fld
	mapfldsit++;
	maphdrit++;
	}
	if (o1.get_id()<=0)
	{
		  rno++;
		  o1.set_id(rno);
	// cout<<"Stop id "<<o1.get_id()<<", Name : "<<o1.get_StopName()<<", Cum Rd Tm "<<o1.get_CRdTm()<<", Historic "<<o1.get_blnHist()<<", Hist Ons "<<o1.get_HistOns()<<", Hist Offs "<<o1.get_HistOffs()<<endl;
	}
	return o1;
}

/*
template <typename o, typename m , typename g>
o& readedge1(string& rec1, o& o1, m& maphdrit,g& gc, char *seps ) //*seps = "\t" 
{
	int i=0; int ib=10; int j=0;
    static long sid=0; // serial id
	//	         class vertexp *vertx[2];
	unsigned long id;  // edge id
	         long eoid; // egde predecessor
			 long evp;  // edge pointer
	         short lbl; // label
         	 double ecost; // edge cost
         	 double scost; // start path cost
         	 double tcost; // total path cost
         	 double slen; // length of shape
         	 double palong; // total path cost
			 long frid;  // edge from id
			 long toid;  // edge to id
			 long orig;  // edge to id
//			 short tway;  // two way or one way indicator
        	string efc; // edge FC pointer
        	string stopoff; //Vertex From FC pointer
        	string stopon; //Vertex To FC pointer
        	string orfc; //Origin FC pointer
        	string ornm; //Origin Name pointer
        	string enote; //Edge Note pointer

// EdgeOId	EdgeSId	EdgeFC	ndStartOId	ndStartFC	ndEndOId	ndEndFC	EdgeWt	PosAlong
//			OriginFC	OriginOId	OriginName	EdgeStartC	EdgeEndCos	Labeled	EdgeNote
//				OBJECTID	shape_Leng
   const int bufsz = 1000; // Buffer size;
//   char  *stop1;
   string  f1, fldval,fldhdr;
//   const char *strx;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
   
  mapflds1 = recoread(rec1);
  mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;
	f1 = maphdrit->second;
//						
//ORIGINFC		ORIGINNAME				
//	
	if (f1 == "EDGEOID" || f1 == "EDGEID" || f1 == "OBJECTID")  // edge id
	{  
		id = fromString<long>(fldval);   //strtol(strx,&stop1,ib);
	  if (eid>0) {
		  o1.set_id(id);
		  o1.set_sid(++sid);
	  }
	  else
	  {
		  return o1;
	  }
	}
	else if (f1 == "EDGESID" || f1 == "ESID")	 // edge SId
	{   
		evp = fromString<long>(fldval); //strtol(strx,&stop1,ib);
		o1.set_evp(evp);
	}
	else if (f1 == "EDGEFC") // edge feature class
	{
		o1.set_efc(efc);
	}
	else if (f1 == "NDSTARTOID" || f1 == "VX1" ) // vertex from id 
	{
		frid = fromString<long>(fldval); //strtol(strx,&stop1,ib);
		o1.set_frid(frid);
	} 
	else if (f1 == "NDSTARTFC") // Start Vertex FC
	{
	   o1.set_stopoff(stopoff);
	}
	else if (f1 == "NDENDOID" || f1 == "VX1") // vertex to id
       {
			toid = fromString<long>(fldval); //strtol(strx,&stop1,ib);
			o1.set_toid(toid);
	  }
	else if (f1 == "EDGEWT" || f1 == "COST") // edge weight(cost)
    {
		ecost = fromString<double>(fldval);   //strtod(strx,&stop1);
		// transform edge cost by the ratio of walkcost to ridescost
		if (gc.get_walkcost()>0 && gc.get_ridecost()>0) {
		    //o1.set_cost(ecost);
		    o1.set_cost(gc.get_walkcost()/gc.get_ridecost()* ecost);
		} else {
		    o1.set_cost(ecost);
		}
	}
	else if (f1 == "POSALONG" || f1 == "PALONG") // position along
	{
		palong = fromString<double>(fldval);   //strtod(strx,&stop1);
		o1.set_palong(palong);
	}
	else if (f1 == "ORIGINOID") // origin id
	{
		orig=fromString<long>(fldval);  //strtol(strx,&stop1,ib);
	    o1.set_orig(orig);
	}
	else if (f1 == "EDGESTARTC" || f1 == "EDGESTARTCOST") // edge start cost
	{
			scost = fromString<double>(fldval); //strtod(strx,&stop1);
		if (gc.get_walkcost()>0 && gc.get_ridecost()>0) {
		    //o1.set_scost(scost);
		    o1.set_scost(gc.get_walkcost()/gc.get_ridecost()* scost);
		} else {
			o1.set_scost(scost);
		}
	}
	else if (f1 == "EDGEENDCOS" || f1 == "EDGEENDCOST" ) // edge end cost
	{
		tcost = fromString<double>(fldval); //strtod(strx,&stop1);
		if (gc.get_walkcost()>0 && gc.get_walkcost()>0) {
		    //o1.set_tcost(tcost);
		    o1.set_tcost(gc.get_walkcost()/gc.get_ridecost()* tcost);
		} else {
			o1.set_tcost(tcost);
		}
	}
	else if (f1 == "LABELED") // label
	{
		o1.set_lbl(lbl=fromString<short>(fldval)); //strtol(strx,&stop1,ib));
	}

	else if (f1 == "EDGENOTE") //Edge Note
	{
		enote = fldval;
	   o1.set_enote(enote);
	}
	else if (f1 == "SHAPE_LENG" || f1 == "LENGTH")	// shape length
	{
		  slen=fromString<double>(fldval); //strtod(strx,&stop1);
		 o1.set_slen(slen);
	}            
   // get the next token
	mapfldsit++;
	maphdrit++;
}
   if (o1.get_orig()>0)
   {
	cout<<"Stop Edge id "<<id<<", orig : "<<orig<<", Cost "<<ecost<<", i "<<frid<<", j "<<toid<<" Ride Time= "<<o1.get_scost()<<endl;
   }

   return o1;
}

*/

template <typename o, typename m, typename k >
o& readparc1(string& rec1, o& o1, m& maphdrit, k& keyFldName, char *seps ) // *seps = "\t" 
{
	  static int rno = 0; //record number

	int i=0; int ib=10; int j=0;
//	         class edgev *edg[2];  //edge pointer for ons & offs
	unsigned long pid;  // o1el id long
        	 string pacid; // Parcel id string
	         long eoid; // egde id
//	         short lbl; // label
//	         short on; // if this is on or off
         	 double dblVal; // a double value
//         	 long lngVal; // a long value
         	 string strVal; // a string value
//         	 double cost; // in/out cost to stop
//         	 double scost; // start path cost
//         	 double tcost; // total path cost
         	 double xc; // x-coordinate
         	 double yc; // x-coordinate
         	 double zc; // x-coordinate
         	 double palong; // total path cost
//	         long evp; // edge pointer id
//			 long frid;  // vertex from id
//			 long toid;  // vertex to id
//			 long origon;  // origin for ons of o1 id
//			 long origoff;  // origin for offs of o1 id
//			 long ons;  // ons strength
//			 long offs;  // offs strength
        	string ptype; // o1el land use type 
        	string pfc; // o1el FC 
        	string pnote; //o1el Note 
//        	o1el* par1;

// PARCELOID,PIDLONG,PTYPE,LANDSF,GROSSAREA,LIVINGAREA,TOTALVAL,LANDVAL,
// BLDGVAL,COMPFACT,PVAL,AVAL,INONS,INOFFS,OUTONS,OUTOFFS,EDGEID,PALONG,XC,YC

   const int bufsz = 1000; // Buffer size;

   string f1, fldval,fldhdr;
   typedef map<int, string, less<int>> mapflds;
// vertex map iterator objects
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;

   mapflds1 = recoread(rec1);
   mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;
	f1 = maphdrit->second;
//	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "PARCELOID" || f1 == "OBJECTID" || f1 == "ID" ) // o1el object id
	{
		if (fldval.length()>0) {
			pid = fromString<long>(fldval); //strtol(val1,&stop1,ib);
			if (pid>0) {
				o1.set_id(pid);
			} else {
				return o1;
			}
		} else {
		  rno++;
		  o1.set_id(rno);
		  pid = o1.get_id();
		}
	}
	else if (f1 == "PID_LONG" || f1 == "PIDLONG" || f1 == "PID" || f1 == "PARCELID") // o1el String Id
	{
		o1.set_pacid(fldval);
	  }
	else if (f1 == "PTYPE" || f1 == "SICTYPE") // o1el land use type - derived from the Land Use KeyFld Name
	{
		o1.set_ptype(fldval);
	  }
	else if (f1 == "LANDSF" || f1 == "LAND_SF") // LAND_SF
	{ 
		dblVal = fromString<double>(fldval); //dblVal = strtod(val1,&stop1);
		o1.set_lndsf(dblVal);
	} 
	else if (f1 == "GROSSAREA" || f1 == "GROSS_AREA" || f1 == "NUMEMPLOYE") // Gross Area - or NumEmployees
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	   o1.set_grarea(dblVal);
	}
	else if (f1 == "LIVINGAREA" || f1 == "LIVING_ARE" || f1 == "NOHOLD" || f1 == "NOHHOLD" || f1 == "HHOLDCOUNT" ) // LIVING_ARE or Num. of persons in a Household
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      o1.set_lvarea(dblVal);
	  }
	else if (f1 == "TOTALVAL" || f1 == "FY2003_TOT") // Total valuation
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_valtot(dblVal);
	}
	else if (f1 == "LANDVAL" || f1 == "FY2003_LAN") // Land valuation
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      o1.set_valand(dblVal);
	  }
	else if (f1 == "BLDGVAL" || f1 == "FY2003_BLD") // Bldg valuation
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_valbld(dblVal);
	  }
	else if (f1 == "COMPFACT") // competition factor
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      o1.set_cfact(dblVal);
		}
	else if (f1 == "PVAL") // Production strength
	{	  dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      o1.set_pval(dblVal);
	  }
	else if (f1 == "AVAL") // o1el attraction strength
	{	  dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      o1.set_aval(dblVal);
	  }
	else if (f1 == "HINONS" || f1 == "HONS" || f1 == "HOUTONS") // in Boardings - inOns
	{		
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_ons(0.0);
		}
	else if (f1 == "HINOFFS" || f1 == "HOFFS" || f1 == "HOUTOFFS")  // Aligtings - inOffs 
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
			o1.set_offs(0.0);
		}
	else if (f1 == "INONS" || f1 == "ONS" || f1 == "OUTONS") // in Boardings - inOns
	{		
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_ons(0.0);
		}
	else if (f1 == "INOFFS" || f1 == "OFFS" || f1 == "OUTOFFS")  // Aligtings - inOffs 
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
			o1.set_offs(0.0);
		}
	else if (f1 == "EDGEID") // Edge Object Id  
	{	
		eoid=fromString<long>(fldval); //strtol(val1,&stop1,ib);
	    o1.set_eoid(eoid);
	  }
	else if (f1 == "PALONG" || f1 == "POSALONG") // position along edge
	{	
		palong = fromString<double>(fldval); //strtod(val1,&stop1);
		o1.set_palong(palong);
	}
	else if (f1 == "XC" || f1 == "XCOORD" || f1 == "X") // X Coordinate 
	{	
		xc = fromString<double>(fldval); 
		o1.set_xc(xc);
	}
	else if (f1 == "YC" || f1 == "YCOORD" || f1 == "Y") // X Coordinate 
	{	
		yc = fromString<double>(fldval); 
		o1.set_yc(yc);
	}
	else if (f1 == "ZC" || f1 == "ZCOORD" || f1 == "Z") // X Coordinate 
	{	
		zc = fromString<double>(fldval); 
		o1.set_zc(zc);
	}
   // get the next token
	mapfldsit++;
	maphdrit++;
}
   if (o1.get_origon()>0)
   {
	// cout<<"Stop o1el id "<<pid<<", orig : "<<orig<<", Cost "<<cost<<", i "<<frid<<", j "<<toid<<endl;
   }

   return o1;
}

#endif // READRECORDATA_H




// WRITE object Data  into a file
#ifndef WRITEOBJECTDATA_H
#define WRITEOBJECTDATA_H

template <typename m, typename  outfile >
 outfile& writeBindblong(m& m1, outfile& pts )
{
m::iterator mit;
long vid;
double mDbl;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     mDbl = mit->first;
     vid = mit->second;
	 pts.write(reinterpret_cast<char *>(&mDbl), sizeof(mDbl));
	 pts.write(reinterpret_cast<char *>(&vid), sizeof(vid));
  }
  pts.close();
  return  pts;
}

template <typename m, typename v1, typename v2, typename  o >
m& writeBindualvar(m& m1,  v1& var1, v2& var2,o& pts)
{
  m::iterator mit;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     var1 = mit->first;
     var2 = mit->second;
	 pts.write(reinterpret_cast<char *>(&var1), sizeof(var1));
	 pts.write(reinterpret_cast<char *>(&var2), sizeof(var2));
  }
  pts.close();
  return  m1;
}


// write an object data using its serialize function
template <typename m, typename o, typename k, typename  outfile >
 outfile& writeBinObjectData(m& m1, o& o1,k& k1, outfile& pOf )
{  
m::iterator mit;

  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     k1 = mit->first;
     o1 = mit->second;
	 o1.serialize(pOf);
  }
  pOf.close();
  return  pOf;
}

// write an object data using its serialize function
template <typename m, typename o, typename q, typename  outfile >
 outfile& writeBinObjectCol(m& m1, o& o1, q& key1, outfile& pOf )
{  
m::iterator mit;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     key1 = mit->first;
	 mit->second.serialize(pOf);
  }
  pOf.close();
  return  pOf;
}

 // write an object data in text fromat using its serializetext function
template < typename m, typename o, typename p, typename  q, typename s >
q& writeTextObjectData(m& m1, o& o1, p& p1, q& pOf, s& txt )
{  
m::iterator mit;
p& key=p1;
if (txt!="") {
	pOf <<txt<<endl;
}
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     key = mit->first;
     //o1 = mit->second;
	 mit->second.serializetext(pOf);
  }
  pOf.close();
  return  pOf;
}

// write multi-map id and related object data using its serialize function
template <typename m, typename o, typename  outfile >
outfile& writeBinIdplusObjectData(m& m1, o& o1, outfile& pOf )
{  
m::iterator mit;
long vid;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     vid = mit->first;
	 pOf.write(reinterpret_cast<char *>(&vid),sizeof(vid));
     o1 = mit->second;
	 o1.serialize(pOf);
  }
  pOf.close();
  return  pOf;
}

// write an the id of an object in text fromat with a text hearder if included  
template <typename m, typename o, typename  outfile,typename s>
outfile& writeTextDPatternData(m& m1, o& o1, outfile& pOf, s& txt)
{  
m::iterator mit;
long id;
if (txt!="") {
	pOf <<txt<<endl;
}
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     id = mit->first;
     o1 = mit->second;
	 pOf <<o1.get_id()<<"\t";
  }
  pOf<<endl;
  return  pOf;
}

// write a multimap data in text format with a text hearder if included
template <typename m, typename o,typename p, typename  outfile,typename s>
outfile& writeTextData(m& m1, o& o1,p& p1, outfile& pOf, s& txt)
{  
m::iterator mit;
if (txt!="") {
	pOf <<txt<<endl;
}
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
     o1 = mit->first;
     p1 = mit->second;
	 pOf <<o1<<"\t"<<p1<<endl;
  }
  pOf<<endl;
  return  pOf;
}


#endif // WRITEOBJECTDATA_H ///:~

// Read object data into the respective object
#ifndef READATAINTOBJECT_H
#define READATAINTOBJECT_H

template <typename m, typename  infile >
 m& readBindblong(m& m1, infile& pts )
{
  
long vid;
double mDbl;
  do
  {
	 pts.read(reinterpret_cast<char *>(&mDbl), sizeof(mDbl));
	 pts.read(reinterpret_cast<char *>(&vid), sizeof(vid));
	 if (!pts.eof()) {
     m1.insert(dblng_Pair (mDbl,vid));
	 }
  } while (!pts.eof());
  pts.close();
  return  m1;
}

template <typename m, typename v1, typename v2, typename  infile >
m& readBindualvar(m& m1,  v1& var1, v2& var2,infile& pts)
{
  typedef pair <v1,v2> obj_pair;
 
  do
  {
	 pts.read(reinterpret_cast<char *>(&var1), sizeof(var1));
	 pts.read(reinterpret_cast<char *>(&var2), sizeof(var2));
	 if (!pts.eof()) {
     m1.insert(obj_pair (var1,var2));
	 }
  } while (!pts.eof());
  pts.close();
  return  m1;
}

// read a vertex object data using variables
template <typename m, typename o, typename  infile >
m& readBinVertexObjectDatax(m& m1, o& o1, infile& pIf)
{
long id;
  do
  {
	 //o1.deserialize(pIf);
	unsigned long pid;  // vertex id
	         int pidp; // predecessor
	         short lbl; // label
         	 double cost;  // cost at vertex
         	 double tcost; // total cost
			 long fid;  //  
			 long orig;  // Origin for vertex least cost path
			 short toStop;  // travel direction toStop=1 (on) or from stop=0 (off) 

	 pIf.read(reinterpret_cast<char *>(&pid), sizeof(pid));
	 pIf.read(reinterpret_cast<char *>(&pidp), sizeof(pidp));
	 pIf.read(reinterpret_cast<char *>(&lbl), sizeof(lbl));
	 pIf.read(reinterpret_cast<char *>(&cost), sizeof(cost));
	 pIf.read(reinterpret_cast<char *>(&tcost), sizeof(tcost));
	 pIf.read(reinterpret_cast<char *>(&fid), sizeof(fid));
	 pIf.read(reinterpret_cast<char *>(&orig),sizeof(orig));
	 pIf.read(reinterpret_cast<char *>(&toStop), sizeof(toStop));
	 o1.vxp(pid,pidp,lbl,cost,fid,orig,tcost,toStop);
	 if (!pIf.eof()) {
		 m1.insert(vx_pair(id,o1));
	 }
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}

// read a vertex object data using its deserialize function
template <typename m, typename o, typename  infile >
m& readBinVertexObjectData(m& m1, o& o1, infile& pIf)
{
  long id;
  vertexp vx;

  do
  {
	 if (!pIf.eof()) {
	 vx.deserialize(pIf);
	 id=vx.get_id();
		 m1.insert(vx_pair(id,vx));
	 }
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}


// read an edge object data using its deserialize function
template <typename m, typename o, typename  infile >
m& readBinEdgeObjectData(m& m1, o& o1, infile& pIf, short onoff=1)
{
long id;
//edgev ev;
  do
  {
	  if (pIf.is_open()) {
		o1.deserialize(pIf);
		//o1.set_toStop(onoff);
		id=o1.get_esid();

		if (!pIf.eof()) {
			 m1.insert(ed_Pair(id,o1));
			}
	  } else {
		  cout<< "file not open!!";
		  break;
	  }
	  
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}

// read an edge object data using its deserialize function
template <typename m, typename o, typename i, typename  infile >
m& readBinIdMap(m& m1, o& o1, i& i1 , infile& pIf)
{
	typedef pair <i,o> objPair;
//edgev ev;
	i id;
	do
	{
	  if (pIf.is_open()) {
		o1.deserialize(pIf);
		id=o1.get_id();

		if (!pIf.eof()) {
			 m1.insert(objPair(id,o1));
			}
	  } else {
		  cout<< "file not open!!";
		  break;
	  }
	  
	} while (!pIf.eof());
	pIf.close();
	return  m1;
}



// read a parcel object data using its deserialize function
template <typename m, typename o, typename  infile, typename h >
m& readBinParcObjectData(m& m1, o& o1, infile& pIf, h& hi)
{
long id;
parcel parc;
  do
  {
	 parc.deserialize(pIf);
	 id= parc.get_eoid();
	 parc.set_hist(hi);
	 if (!pIf.eof()) {
		 if (id<0) {
			 id = 0;
		 }
		 m1.insert(PE_Pair(id,parc));
	 }
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}

// read a multi-map long id and the associated parcel object data using its deserialize function
template <typename m, typename o, typename  infile >
m& readBinIdParcObjectData(m& m1, o& o1, infile& pIf)
{
long id;
parcel parc;
  do
  {
	 pIf.read(reinterpret_cast<char *>(&id), sizeof(id));
	 if (!pIf.eof()) {
		 parc.deserialize(pIf);
	 }
//	 id= parc.get_eoid();
	 if (!pIf.eof()) {
		 if (id<0) {
			 id = 0;
		 }
		 m1.insert(PE_Pair(id,parc));
	 }
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}


// read a vertex object data using its deserialize function
template <typename m, typename o, typename  infile >
m& readBinStopObjectData(m& m1, o& o1, infile& pIf)
{
  long id;
  tstop stp;

  do
  {
	 if (!pIf.eof()) {
	 stp.deserialize(pIf);
	 id=stp.get_id();
		 m1.insert(tsidpair(id,stp));
	 }
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}

// read a dp stops object data using its deserialize function
template <typename m, typename o,typename q, typename  infile >
m& readBinDPStopObjectData(m& m1, o& o1, q& q1, infile& pIf)
{
  typedef pair <q,o> qopair;

  do
  {
	 if (!pIf.eof()) {
	 o1.deserialize(pIf);
	 q1=o1.get_dpkey();
		 m1.insert(qopair(q1,o1));
	 }
  } while (!pIf.eof());
  pIf.close();
  return  m1;
}

// read a dp stops Immediate Cost data from the ArcGIS VB output file
template <typename m, typename o,typename q, typename  infile >
m& readtstopMapData(m& m1, o& o1, q& q1, infile& pIf,char *seps)
{
	o* op1;
    typedef pair <q,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
// read stop data file
	if (!pIf.eof()) {
	    getline(pIf,rec1);
        mapfldshdr = readatahdr (rec1,seps);
	}


while(!pIf.eof())
{ //1w
    getline(pIf,rec1);
	if (rec1.length() >0)
	{
        fld1map_Iter = mapfldshdr.begin();
		o o2;	
		op1=&o2;
		op1 = &readtstop1(rec1,*op1,fld1map_Iter,seps);
     	ip++;
		if (op1->get_id()<=0) {
			q1=ip;
			op1->set_id(ip);
		} else {
			q1=op1->get_id();
		}
		m1.insert(qopair(q1,o2));
	}

}
  pIf.close();
  return  m1;

}

// read a edge data into a map container
template <typename m, typename o,typename q, typename  infile, typename g >
m& readedgeMapData(m& m1, o& o1, q& q1, infile& pIf,g& gc,char *seps)
{
	o* op1;
	op1=&o1;
    typedef pair <q,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
	double ecost=0;
// read stop data file
	if (!pIf.eof()) {
	    getline(pIf,rec1);
        mapfldshdr = readatahdr (rec1,seps);
	}


while(!pIf.eof())
{ //1w
    getline(pIf,rec1);
	if (rec1.length() >0)
	{
        fld1map_Iter = mapfldshdr.begin();
		op1 = &readedge1(rec1,*op1,fld1map_Iter,gc,seps);
		if ( op1->get_slen()>0) {
				if (gc.get_walkspd() != 0) {
					ecost = op1->get_slen()/(gc.get_walkspd()*60); // convert the walk length to minutes
					if (gc.get_walkcost()> 0 && gc.get_ridecost()>0) {
						ecost *= (gc.get_walkcost()/gc.get_ridecost()); // multiply the ride to walk factor
					}
						op1->set_cost(ecost);
				}
		}
     	ip++;
		q1=op1->get_id();
		m1.insert(qopair(q1,o1));
	}

}
  pIf.close();
  return  m1;

}


template <typename a, typename q, typename k,typename o,typename m,  typename d, typename g >
m& readedgeTableMapData(a& netdb, q& q1, k& k1,o& o1, m& m1, d& logFile,g& gc)
{
	o* op1;
	op1=&o1;
    typedef pair <k,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
	double ecost=0;
	int ret;
	int i=0, j=0;
	sqlite3_stmt *stmt = NULL;
	char **results;
	int rows=0;
	int columns=0 , col=0;
	char *zErrMsg = NULL;
	string pkName="", sstopi="", sFldVal="";
	int pkCount = 0;
	int fldNo = 0;
	int iFldVal=0, stopid=0,stopi=0;
	long lFldVal=0;
	double dFldVal=0,mcost=0,dschlTime=0;
	//ret = sqlite3_get_table( netdb, q1.c_str(), &results, &rows, &columns, &zErrMsg );
	if ( sqlite3_prepare_v2( netdb, q1.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
	{
	  // some error occurred
	  fprintf(logFile, "SQLite error: %s\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
	  fprintf(stderr, "SQLite error: %s\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
	  return m1;
	}
	m1.clear();
	// query the stop detail data from the stop table

	while ( sqlite3_step( stmt ) == SQLITE_ROW )
	{
		rows++;
		op1 = &o1;
		// query the stop detail data from the stop table
		// j=0 // edgeid
		col=0;
		iFldVal = sqlite3_column_int(stmt, col);
		if (iFldVal>0) {  // id
			op1->set_id(iFldVal);
		}
		// (j==1) // Edge serial Id
		col=1;
		iFldVal = sqlite3_column_int(stmt, col);
		if (iFldVal>0) {  // id
			op1->set_esid(iFldVal);
		}
		// (j==2) // vx from 
		col=2;
		iFldVal = sqlite3_column_int(stmt, col);
		if (iFldVal>0) {  // id
			op1->set_frid(iFldVal);
		}
		// (j==3) // Vx to
		col=3;
		iFldVal = sqlite3_column_int(stmt, col);
		if (iFldVal>0) {  // id
			op1->set_toid(iFldVal);
		}
		// (j==4) // length
		col=4;
		dFldVal = sqlite3_column_double(stmt, col);
		if (dFldVal>0) {  // id
			op1->set_slen(dFldVal);
		}
		//  (j==5) // Street Name
		col=5;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
			op1->set_ornm(sFldVal);
		}
		// (j==6) // Cost 
		col=6;
		dFldVal = sqlite3_column_double(stmt, col);
		if (dFldVal>0) {  // cost
			op1->set_cost(dFldVal);
		}
		//  (j==7) // speed limit
		col=7;
		dFldVal = sqlite3_column_double(stmt, col);
		// 	(j==8) // comment 
		col=8;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
			op1->set_enote(sFldVal);
		}
		//  (j==9) // class
		col=9;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
		}
		// (j==10) // ROW Width
		col=10;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
		}
		// (j==11) // Face of Curb width 
		col=11;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
		}
		// (j==12) // one way  
		col=12;
		iFldVal = sqlite3_column_int(stmt, col);
		op1->set_tway(iFldVal);
		// 	(j==13) // Private road? 
		col=13;
		iFldVal = sqlite3_column_int(stmt, col);
		//	(j==14) // Functional Class   
		col=14;
		iFldVal = sqlite3_column_int(stmt, col);
		// (j==15) // Facility Type (one way or two way 
		col=15;
		iFldVal = sqlite3_column_int(stmt, col);
		// 	(j==16) // Path width (bike path)
		col=16;
		iFldVal = sqlite3_column_int(stmt, col);
		// 	(j==17) // Parking Width 
		col=17;
		dFldVal = sqlite3_column_double(stmt, col);
		//  (j==18) // Number Lane 
		col=18;
		iFldVal = sqlite3_column_int(stmt, col);
		//  (j==19) // Illegal Parking Act 
		col=19;
		iFldVal = sqlite3_column_int(stmt, col);
		// (j==20) // ADT  
		col=20;
		iFldVal = sqlite3_column_int(stmt, col);
		//  (j==21) // Right Lane count 
		col=21;
		iFldVal = sqlite3_column_int(stmt, col);
		// (j==22) // Signal
		col=22;
		iFldVal = sqlite3_column_int(stmt, col);
		// (j==23) // Median Width
		col=23;
		dFldVal = sqlite3_column_double(stmt, col);
		//	(j==24) // Shared Space width
		col=24;
		dFldVal = sqlite3_column_double(stmt, col);
		// (j==25) // Pressence of Center Line
		col=25;
		iFldVal = sqlite3_column_int(stmt, col);
		//  (j==26) // Traffic Stress value
		col=26;
		iFldVal = sqlite3_column_int(stmt, col);
		op1->set_lts(iFldVal);
		// (j==27) // Overide Value for Stress
		col=27;
		iFldVal = sqlite3_column_int(stmt, col);
		if (iFldVal > 0) {
			op1->set_lts(iFldVal);
		}
		// (j==2) // Geometry
		if ( op1->get_slen()>0) {
				if (gc.get_walkspd() != 0) {
					ecost = op1->get_slen()/(gc.get_walkspd()*60); // convert the walk length to minutes
					if (gc.get_walkcost()> 0 && gc.get_ridecost()>0) {
						ecost *= (gc.get_walkcost()/gc.get_ridecost()); // multiply the ride to walk factor
					}
						op1->set_cost(ecost);
				}
		}
		ip++;
		if (op1->get_id())
		{
			m1.insert(qopair(op1->get_esid(),o1));
		}
	} // loop over rows
	sqlite3_finalize(stmt );

	return  m1;

}


template <typename a, typename q, typename k,typename o,typename m,  typename d >
m& readvertexTableMapData(a& netdb, q& q1, k& k1,o& o1, m& m1, d& logFile)
{
	o* op1;
	op1=&o1;
    typedef pair <k,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
	double ecost=0;
	int ret;
	int i=0, j=0;
	sqlite3_stmt *stmt = NULL;
	char **results;
	int rows=0;
	int columns=0,col=0;
	char *zErrMsg = NULL;
	string pkName="", sstopi="", sFldVal="";
	int pkCount = 0;
	int fldNo = 0;
	int iFldVal=0, stopid=0,stopi=0;
	long lFldVal=0;
	double dFldVal=0,mcost=0,dschlTime=0;
//	ret = sqlite3_get_table( netdb, q1.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, q1.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(logFile, "\nSQLite error: %s\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
		  return m1;
		}
		m1.clear();
		// query the stop detail data from the stop table

		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			rows++;
			op1 = &o1;
			// query the vertex detail data from the vertex table
			//		if (j==0) // Vertex Id
			col=0;
			iFldVal = sqlite3_column_int(stmt, col);
			if (iFldVal>0) {  // id
				op1->set_id(iFldVal);
			}
			//  (j==1) // Edge Id
			col=1;
			iFldVal = sqlite3_column_int(stmt, col);
			// (j==2) // Idp  
			col=2;
			iFldVal = sqlite3_column_int(stmt, col);
			op1->set_idp(iFldVal);
			// (j==3) // vertex cost
			col=3;
			dFldVal = sqlite3_column_double(stmt, col);
			op1->set_cost(inf);
			op1->set_tcost(inf);
			// (j==4) // signal presence value > 0
			col=4;
			iFldVal = sqlite3_column_int(stmt, col);
			// (j==5) // x coordinates
			col=5;
			dFldVal = sqlite3_column_double(stmt, col);
			op1->set_x(dFldVal);
			// (j==6) // y coordinate
			col=6;
			dFldVal = sqlite3_column_double(stmt, col);
			op1->set_y(dFldVal);
			op1->set_z(0);
     			ip++;
			if (op1->get_id())
			{
				Point p1(op1->get_x(),op1->get_y());
				op1->set_pt(p1);
				m1.insert(qopair(op1->get_id(),o1));
			}
		} // loop over rows
		sqlite3_finalize( stmt );

  return  m1;

}



// read Vertex data into map container
template <typename m, typename o,typename q, typename  infile >
m& readvertexMapData(m& m1, o& o1, q& q1, infile& pIf,char *seps)
{
	o* op1;
	op1=&o1;
    typedef pair <q,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
// read stop data file
	if (!pIf.eof()) {
	    getline(pIf,rec1);
        mapfldshdr = readatahdr (rec1,seps);
	}


while(!pIf.eof())
{ //1w
    getline(pIf,rec1);
	if (rec1.length() >0)
	{
        fld1map_Iter = mapfldshdr.begin();
		op1 = &readvertex1(rec1,*op1,fld1map_Iter,seps);
     	ip++;
		q1=o1.get_id();
		m1.insert(qopair(q1,o1));
	}

}
  pIf.close();
  return  m1;

}

// read a edge data into a map container
template <typename m, typename o,typename q, typename  infile,typename k >
m& readparcelMapData(m& m1, o& o1, q& q1, infile& pIf,k& keyFld,char *seps)
{
	o* op1;
	op1=&o1;
    typedef pair <q,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
// read stop data file
	if (!pIf.eof()) {
	    getline(pIf,rec1);
        mapfldshdr = readatahdr (rec1,seps);
	}


while(!pIf.eof())
{ //1w
    getline(pIf,rec1);
	if (rec1.length() >0)
	{
        fld1map_Iter = mapfldshdr.begin();
		op1 = &readparc1(rec1,*op1, fld1map_Iter,keyFld,seps);
     	ip++;
		q1=op1->get_id();
		if (q1>0) {
			q1=op1->get_eoid();
			m1.insert(qopair(q1,o1));
		}
	}

}
  pIf.close();
  return  m1;

}



#endif // READATAINTOBJECT_H ///:~// Read object data into the respective object


#ifndef EDGEOBJECTID2EDGEIDANDVISEVERSA_H
#define EDGEOBJECTID2EDGEIDANDVISEVERSA_H

// make a new map of any object (with a get_id function) using the id as the key  
template <typename m, typename n, typename o, typename r >
n& remapObj2Id( m& m1, n& n1, o& o1,r& r1)
{
  m::iterator mit;
  typedef pair <r,o> obj_pair;
  n1.clear();
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  r1 = o1.get_id();
	 n1.insert(obj_pair(r1,o1));
  } 
  return  n1;
}

// make a new map the edge object id as the key of the collection
template <typename m, typename n, typename o, typename p, typename r >
n& remapEdgeId2ObjId( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <r,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  r1 = o1.get_eoid();
	 n1.insert(obj_pair(r1,o1));
  } 
  return  n1;
}

// make a new map using the edge id as the key for the edge object id
template <typename m, typename n, typename o, typename p, typename r >
n& funcEdgeId2ObjId( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <p,r> obj_pair;
  n1.clear();
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  r1 = o1.get_esid();
	  p1 = o1.get_id();
	 n1.insert(obj_pair(p1,r1));
  } 
  return  n1;
}

// make a new map using the edge object id as the key for the edge id
template <typename m, typename n, typename o, typename p, typename r >
n& funcEdgeObjId2Id( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <r,p> obj_pair;
  n1.clear();
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  r1 = o1.get_esid();
	  p1 = o1.get_id();
	 n1.insert(obj_pair(r1,p1));
  } 
  return  n1;
}

// make a new map for the edge vertices from and to ids with the edge object id 
template <typename m, typename n, typename o, typename p, typename r >
n& funcV1V2EdgeId( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <r,p> obj_pair;
  n1.clear();
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_id();
	  r1 = o1.get_frid();
	  n1.insert(obj_pair(r1,p1));
	  r1 = o1.get_toid();
	  n1.insert(obj_pair(r1,p1));
  } 
  return  n1;
}


// make a new map for the edge vertices from and to ids with the edge object id 
template <typename m, typename n, typename o, typename p, typename r >
n& funcEdgeIdV1V2( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <r,p> obj_pair;
  n1.clear();
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_id();
	  r1 = o1.get_frid();
	  n1.insert(obj_pair(p1,r1));
	  r1 = o1.get_toid();
	  n1.insert(obj_pair(p1,r1));
  } 
  return  n1;
}

// make a new map for the edge vertices with from and to ids and the edge object id 
template <typename m, typename n, typename o, typename p, typename r >
n& funcV1V2EOId( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <r,p> obj_pair;
  n1.clear();
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_esid();
	  r1 = o1.get_frid();
	  n1.insert(obj_pair(r1,p1));
	  r1 = o1.get_toid();
	  n1.insert(obj_pair(r1,p1));
  } 
  return  n1;
}



// make a new map by of from and to id's of an edge
template <typename m, typename n, typename o, typename p, typename r >
n& funcV1V2( m& m1, n& n1, o& o1,p& p1, r& r1)
{
  m::iterator mit;
  typedef pair <r,p> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_frid();
	  r1 = o1.get_toid();
	  n1.insert(obj_pair(r1,p1));
	  n1.insert(obj_pair(p1,r1));
  } 
  return  n1;
}

// make a new map by changing the vertex id key to vertex orig key of the original map
template <typename m, typename n, typename o, typename p >
n& remapVertId2Orig( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_orig();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

// make a new map by changing the vertex id key to vertex index of the original map
template <typename m, typename n, typename o, typename p >
n& remapVertId2Index( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_index();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

// make a new map by changing the id of a map list to origin ons
template <typename m, typename n, typename o, typename p >
n& remaParc2OrigOn( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_origon();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

// make a new map by changing the id of a map list to origin offs
template <typename m, typename n, typename o, typename p >
n& remaParc2OrigOff( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.get_origoff();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

// make a new map by changing the trip time to be the sorting key
template <typename m, typename n, typename o, typename p >
n& remap2TripTime( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.tripTime();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

// make a new map by changing the trip time to be the sorting key
template <typename m, typename n, typename o, typename p >
n& remap2TripNumber( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.tripNumber();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

// make a new map by changing the trip time to be the sorting key
template <typename m, typename n, typename o, typename p >
n& remap2TripCount( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <p,o> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	  o1 = mit->second;
	  p1 = o1.stopCount();
	 n1.insert(obj_pair(p1,o1));
  } 
  return  n1;
}

#endif // EDGEOBJECTID2EDGEIDANDVISEVERSA_H ///:~



#ifndef REVERSEMAPKEYS_H
#define REVERSEMAPKEYS_H

// make a new map by changing the key of the original map
template <typename m, typename n, typename o, typename p >
n& reverseMapKeys( m& m1, n& n1, o& o1,p& p1)
{
  m::iterator mit;
  typedef pair <o,p> obj_pair;
  for(mit=m1.begin();mit!=m1.end();mit++)
  {
	 n1.insert(obj_pair(mit->second,mit->first));
  } 
  return  n1;
}

#endif // REVERSEMAPKEYS_H ///:~

// Aggregate parcel data into a DEQUE array based on stop id.
#ifndef PARCELAGGREGATE_H
#define PARCELAGGREGATE_H
//maplngstop-s,mmaplngpar-p,ParcelOffAccum-t
template < typename S, typename P, typename T> 
deque<T> pagwalkstop (S& tsidmap, P& mmaparStop, T& pa)
{
    long ip;
	tstop* pstop;
    parcel* parp;
 deque<parcel> vpar;
 deque<T> vpa;
 deque<T>::iterator vpait;
 tsidmit = tsidmap.begin();

// aggregate parcel walk time impact by transit stop

///*
 while (tsidmit !=tsidmap.end()) {
	 pstop = &(tsidmit->second);
	 ip = pstop->get_id();
// get the range of parcel objects of the current stop id
 		pair<mmapE_AIter, mmapE_AIter> pastrange = mmaparStop.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mmAPEi=pastrange.first; mmAPEi!=pastrange.second;mmAPEi++)
	    {
    	   parp = &(mmAPEi->second);
   // compute the on/off strengths per parcel
				vpar.push_back((*parp));
		}
		
  vpa.push_back(for_each(vpar.begin(),vpar.end(), T()));
  vpar.clear();
  tsidmit++;
 }
return vpa;
}

#endif // PARCELAGGREGATE_H ///:~


#ifndef PARCELPACALCONOFF_H
#define PARCELPACALCONOFF_H
//maplngstop-s,mmaplngpar-p,ParcelOffAccum-t
template < typename S, typename P, typename lu, typename gc,typename onf, typename hst> 
P& pacalc (S& stops, P& mmaparStop, lu* luci, gc* gc1,onf& onoff,
					  hst& blnHist)
{
    long ip;
	tstop* pstop;
    parcel* parp;
 double propoff=0,propon=0;
 tsidmit = stops.begin();
// aggregate parcel walk time impact by transit stop

 while (tsidmit !=stops.end()) {
		pstop = &(tsidmit->second);
		ip = pstop->get_id();
// get the range of parcel objects of the current stop id
 		pair<mmapP_AIter, mmapP_AIter> pastrange = mmaparStop.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mmAPEi=pastrange.first; mmAPEi!=pastrange.second;mmAPEi++)
	    {
    	   parp = &(mmAPEi->second);
   // compute the on/off strengths per parcel
		   if (blnHist) 
		   {
			if (onoff == 0 && parp->get_hwkoffcost()>=0) {
					propoff = exp(-gc1->get_propensity()/300*parp->get_hwkoffcost()*gc1->get_walkspd()*60);
			} else	if (onoff == 1 && parp->get_hwkoncost()>=0) {
					propon = exp(-gc1->get_propensity()/300*parp->get_hwkoncost()*gc1->get_walkspd()*60);
			} else {
				propoff=propon=0;
			}
			if (propon!=0 || propoff!=0) {
			   string strPtype = parp->get_ptype();
		       lucmit = lumap.find(strPtype);
			   if (lucmit!=lumap.end()) {
			     luci = &(lucmit->second);
				 parp->set_lucd(luci->get_LUC());
	 		     if (luci->get_KeyProp()=="Gross_Area" || luci->get_KeyProp()=="GROSS_AREA"
					 || luci->get_KeyProp()=="GROSSAREA" || luci->get_KeyProp()=="NumEmploye" 
					 || luci->get_KeyProp()=="NumberOfEmployees" ) {
						if (onoff) {
							propon=(luci->get_OnCoeff()*parp->get_grarea()*parp->get_cfact()* propon);
						}
						else 
						{
							propoff=(luci->get_OffCoeff()*parp->get_grarea()*parp->get_cfact()*propoff);
						}
				 }
				 else if (luci->get_KeyProp()=="Living_are" || luci->get_KeyProp()=="LIVINGAREA" 
					 || luci->get_KeyProp()=="LIVING_ARE" || luci->get_KeyProp() == "NOHOLD" 
					 || luci->get_KeyProp() == "NoHold" || luci->get_KeyProp() == "OCCCOUNT" 
					 || luci->get_KeyProp() == "POP_2000") {
						if (onoff) {
							propon = (luci->get_OnCoeff() * parp->get_lvarea()*parp->get_cfact()* propon);
						}
						else {
							propoff = (luci->get_OffCoeff() * parp->get_lvarea()*parp->get_cfact()* propoff);
						}
				 }  else if ( luci->get_KeyProp() == "POP") {  // blocks with population and gross area as scale factor
						if (onoff) {
							propon = (luci->get_OnCoeff() * parp->get_grarea()*parp->get_cfact()* propon);
						}
						else {
							propoff = (luci->get_OffCoeff() * parp->get_grarea()*parp->get_cfact()* propoff);
						}
				 }
				if (onoff) {
					parp->set_pval(propon);
				} else {
					parp->set_aval(propoff);
				}
			   } else  { // skip parcel since land use code data is not found
				cout<<"\nLand use code data not found for Type "<<strPtype<< " for parcel "<<parp->get_pacid();
			   }
			}
		   }
		}
  // aggregate the parcel data		
  tsidmit++;
 }
return mmaparStop;
}

#endif // PARCELPACALCONOFF_H ///:~

#ifndef PARCELAGGREGATEONOFF_H
#define PARCELAGGREGATEONOFF_H
//maplngstop-s,mmaplngpar-p,ParcelOffAccum-t
template < typename S, typename P, typename T> 
deque<T> pagstop (S& tsidmap, P& mmaparStop, T& pa)
{
 long ip;
 tstop* pstop;
 parcel* parp;
 deque<parcel> vpar;
 deque<T> vpa;
 deque<T>::iterator vpait;
 S::iterator sit;
 sit = tsidmap.begin();
// aggregate parcel walk time, on/off demand, pval/aval by transit stop

///*
 while (sit !=tsidmap.end()) {
		pstop = &(sit->second);
		ip = pstop->get_id();
// get the range of parcel objects of the current stop id
 		pair<mmapP_AIter, mmapP_AIter> pastrange = mmaparStop.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mmAPEi=pastrange.first; mmAPEi!=pastrange.second;mmAPEi++)
	    {
    		parp = &(mmAPEi->second);
			vpar.push_back((*parp)); // if alternative run include the parcel
		}
  // aggregate the parcel data		
  vpa.push_back(for_each(vpar.begin(),vpar.end(),T()));
  vpar.clear();
  sit++;
 }
return vpa;
}


template < typename S, typename M, typename T,typename N> 
deque<T> aggrebystop (S& s1, M& m1, T& t1, N& n1)
{
 long ip;
 T t2;
 T* pt2=&t2;
 N* pn2;
 typedef S::iterator sit;
 sit sit1;
 typedef M::iterator mit;
 mit mit1;

 deque<N> vpar;
 deque<T> vpa;
 deque<T>::iterator vpait;
 sit1 = s1.begin();
// aggregate parcel walk time, on/off demand, pval/aval by transit stop

///*
 while (sit1 !=s1.end()) {
		pt2 = &(sit1->second);
		ip = pt2->get_id();
// get the range of parcel objects of the current stop id
 		pair<mit, mit> pastrange = m1.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mit1=pastrange.first; mit1!=pastrange.second;mit1++)
	    {
    		pn2 = &(mit1->second);
			vpar.push_back((*pn2)); // if alternative run include the parcel
		}
  // aggregate the parcel data		
  vpa.push_back(for_each(vpar.begin(),vpar.end(),T()));
  vpar.clear();
  sit1++;
 }
return vpa;
}

template < typename S, typename M, typename T,typename N, typename A> 
deque<A> aggbystop (S& s1, M& m1, T& t1, N& n1,A& a1)
{
 long ip;
 T t2;
 T* pt2=&t2;
 N* pn2;
 typedef S::iterator sit;
 sit sit1;
 typedef M::iterator mit;
 mit mit1;

 deque<N> vpar;
 deque<A> vpa;
 deque<A>::iterator vpait;
 sit1 = s1.begin();
// aggregate vertex by transit stop
 while (sit1 !=s1.end()) {
		pt2 = &(sit1->second);
		ip = pt2->get_id();
// get the range of parcel objects of the current stop id
 		pair<mit, mit> pastrange = m1.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mit1=pastrange.first; mit1!=pastrange.second;mit1++)
	    {
    		pn2 = &(mit1->second);
			vpar.push_back((*pn2)); // if alternative 
		}
  // aggregate the parcel data		
	vpa.push_back(for_each(vpar.begin(),vpar.end(),A()));
	vpar.clear();
	sit1++;
 }
// get the range of objects with a -1 key if present in multimap collection
		ip = -1;
 		pair<mit, mit> pastrange = m1.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		if (j>1) {
			for (mit1=pastrange.first; mit1!=pastrange.second;mit1++)
			{
    			pn2 = &(mit1->second);
				vpar.push_back((*pn2)); // if alternative 
			}
		  // aggregate the collection data		
		    vpa.push_back(for_each(vpar.begin(),vpar.end(),A()));
			vpar.clear();
		}
 return vpa;
}


template < typename H,typename L, typename M,typename N, typename A> 
deque<A> vertaggbyMaxKey (H& h1, L& l1, M& m1, N& n1,A& a1)
{
 N* pn2;
 typedef M::iterator mit;
 mit mit1;

 deque<N> vpar;
 deque<A> vpa;
 deque<A>::iterator vpait;
// aggregate vertex by a key value
///*
 for (int i=h1; i>= (int) l1;i--) {
// get the range of parcel objects of the current stop id
 		pair<mit, mit> pastrange = m1.equal_range(i);
	   	size_t j = distance(pastrange.first,pastrange.second);
		if (j>0) 
		{
			for (mit1=pastrange.first; mit1!=pastrange.second;mit1++)
			{
    			pn2 = &(mit1->second);
				vpar.push_back((*pn2)); // if alternative 
			}
  // aggregate the parcel data		
			vpa.push_back(for_each(vpar.begin(),vpar.end(),A()));
			vpar.clear();
		}
 }
return vpa;
}


#endif // AGGREGATEBYKEY_H ///:~

// geenrate the dp pattern
#ifndef DPATTERNGENERATOR_H
#define DPATTERNGENERATOR_H
// m0 - multimap, o0 - maximum skip value  
// n0 - binary mapdata input filename , f0 - pattern output filename 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGenerator(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;
m0::iterator m0mit;
tstop stop0;
//tstop* stop1=&stop0;
i=-1;
fileName(n1,"dpat.txt",f1,o1,i);
ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"q"<<"\t"<<"l"<<"\t"<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"m"<<"\t"<<"DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
fileName(n1,".bin",f1,i,i);
ifstream infile(f1, ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	for (i=1;i<=o1;i++)
	{
		for (j=i;j<=o1;j++)
		{
			if (i==j) 
			{k=i;}
			else 
			{k=i+j;} // k - pattern repetition for a specific i,j combination
			// make stop pattern for this set skip1=i,skip2=j starting at i from starting point m from 1 to k
			for (m=1;m<=k;m++) 
			{
				n = m0m.size();
	
// remove the stops between the begining stop and i distance away, and then all the way to m
// if m > i and < i+j remove as much as can be upto j starting at + i 
				o=2;
				for (p=m-1;p>1 && p-j!=0 && p<m;p--) 
				{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
				}
				o=p;
				for (p=o-j-1;p>1 && p<m;p--) 
				{ //i - remove equal number of stops
					m0mit = m0m.find(p);
					m0m.erase(m0mit);
				}
				o=m+1;

				while (o<n)
				{
					for (p=o;p<(o+i-1) && p<n;p++) 
					{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
					}
					o=p+1;
					for (p=o;p<(o+j-1) && p<n;p++) 
					{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
					}
					o=p+1;
				}
				q++;
				outfile<<q<<"\t"<<l<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<m<<"\t";
				string txthdr = ""; 
				writeTextDPatternData(m0m,stop0,outfile,txthdr);
				// open and read stop data into multimap
				m0m.clear();

				infile.open(f1, ios::in | ios::binary);
				readBinStopObjectData(m0m,stop0,infile);
				if (infile.is_open()) {
					infile.close();
				}
				infile.clear();
			}
		}
	}
return m0m;
}

// o1 - the maximum difference in stops between between successive stops 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGenerator3d(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;
m0::iterator m0mit;
tstop stop0;
//tstop* stop1=&stop0;
i=-1;
f1 = n1  + "_" + to_string<long>(o1) +"Mdpat5d.tab";

//fileName(n1,"Mdpat3d.tab",f1,o1,i);
ofstream outfile(f1.c_str(), ios::out |ios::trunc);

//ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"ScId"<<"\t"<<"i-j"<<"\t"<<"j-k"<<"\t"<<"MxPat"<<"\t"<<"PatNo"<<"\t"<<" DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
//fileName(n1,".bin",f1,i,i);
f1 = n1  + ".bin";
ifstream infile(f1.c_str(), ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	for (i=1;i<=o1;i++)
	{
		for (j=i;j<=o1;j++)
		{
			if (i==j) 
			{k=i;}
			else 
			{k=i+j;} // k - pattern repetition for a specific i,j,k combination
			// make stop pattern for this set skip1=i-1,skip2=j-1 from 1 to k starting at m 
			for (m=1;m<=k;m++) 
			{
				n = m0m.size();
	
// remove the stops between the begining stop and i distance away, and then all the way to m
// if m > i and < i+j remove as much as can be upto j starting at + i 
				o=2;
				for (p=m-1;p>1 && p>m-j;p--) 
				{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
				}
				o=p;
				for (p=m-j-1;p>1 && m-j-i>1;p--) 
				{ //i - remove equal number of stops
					m0mit = m0m.find(p);
					m0m.erase(m0mit);
				}
				o=m+1;

				while (o<n)
				{
					for (p=o;p<(o+i-1) && p<n;p++) 
					{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
					}
					o=p+1;
					for (p=o;p<(o+j-1) && p<n;p++) 
					{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
					}
					o=p+1;
				}
				q++;
				outfile<<q<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<m<<"\t";
				string txthdr = ""; 
				writeTextDPatternData(m0m,stop0,outfile,txthdr);
				// open and read stop data into multimap
				m0m.clear();

				infile.open(f1.c_str(), ios::in | ios::binary);
				readBinStopObjectData(m0m,stop0,infile);
				if (infile.is_open()) {
					infile.close();
				}
				infile.clear();
			}
		}
	}
return m0m;
}

// o1 - the maximum difference in stops between between successive stops 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGenerator3dcut(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;
m0 m1m;
m0::iterator m0mit;
tstop stop0;
typedef pair <long, tstop> o_Pair; 

//tstop* stop1=&stop0;
i=-1;
fileName(n1,"Mdpat3dcut.tab",f1,o1,i);
ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"ScId"<<"\t"<<"i-j"<<"\t"<<"j-k"<<"\t"<<"MxPat"<<"\t"<<"PatNo"<<"\t"<<" DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
fileName(n1,".bin",f1,i,i);
ifstream infile(f1, ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	 // 3d (i,j,k) then make DP pattern 
	for (i=o1+1;i>=1;i--)
	{
		for (j=i+1;j<=i+o1 && j<=n;j++)
		{
			for (k=j+1;k<=j+o1 && k<=n;k++)
			{
						q++; //pattern serial id
						// first insert the i,j,k,l,m stops into the new set - m1m 
						m0mit = m0m.find(i);
						if (m0mit!=m0m.end()) {
							m1m.insert(o_Pair(i,m0mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<i;
						}
						m0mit = m0m.find(j);
						if (m0mit!=m0m.end()) {
							m1m.insert(o_Pair(j,m0mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<j;
						}
						m0mit = m0m.find(k);
						if (m0mit!=m0m.end()) {
							m1m.insert(o_Pair(k,m0mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<k;
						}
						// move forward until m = n according to the pattern by
						// repeating the i,j,k pattern itself  
						p=k;
						while (o<n)
						{
							o=o+j-i;
							if (o>=n) {
								m0mit = m0m.find(n);
								m1m.insert(o_Pair (n,m0mit->second));
								break;
							}
							m0mit = m0m.find(o);
							if (m0mit!=m0m.end()) {
								m1m.insert(o_Pair (o,m0mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o+k-j;
							if (o>=n) {
								m0mit = m0m.find(n);
								m1m.insert(o_Pair (n,m0mit->second));
								break;
							}
							m0mit = m0m.find(o);
							if (m0mit!=m0m.end()) {
								m1m.insert(o_Pair (o,m0mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
						}
						o=p=i;
						while (o>0)
						{
							o=o-k+j;
							if (o<=1) {
								p=1;
								m0mit = m0m.find(p);
								m1m.insert(o_Pair (p,m0mit->second));
								break;}
							m0mit = m0m.find(o);
							if (m0mit!=m0m.end()) {
								m1m.insert(o_Pair (o,m0mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o-j+i;
							if (o<=1) {
								p=1;
								m0mit = m0m.find(p);
								m1m.insert(o_Pair (p,m0mit->second));
								break;
							}
							m0mit = m0m.find(o);
							if (m0mit!=m0m.end()) {
								m1m.insert(o_Pair (o,m0mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
						}

				outfile<<q<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<m<<"\t";
				string txthdr = ""; 
				writeTextDPatternData(m0m,stop0,outfile,txthdr);
				// open and read stop data into multimap
				m0m.clear();

				infile.open(f1, ios::in | ios::binary);
				readBinStopObjectData(m0m,stop0,infile);
				if (infile.is_open()) {
					infile.close();
				}
				infile.clear();
			}
		}
	}
return m0m;
}


// o1 - the maximum difference in stops between between successive stops 
// n1 - the name of the stop file 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGenerator4d(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;
m0::iterator mit;
tstop stop0;
//tstop* stop1=&stop0;
i=-1;
fileName(n1,"Mdpat4d.tab",f1,o1,i);
ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"ScId"<<"\t"<<"i-j"<<"\t"<<"j-k"<<"\t"<<"k-l"<<"\t"<<"MxPat"<<"\t"<<"PatNo"<<"\t"<<"DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
fileName(n1,".bin",f1,i,i);
ifstream infile(f1, ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	for (i=1;i<=o1;i++)
	{
		for (j=i;j<=o1;j++)
		{
			for (l=j;l<=o1;l++)
			{
					if (i==j && j==l ) {
						k=i;}
					else {
						k=i+j+l;} // k - pattern repetition for a specific i,j,l combination
				// make stop pattern for this set skip1=i-1,skip2=j-1,skip3=l-1 from 1 to k starting at m  
					for (r=1;r<=k;r++) 
					{
						n = (int) m0m.size();
	// remove the stops between the begining stop and i distance away, and then all the way to m
	// if m > i and < i+j remove as much as can be upto j starting at + i 
						o=r-1;
						for (p=o;p>1 && p>o+1-l;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-l-1;
						for (p=o;p>1 && p>o+1-j;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-l-j-1;
						for (p=o;p>1 && o-i>0;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r+1;

						while (o<n)
						{
							for (p=o;p<(o+i-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+j-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+l-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
						}
						q++;
						outfile<<q<<"\t"<<i<<"\t"<<j<<"\t"<<l<<"\t"<<k<<"\t"<<r<<"\t";
						string txthdr = ""; 
						writeTextDPatternData(m0m,stop0,outfile,txthdr);
						// open and read stop data into multimap
						m0m.clear();
	
						infile.open(f1, ios::in | ios::binary);
						readBinStopObjectData(m0m,stop0,infile);
						if (infile.is_open()) {
								infile.close();
						}
						infile.clear();
					} // r 1..k (where k-is the pattern length)
			} // l
		} // j
	} // i
return m0m;
}

// o1 - the maximum difference in stops between between successive stops 
// n1 - the name of the stop file 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGenerator5dold(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;
m0::iterator mit;
tstop stop0;
//tstop* stop1=&stop0;
i=-1;
fileName(n1,"Mdpat5d.tab",f1,o1,i);
ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"ScId"<<"\t"<<"i-j"<<"\t"<<"j-k"<<"\t"<<"k-l"<<"\t"<<"l-m"<<"\t"<<"MxPat"<<"\t"<<"PatNo"<<"\t"<<" DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
fileName(n1,".bin",f1,i,i);
ifstream infile(f1, ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	for (i=1;i<=o1;i++)
	{
		for (j=i;j<=o1;j++)
		{
			for (l=j;l<=o1;l++)
			{
				for (m=l;m<=o1;m++)
				{
					if (i==j && j==l && l==m) {
						k=i;}
					else {
						k=i+j+l+m;} // k - pattern repetition for a specific i,j,l,m combination
				// make stop pattern for this set skip1=i-1,skip2=j-1,skip3=l-1,skip4=m-1 from 1 to k starting at m  
					for (r=1;r<=k;r++) 
					{
						n = (int) m0m.size();
	// remove the stops between the begining stop and i distance away, and then all the way to m
	// if m > i and < i+j remove as much as can be upto j starting at + i 
						o=r-1;
						for (p=o;p>1 && p>o+1-m;p--) 
						{ //m - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-m-1;
						for (p=o;p>1 && p>o+1-l;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-m-l-1;
						for (p=o;p>1 && p>o+1-j;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-m-l-j-1;
						for (p=o;p>1 && o-i>0;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r+1;

						while (o<n)
						{
							for (p=o;p<(o+i-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+j-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+l-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+m-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
						}
						q++;
						outfile<<q<<"\t"<<i<<"\t"<<j<<"\t"<<l<<"\t"<<m<<"\t"<<k<<"\t"<<r<<"\t";
						string txthdr = ""; 
						writeTextDPatternData(m0m,stop0,outfile,txthdr);
						// open and read stop data into multimap
						m0m.clear();
	
						infile.open(f1, ios::in | ios::binary);
						readBinStopObjectData(m0m,stop0,infile);
						if (infile.is_open()) {
								infile.close();
						}
						infile.clear();
					} // r 1..k (where k-is the pattern length)
				} // m
			} // l
		} // j
	} // i
return m0m;
}


// o1 - the maximum difference in stops between between successive stops 
// n1 - the name of the stop file 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGeneratorxd(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;
m0::iterator mit;
tstop stop0;
//tstop* stop1=&stop0;
i=-1;
fileName(n1,"Mdpatxd.tab",f1,o1,i);
ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"ScId"<<"\t"<<"i-j"<<"\t"<<"j-k"<<"\t"<<"k-l"<<"\t"<<"l-m"<<"\t"<<"MxPat"<<"\t"<<"PatNo"<<"\t"<<" DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
fileName(n1,".bin",f1,i,i);
ifstream infile(f1, ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	for (i=1;i<=o1;i++)
	{
		for (j=1;j<=o1;j++)
		{
			for (l=1;l<=o1;l++)
			{
				for (m=l;m<=o1;m++)
				{
					if (i==j && j==l && l==m) {
						k=i;}
					else {
						k=i+j+l+m;} // k - pattern repetition for a specific i,j,l,m combination
				// make stop pattern for this set skip1=i-1,skip2=j-1,skip3=l-1,skip4=m-1 from 1 to k starting at m  
					for (r=1;r<=k;r++) 
					{
						n = (int) m0m.size();
	// remove the stops between the begining stop and i distance away, and then all the way to m
	// if m > i and < i+j remove as much as can be upto j starting at + i 
						o=r-1;
						for (p=o;p>1 && p>o+1-m;p--) 
						{ //m - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-m-1;
						for (p=o;p>1 && p>o+1-l;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-m-l-1;
						for (p=o;p>1 && p>o+1-j;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r-m-l-j-1;
						for (p=o;p>1 && o-i>0;p--) 
						{ //i - remove equal number of stops
							mit = m0m.find(p);
							m0m.erase(mit);
						}
						o=r+1;

						while (o<n)
						{
							for (p=o;p<(o+i-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+j-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+l-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
							for (p=o;p<(o+m-1) && p<n;p++) 
							{ //j - remove equal number of stops
								mit = m0m.find(p);
								m0m.erase(mit);
							}
							o=p+1;
						}
						q++;
						outfile<<q<<"\t"<<i<<"\t"<<j<<"\t"<<l<<"\t"<<m<<"\t"<<k<<"\t"<<r<<"\t";
						string txthdr = ""; 
						writeTextDPatternData(m0m,stop0,outfile,txthdr);
						// open and read stop data into multimap
						m0m.clear();
	
						infile.open(f1, ios::in | ios::binary);
						readBinStopObjectData(m0m,stop0,infile);
						if (infile.is_open()) {
								infile.close();
						}
						infile.clear();
					} // r 1..k (where k-is the pattern length)
				} // m
			} // l
		} // j
	} // i
return m0m;
}


// o1 - the maximum difference in stops between between successive stops 
// n1 - the name of the stop file 
template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGenerator5d(m0& m1, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,z=o1+1;
m0::iterator mit;
m0 m2;

typedef pair <long, tstop> mPair;

tstop stop0;
//tstop* stop1=&stop0;
i=-1;
f1 = n1  + "_" + to_string<long>(o1) +"Mdpat5d.tab";

//fileName(n1,"Mdpat5d.tab",f1,o1,i);
ofstream outfile(f1.c_str(), ios::out |ios::trunc);
outfile<<"ScId"<<"\t"<<"i-j"<<"\t"<<"j-k"<<"\t"<<"k-l"<<"\t"<<"l-m"<<"\t"<<"MxPat"<<"\t"<<"PatNo"<<"\t"<<" DPattern"<<endl;

// open and read stop data into multimap
//m1.clear();
//fileName(n1,".bin",f1,i,i);
//ifstream infile(f1, ios::in | ios::binary);
//readBinStopObjectData(m1,stop0,infile);
//if (infile.is_open()) {
//	infile.close();
// }
//infile.clear();
	// pattern making 
	n = (int) m1.size();
	// to generate the arcs for the min-cut start with the boundary stage at 
	// z= o1+1 , where o1 is equal to maximum stops skipped plus one (o1=s+1) 
	// and move backwards until the begining of the stop. Each quintuplet formulated
	// defines a pattern that will be built into a scenario for the whole stop set.
	// since the arcs are unique (part of being a min-cut that covers the arc set)
	// we have a proof that our generated senarios produce the complete immediate costs
	for (i=z;i>=1;i--)
	{
		// start always with z since it is the cut boundary stage
		for (j=z+1;j<=i+o1 && j<=n;j++)
		{
			for (k=j+1;k<=j+o1 && k<=n;k++)
			{
				for (l=k+1;l<=k+o1 && l<=n;l++)
				{
					for (m=l+1;m<=l+o1 && m<=n;m++) 
					{
						q++;
						// first insert the i,j,k,l,m stops into the new set - m2 
						mit = m1.find(i);
						if (mit!=m1.end()) {
							m2.insert(mPair(i,mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<i;
						}
						mit = m1.find(j);
						if (mit!=m1.end()) {
							m2.insert(mPair(j,mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<j;
						}
						mit = m1.find(k);
						if (mit!=m1.end()) {
							m2.insert(mPair(k,mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<k;
						}
						mit = m1.find(l);
						if (mit!=m1.end()) {
							m2.insert(mPair(l,mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<l;
						}
						mit = m1.find(m);
						if (mit!=m1.end()) {
							m2.insert(mPair(m,mit->second));
						} else {
							outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<m;
						}
						// move forward until m = n according to the pattern by
						// repeating the i,j,k,l,m pattern itself  
						o=p=m;
						while (o<n)
						{
							o=o+j-i;
							if (o>=n) {
								mit = m1.find(n);
								m2.insert(mPair (n,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o+k-j;
							if (o>=n) {
								mit = m1.find(n);
								m2.insert(mPair (n,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o+l-k;
							if (o>=n) {
								mit = m1.find(n);
								m2.insert(mPair (n,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o+m-l;
							if (o>=n) {
								mit = m1.find(n);
								m2.insert(mPair (n,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
						}
						o=p=i;
						while (o>0)
						{
							o=o-m+l;
							if (o<=1) {
								p=1;
								mit = m1.find(p);
								m2.insert(mPair (p,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o-l+k;
							if (o<=1) {
								p=1;
								mit = m1.find(p);
								m2.insert(mPair (p,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o-k+j;
							if (o<=1) {
								p=1;
								mit = m1.find(p);
								m2.insert(mPair (p,mit->second));
								break;}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
							o=o-j+i;
							if (o<=1) {
								p=1;
								mit = m1.find(p);
								m2.insert(mPair (p,mit->second));
								break;
							}
							mit = m1.find(o);
							if (mit!=m1.end()) {
								m2.insert(mPair (o,mit->second));
							} else {
								outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
							}
						}
						outfile<<q<<"\t"<<j-i<<"\t"<<k-j<<"\t"<<l-k<<"\t"<<m-l<<"\t"<<p<<"\t"<<o<<"\t";
						string txthdr = ""; 
						writeTextDPatternData(m2,stop0,outfile,txthdr);
						m2.clear();
					} // r 1..k (where k-is the pattern length)
				} // m
			} // l
		} // j
	} // i
return m2;
}


template <typename m0 , typename o0, typename n0, typename f0>
m0& PatternGeneratorq(m0& m0m, o0& o1, n0& n1, f0& f1)
{
int i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;
m0::iterator m0mit;
tstop stop0;
//tstop* stop1=&stop0;
i=-1;
fileName(n1,"dpatq.txt",f1,o1,i);
ofstream outfile(f1, ios::out |ios::trunc);
outfile<<"q"<<"\t"<<"l"<<"\t"<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"m"<<"\t"<<"DPattern"<<endl;

// open and read stop data into multimap
m0m.clear();
fileName(n1,".bin",f1,i,i);
ifstream infile(f1, ios::in | ios::binary);
readBinStopObjectData(m0m,stop0,infile);
if (infile.is_open()) {
	infile.close();
 }
infile.clear();
	// pattern making 
	for (i=1;i<=o1;i++)
	{
		for (j=i;j<=o1;j++)
		{
			if (i==j) 
			{k=i;}
			else 
			{k=i+j;} // k - pattern repetition for a specific i,j combination
			// make stop pattern for this set skip1=i,skip2=j starting at i from starting point m from 1 to k
			for (m=1;m<=k;m++) 
			{
				n = m0m.size();
	
// remove the stops between the begining stop and i distance away, and then all the way to m
// if m > i and < i+j remove as much as can be upto j starting at + i 
				o=2;
				for (p=m-1;p>1 && p>m-j;p--) 
				{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
				}
				o=p;
				for (p=m-j-1;p>1 && m-j-i>1;p--) 
				{ //i - remove equal number of stops
					m0mit = m0m.find(p);
					m0m.erase(m0mit);
				}
				o=m+1;

				while (o<n)
				{
					for (p=o;p<(o+i-1) && p<n;p++) 
					{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
					}
					o=p+1;
					for (p=o;p<(o+j-1) && p<n;p++) 
					{ //j - remove equal number of stops
						m0mit = m0m.find(p);
						m0m.erase(m0mit);
					}
					o=p+1;
				}
				q++;
				outfile<<q<<"\t"<<l<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<m<<"\t";
				string txthdr = ""; 
				writeTextDPatternData(m0m,stop0,outfile,txthdr);
				// open and read stop data into multimap
				m0m.clear();

				infile.open(f1, ios::in | ios::binary);
				readBinStopObjectData(m0m,stop0,infile);
				if (infile.is_open()) {
					infile.close();
				}
				infile.clear();
			}
		}
	}
return m0m;
}


#endif // DPATTERNGENERATOR_H ///:~


#ifndef DP_OPTIMIZATION_H
#define DP_OPTIMIZATION_H

// forward dp run - inputs multimap with dp objects,  
template <typename m, typename o,typename n, typename  M,typename s>
m& dpforun(m& m1, o& o1,n& n1, M& M1, s& txt)
{
m::iterator mit;
 mit = m1.begin();
long i=0,j=0,k=0,l=0;
 //at stop 1 , f*(0) = 0
	double icost=0, fcost,fcostar;
//	mapDPRes.insert(lngdbl_Pair(1,0.0)); 
//	mmapStrDPRes.insert(strdbl_Pair("-1112-1",0.0)); 
//	mmapStopDPRes.insert(lngstr_Pair(1,"-1112-1"));
// run dp routine
	for (j=1;j<(n1);j++)
	{
			k=j+1;
		for (i=j;i>j-M1 && i>0;i--)
		{
			
				strkey = "-1" + to_string(i)+ to_string(j)+ to_string(k)+ "-1";
				mit=m1.find(strkey);
				if (mit!=m1.end()) 
				{
					icost = mit->second.get_tstop().get_TCost();
					if (icost<=0) { icost=inf;}
					mapDPRit = mapDPRes.find(i);
					if ( mapDPRit!=mapDPRes.end() )
					{
						fcost = icost + mapDPRit->second;
					}
					else
					{
						fcost = icost;
					}
					mmapIJPResICost.insert(lngdbl_Pair(j,fcost));
					mmapKeyDPICost.insert(strdbl_Pair(strkey,fcost));
					mmapICostDPIJ.insert(dblstr_Pair(fcost,strkey));
				}
				else
				{
					cout<<"key " <<strkey<<" does not have an immediate cost"<<endl;
				}

		}
					mmapICDPIJRit=mmapICostDPIJ.begin();
					fcostar =mmapICDPIJRit->first;
					strkey =mmapICDPIJRit->second;

					mapDPRes.insert(lngdbl_Pair(j,fcostar)); 
					mmapStrDPRes.insert(strdbl_Pair(strkey,fcostar)); 
					mmapStopDPRes.insert(lngstr_Pair(j,strkey));
					mmapICostDPIJ.clear();
					mmapIJPResICost.clear();
					mmapKeyDPICost.clear();

	}


 q=-1;
 fileName(tstopiname,"dpts.txt",outfilename,q);
 outfile.open(outfilename, ios::out  );

				{// write the dp result data to a file
					writeTextData(mapDPRes,j,fcost,outfile,"Stage|OptimalCost");
					writeTextData(mmapStrDPRes,strkey,fcost,outfile,"KeyPath|OptimalCost");
					writeTextData(mmapStopDPRes,j,strkey,outfile,"Stage|KeyPath");
				}

 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 }


// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename m, typename o,typename n, typename  M,typename s,typename f>
m& dpbarun(m& m1, o& o1,n& n1, M& M1, s& s1,f& f1)
{
 m::iterator mit;
 mit = m1.begin();
 long i=0,j=0,k=0,l=0;
 string strijk, strij,strjk;
 //at stop 1 , f*(0) = 0
	double icost=0, fcost,fcostar;
//	mapDPRes.insert(lngdbl_Pair(1,0.0)); 
//	mmapStrDPRes.insert(strdbl_Pair("-1112-1",0.0)); 
//	mmapStopDPRes.insert(lngstr_Pair(1,"-1112-1"));
// run dp routine
	strjk = to_string(n1) + "_" + to_string(n1);
	mmapIJDPK.insert(strlng_Pair(strjk,n1)); // predecessor  
	mmapIJDPRes.insert(strdbl_Pair(strjk,fcostar)); 
	mmapStopDPRes.insert(lngstr_Pair(n1,strjk));
	mapDPRes.insert(lngdbl_Pair(n1,fcostar));


	for (j=n1;j>=(1);j--)
	{

		for (i=j;i>=j-M1 && i>0;i--)
		{
				strij = to_string(i) + "_" + to_string(j);
			for (k=j;k<j+M1 && k<=n1;k++)
			{
				strjk = to_string(j) + "_" + to_string(k);
				strijk = "-1" + to_string(i)+ to_string(j)+ to_string(k)+ "-1";
				mit=m1.find(strijk);
				if (mit!=m1.end()) 
				{
					icost = mit->second.get_tstop().get_TCost();
					if (icost<0) { icost=inf;}
					mmapIJKDPRit = mmapIJDPRes.find(strjk);
					mapDPRit = mapDPRes.find(k);
					if ( mmapIJKDPRit!=mmapIJDPRes.end() )
					{
						fcost = icost + mmapIJKDPRit->second;
					}
					else
					{
						fcost = icost;
					}
					mmapIJKDPICost.insert(strdbl_Pair(strij,fcost));
					mmapICostDPIJ.insert(dblstr_Pair(fcost,strjk));
					mmapDPIC.insert(dblng_Pair(fcost,k));
					mmapIJDPResICost.insert(strlng_Pair(strij,k));
				}
				else
				{
					cout<<"key " <<strkey<<" does not have an immediate cost"<<endl;
				}
			} // k loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapDPICit!=mmapDPIC.end() )
					{
						fcostar =mmapDPICit->first;
						k =mmapDPICit->second;
						strjk =mmapICDPIJRit->second;
							mmapIJDPK.insert(strlng_Pair(strij,k)); // predecessor  
							mmapStopDPRes.insert(lngstr_Pair(k,strij));
							mmapDPredK.insert(lng_Pair(j,k));

						mmapIJDPRes.insert(strdbl_Pair(strij,fcostar)); 
						mapDPRes.insert(lngdbl_Pair(j,fcostar));
					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapIJKDPICost.clear();
						mmapDPIC.clear();
		}// i loop

	} // j loop


// add the cost at the start of the route
	i=1; 
 // trace the optimal path from the above result

 q=-1;
 fileName(tstopiname,"dptresult.txt",outfilename,q);
 outfile.open(outfilename, ios::out  );

				{// write the dp result data to a file
					writeTextData(mmapIJDPK,strkey,j,outfile,"Stage\tK");
					writeTextData(mmapIJDPRes,strkey,fcost,outfile,"KeyPath\tOptimalCost");
					writeTextData(mmapStopDPRes,j,strkey,outfile,"Stage\tKeyPath");
					writeTextData(mmapIJDPK,strkey,j,outfile,"i-j\tk");
					writeTextData(mmapDPredK,j,i,outfile,"i\tj");
				}
			


 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 }

// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename m, typename o,typename n, typename  M,typename s,typename f,typename t>
m& dpbarun2(m& m1, o& o1,n& n1, M& M1, s& s1,f& f1,t& t1)
{
 m::iterator mit;
 typedef map<s, o> mDP;
 mDP mDP1;
 typedef pair <s, o> o_Pair; 
 mit = m1.begin();
 long i=0,j=0,k=0,l=0,m=0;
 char outfilename[ MaxStrLen +1]=""; // data output file 
 string strijk, strij,strjk;
 //at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
// run dp routine
	strjk = to_string(n1) + "_" + to_string(n1);
	mmapIJDPK.insert(strlng_Pair(strjk,n1)); // predecessor  
	mmapIJDPRes.insert(strdbl_Pair(strjk,fcostar)); 
	mmapStopDPRes.insert(lngstr_Pair(n1,strjk));
	mapDPRes.insert(lngdbl_Pair(n1,fcostar));
 //mDP.clear();
	l=-1;
	m=-1;
	for (j=n1;j>=(1);j--)
	{

		for (i=j;i>=j-M1 && i>0;i--)
		{
				strij = to_string(i) + "_" + to_string(j);
			for (k=j;k<j+M1 && k<=n1;k++)
			{
				strjk = to_string(j) + "_" + to_string(k);
				strijk = to_string(l) + *t1 + to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1+ to_string(m);
				mit=m1.find(strijk);
				if (mit!=m1.end()) 
				{
					o1 = mit->second;
					icost = o1.get_TCost();
					if (icost<0) { icost=inf;}
					mmapIJKDPRit = mmapIJDPRes.find(strjk);
					mapDPRit = mapDPRes.find(k);
					if ( mmapIJKDPRit!=mmapIJDPRes.end() )
					{
						fcost = icost + mmapIJKDPRit->second;
					}
					else
					{
						fcost = icost + inf;
					}
					mmapICostDPIJ.insert(dblstr_Pair(fcost,strjk));
					mmapDPIC.insert(dblng_Pair(fcost,k));
					mmapIJDPResICost.insert(strlng_Pair(strij,k));
				}
				else
				{
					cout<<"key " <<strijk<<" does not have an immediate cost"<<endl;
				}
			} // k loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapDPICit!=mmapDPIC.end() )
					{
						fcostar =mmapDPICit->first;
						k =mmapDPICit->second;
						strjk =mmapICDPIJRit->second;
							mmapIJDPK.insert(strlng_Pair(strij,k)); // predecessor  
							mmapStopDPRes.insert(lngstr_Pair(k,strij));
							mmapDPredK.insert(lng_Pair(j,k));

						mmapIJDPRes.insert(strdbl_Pair(strij,fcostar)); 
						mapDPRes.insert(lngdbl_Pair(j,fcostar));
					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
		}// i loop

	} // j loop


 // trace the optimal path from the above result
 // find the i1=1_j1=1 (i.e. "1_1") key value of k1 value then find the next i2=j1_j2=k1 value k2 recurse  
	i=1; 
	j=1;
	size_t pos;
	while (k<n1) 
	{
		strij = to_string(i) + "_" + to_string(j);
		mmapIJDPRit = mmapIJDPK.find(strij);
		k=mmapIJDPRit->second;
		mmapIJDPResICost.insert(strlng_Pair(strij,k));
	    pos = strij.find("_");    // position of "live" in str
		strjk = strij.substr(pos+1); //get from "_" to the end
		strijk = to_string(l) + *t1 + to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1+ to_string(m);
		mit = m1.find(strijk);
		if (mit!=m1.end()) {
			o1 = mit->second;
			mDP1.insert (o_Pair(strij,o1));
		}
		i=  from_String<long>(strjk);
		j=k;
	}
// write the resultant dp path
 int q=-1;
 fileName(f1,"dpobjectrace.txt",outfilename,q,-1);
 outfile.open( outfilename, ios::out );
				{// write the dp result data to a file
	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
					writeTextObjectData(mDP1,o1,strijk,outfile,txthdr);
				}
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 fileName(f1,"dptrace.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );
				{// write the dp result data to a file
					writeTextData(mmapIJDPResICost,strijk,j,outfile,"Stij\tK");
				}
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
// add the cost at the start of the route
 fileName(f1,"dptres.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );

				{// write the dp result data to a file
					writeTextData(mmapIJDPK,strijk,j,outfile,"Stage\tK");
					writeTextData(mmapIJDPRes,strijk,fcost,outfile,"KeyPath\tOptimalCost");
					writeTextData(mmapStopDPRes,j,strijk,outfile,"Stage\tKeyPath");
					writeTextData(mmapIJDPK,strijk,j,outfile,"i-j\tk");
					writeTextData(mmapDPredK,j,i,outfile,"i\tj");
				}
			


 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
return m1;
 }

// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename m, typename o,typename n, typename  M,typename s,typename f,typename t>
m& dpbarun3(m& m1, o& o1,n& n1, M& M1, s& s1,f& f1,t& t1)
{
 m::iterator mit;
 typedef map<s, o> mDP;
 mDP mDP1;
 tstop ts;
 typedef pair <s, o> o_Pair; 
 typedef pair <s,tstop> tdPair;
 mit = m1.begin();
 long i=0,j=0,k=0,l=0,m=0;
 char outfilename[ MaxStrLen +1]=""; // data output file 
 string strijk, strij,strjk;
 //at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
// run dp routine
	strjk = to_string(n1) + *t1 + to_string(n1);
	mmapIJDPK.insert(strlng_Pair(strjk,(long) n1)); // predecessor  
	strjk = to_string(n1) + *t1 + to_string(n1);
	mmapIJDPRes.insert(strdbl_Pair(strjk,fcostar)); 
	mmapIJDPK.insert(strlng_Pair(strjk,(long) n1)); // predecessor  
	mmapStopDPRes.insert(lngstr_Pair((long) n1,strjk));
	mapDPRes.insert(lngdbl_Pair((long) n1,fcostar));
 //mDP.clear();
	l=-1;
	m=-1;
	for (i= (long) n1;i>=(1);i--)
	{

		for (j=i;j>=i-M1 && j>0;j--)
		{
				strij = to_string(i) + *t1 + to_string(j);
			for ( k = j ; k < (j + (long) M1) && k <= (long) n1; k++)
			{
				strjk = to_string(j) + *t1 + to_string(k);
				strijk = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) ;
				mit=m1.find(strijk);
				if (mit!=m1.end()) 
				{
					o1 = mit->second;
					icost = mit->second.get_tstop().get_TCost();
					if (icost<0) { icost=inf;}
					mmapIJKDPRit = mmapIJDPRes.find(strjk);
					mapDPRit = mapDPRes.find(k);
					if ( mmapIJKDPRit!=mmapIJDPRes.end() )
					{
						fcost = icost + mmapIJKDPRit->second;
					}
					else
					{
						fcost = icost + inf;
					}
					mmapICostDPIJ.insert(dblstr_Pair(fcost,strjk));
					mmapDPIC.insert(dblng_Pair(fcost,k));
					mmapIJDPResICost.insert(strlng_Pair(strij,k));
				}
				else
				{
					cout<<"key " <<strijk<<" does not have an immediate cost"<<endl;
				}
			} // k loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapDPICit!=mmapDPIC.end() )
					{
						fcostar =mmapDPICit->first;
						k =mmapDPICit->second;
						strjk =mmapICDPIJRit->second;
							mmapIJDPK.insert(strlng_Pair(strij,k)); // predecessor  
							mmapStopDPRes.insert(lngstr_Pair(k,strij));
							mmapDPredK.insert(lng_Pair(j,k));

						mmapIJDPRes.insert(strdbl_Pair(strij,fcostar)); 
						mapDPRes.insert(lngdbl_Pair(j,fcostar));
					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
		}// j loop

	} // i loop


 // trace the optimal path from the above result
// write the resultant dp path
 int q=-1;
 fileName(f1,"dpobjectrace.txt",outfilename,q,-1);
 outfile.open( outfilename, ios::out );
	// write the dp result data to a file
	o1.serializetexthdr(outfile);
// find the i1=1 + *t1 + j1=1 (i.e. "1|1") key value of k1 value then find the next i2=j1 + *t1 + j2=k1 value k2 recurse  
	i=1; 
	j=1;
	size_t pos=0;
	while (k < (long) n1) 
	{
		strij = to_string(i) + *t1 + to_string(j);
		mmapIJDPRit = mmapIJDPK.find(strij);
		k=mmapIJDPRit->second;
		mmapIJDPResICost.insert(strlng_Pair(strij,k));
	    pos = strij.find(*t1);    // position of "live" in str
		strjk = strij.substr(pos+1); //get from "_" to the end
		strijk = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k);
		mit = m1.find(strijk);
		if (mit!=m1.end()) {
			mit->second.serializetext(outfile);
			o1 = mit->second;
			ts = o1.get_tstop();
			mDP1.insert (o_Pair(strij,o1));
//			mStop.insert (tdPair(strij,ts));
		}
		i=  fromString<long>(strjk);
		j=k;
	}
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 fileName(f1,"dptrace.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );
				{// write the dp result data to a file
					writeTextData(mmapIJDPResICost,strijk,j,outfile,"Stij\tK");
				}
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
// add the cost at the start of the route
 fileName(f1,"dptres.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );

				{// write the dp result data to a file
					writeTextData(mmapIJDPK,strijk,j,outfile,"Stage\tK");
					writeTextData(mmapIJDPRes,strijk,fcost,outfile,"KeyPath\tOptimalCost");
					writeTextData(mmapStopDPRes,j,strijk,outfile,"Stage\tKeyPath");
					writeTextData(mmapIJDPK,strijk,j,outfile,"i-j\tk");
					writeTextData(mmapDPredK,j,i,outfile,"i\tj");
				}
			


 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
return m1;
 }

// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename m, typename o,typename n, typename  M,typename s,typename f,typename t>
m& dpbarun5dx(m& m1, o& o1,n& n1, M& M1, s& s1,f& f1,t& t1)
{
	m::iterator mit;
	typedef map<s, o> mDP;
	mDP mDP1;
	tstop ts;
	typedef pair <s, o> o_Pair; 
	typedef pair <s,tstop> tdPair;
	mapstop mStop;
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
	char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijklm, strijkl,strjklm,strjkl,strklm,strijk,str1;
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
// run dp routine
	strjklm =  to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
	strjkl =  to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
	mmapIJDPK.insert(strlng_Pair(strjkl,n1)); // predecessor  
	mmapIJDPRes.insert(strdbl_Pair(strjkl,fcostar)); 
	mmapStopDPRes.insert(lngstr_Pair(n1,strjkl));
	mapDPRes.insert(lngdbl_Pair(n1,fcostar));

	int q=-1;
	fileName(f1,"dpmissed.txt",outfilename,q,-1);
	ofstream missfile( outfilename, ios::out | ios::trunc );

	//mDP.clear();
	for (k=n1;k>0;k--)
	{
		for (j=k;j>=k-M1 && j>0;j--)
		{
			for (i=j;i>=j-M1 && i>0;i--)
			{
						strijk = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) ;
				for (l=k;l<=k+M1 && l<=n1;l++)
				{
						strijkl = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l);
						strjkl = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l);

					for (m=l;m<=l+M1 && m<=n1;m++)
					{
						strjklm = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1+ to_string(m);
						strijklm = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1+ to_string(m);
						mit=m1.find(strijklm);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_tstop().get_TCost();
							if (icost<0) { icost=inf;}
							mmapKLMDPRit = mmapIJDPRes.find(strjklm);
							mapDPRit = mapDPRes.find(m);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (m==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strijkl));
							mmapDPIC.insert(dblng_Pair(fcost,m));
							mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
						}
						else
						{	
							if (i==j||j==k||k==l||l==m) {
								continue;
							} else {
								cout<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
							}
						}
					} // m - loop
				} // l - loop
			} // i loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						m =mmapDPICit->second;
						strijkl =mmapICDPIJRit->second;
						strijk = strijkl.substr(strijkl.find_last_of(*t1)-1);
						mmapIJDPK.insert(strlng_Pair(strijkl,m)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(m,strijkl));
						mmapDPredK.insert(lng_Pair(k,m));

						mmapIJDPRes.insert(strdbl_Pair(strijkl,fcostar)); 
						mapDPRes.insert(lngdbl_Pair(k,fcostar));
					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
		}// j loop

	} // k loop

 if (missfile.is_open()) {
	 missfile.close();
 }
 missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
 fileName(f1,"dpobjtrc5d.txt",outfilename,q,-1);
 outfile.open( outfilename, ios::out );
 fileName(f1,"dptrc.txt",outfilename,q,-1);
 missfile.open(outfilename, ios::out  );
				// write the dp result data to a file
	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
//	outfile<<txthdr<<endl;
	o1.serializetexthdr(outfile);
//					writeTextObjectData(mDP1,o1,strijklm,outfile,txthdr);
//					writeTextObjectData(mStop,ts,strijklm,outfile,txthdr);

 // find the i1=1_j1=1 (i.e. "1_1") key value of k1 value then find the next i2=j1_j2=k1 value k2 recurse  
	i=1; 
	j=1;
	k=1;
	l=2;
	m=3;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
	while (mmapIJDPKit!=mmapIJDPK.end())
	{
		strijkl =mmapIJDPKit->first;
		strijk = strijkl.substr(0,strijkl.find_last_of(*t1));
		if (strijk== (to_string(i) + *t1 + to_string(j)+ *t1 + to_string(k))) {
			break;
		}
		mmapIJDPKit++;
	}

	mmapIJDPKit=mmapIJDPK.find(strijkl);
	if (mmapIJDPKit!=mmapIJDPK.end())
	{
			// find first *t1 
			strijk =strijkl; 
			i = getLoc(strijk,*t1,i);
//		i = from_String<long> strijkl.substr(1,i1-1);
		// get the substring beyound i
			strijk=getSubstr(strijk,*t1);
			// find next *t1 
			j = getLoc(strijk,*t1,j);
			// get the substring beyound j
			strijk=getSubstr(strijk,*t1);
			// find next *t1 
			k= getLoc(strijk,*t1,k);
			// get the substring beyound j
			strijk=getSubstr(strijk,*t1);
			// find last *t1 
			l= getLoc(strijk,*t1,l);		
		while (j<n1) 
		{
			mmapIJDPKit = mmapIJDPK.find(strijkl);
			if (mmapIJDPKit!=mmapIJDPK.end())
			{
				m=mmapIJDPKit->second;
				mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
			}
			else 
			{  if (m==n1)
				{
					strijkl =  to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;
					mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
				}

			}
			strijklm =  to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1+ to_string(m);
			mit = m1.find(strijklm);
			if (mit!=m1.end()) {
				mit->second.serializetext(outfile);
				o1.serializetext2(missfile);
				mit->second.serializetext(missfile);
				o1 = mit->second;
				ts = o1.get_tstop();
				mDP1.insert (o_Pair(strijkl,mit->second));
//			mStop.insert (tdPair(strijkl,ts));
			}
			i=j;
			j=k;
			k=l;
			l=m;
			strijkl =  to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;

		}
	}
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

/* 
	fileName(f1,"dptrc.txt",outfilename,q,-1);
	outfile.open(outfilename, ios::out  );
				{// write the dp result data to a file
					writeTextData(mmapIJDPResICost,strijklm,j,outfile,"Stij\tK");
				}

if (outfile.is_open()) {
	 outfile.close();
}
 outfile.clear();
*/
 fileName(f1,"dpPatrc.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );
	o1.serializetexthdr(outfile);
				{// write the dp result data to a file
					writeTextObjectData(mDP1,o1,strijklm,outfile,"");
				}
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 // add the cost at the start of the route
 fileName(f1,"dptrest.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );

				{// write the dp result data to a file
					writeTextData(mmapIJDPK,strijklm,j,outfile,"Stage\tK");
					writeTextData(mmapIJDPRes,strijklm,fcost,outfile,"KeyPath\tOptimalCost");
					writeTextData(mmapStopDPRes,j,strijklm,outfile,"Stage\tKeyPath");
					writeTextData(mmapIJDPK,strijklm,j,outfile,"i-j\tk");
					writeTextData(mmapDPredK,j,i,outfile,"i\tj");
				}
			


 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
return m1;
 }


// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename m, typename o,typename n, typename  M,typename s,typename f,typename t>
m& dpbarunxd(m& m1, o& o1,n& n1, M& M1, s& s1,f& f1,t& t1)
{
	m::iterator mit;
	typedef map<s, o> mDP;
	mDP mDP1;
	tstop ts;
	typedef pair <s, o> o_Pair; 
	typedef pair <s,tstop> tdPair;
	mapstop mStop;
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
	char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijklm, strijkl,strjklm,strjkl,strklm,strijk,str1;
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
// run dp routine

	for (j=n1;j>=n1-M1;j--)
	{
		strjklm =  to_string(j) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
		mmapIJDPK.insert(strlng_Pair(strjklm,n1)); // predecessor  
		mmapIJDPRes.insert(strdbl_Pair(strjklm,fcostar)); 
		mmapStopDPRes.insert(lngstr_Pair(n1,strjklm));
		mapDPRes.insert(lngdbl_Pair(n1,fcostar));
	}
	int q=-1;
	fileName(f1,"dpmissedic.txt",outfilename,q,-1);
	ofstream missfile( outfilename, ios::out | ios::trunc );

	// mDP.clear();
	for (m=n1;m>0;m--)
	{
		for (l=m;l>=m-M1 && l>0;l--)
		{
			for (k=l;k>=l-M1 && k>0;k--)
			{
				for (j=k;j>=k-M1 && j>0;j--)
				{
						strijkl = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l);
					for (i=j;i>=j-M1 && i>0;i--)
					{
						strjklm = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1+ to_string(m);
						strijklm = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1+ to_string(m);
						mit=m1.find(strijklm);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_tstop().get_TCost();
							if (icost<0) { icost=inf;}
							mmapKLMDPRit = mmapIJDPRes.find(strjklm);
							mapDPRit = mapDPRes.find(m);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (m==n1 && l==n1 && k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strijkl));
							mmapDPIC.insert(dblng_Pair(fcost,m));
							mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
						}
						else
						{	
							if (i==j||j==k||k==l||l==m) {
								continue;
							} else {
								cout<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
							}
						}
					} // i - loop
				} // j - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						m =mmapDPICit->second;
						strijkl =mmapICDPIJRit->second;
						strijk = strijkl.substr(strijkl.find_last_of(*t1)-1);
						mmapIJDPK.insert(strlng_Pair(strijkl,m)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(m,strijkl));
						mmapDPredK.insert(lng_Pair(k,m));

						mmapIJDPRes.insert(strdbl_Pair(strijkl,fcostar)); 
						mapDPRes.insert(lngdbl_Pair(k,fcostar));
					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
			} // k - loop
		} // l loop

	} // m loop

 if (missfile.is_open()) {
	 missfile.close();
 }
 missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
 fileName(f1,"dpobjtrc5d.txt",outfilename,q,-1);
 outfile.open( outfilename, ios::out );
 fileName(f1,"dptrc.txt",outfilename,q,-1);
 missfile.open(outfilename, ios::out  );
				// write the dp result data to a file
	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
//	outfile<<txthdr<<endl;
	o1.serializetexthdr(outfile);
//					writeTextObjectData(mDP1,o1,strijklm,outfile,txthdr);
//					writeTextObjectData(mStop,ts,strijklm,outfile,txthdr);

 // find the i1=1_j1=1 (i.e. "1_1") key value of k1 value then find the next i2=j1_j2=k1 value k2 recurse  
	i=1; 
	j=1;
	k=1;
	l=2;
	m=3;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
	while (mmapIJDPKit!=mmapIJDPK.end())
	{
		strijkl =mmapIJDPKit->first;
		strijk = strijkl.substr(0,strijkl.find_last_of(*t1));
		if (strijk== (to_string(i) + *t1 + to_string(j)+ *t1 + to_string(k))) {
			break;
		}
		mmapIJDPKit++;
	}

	mmapIJDPKit=mmapIJDPK.find(strijkl);
	if (mmapIJDPKit!=mmapIJDPK.end())
	{
			// find first *t1 
			strijk =strijkl; 
			i = getLoc(strijk,*t1,i);
//		i = from_String<long> strijkl.substr(1,i1-1);
		// get the substring beyound i
			strijk=getSubstr(strijk,*t1);
			// find next *t1 
			j = getLoc(strijk,*t1,j);
			// get the substring beyound j
			strijk=getSubstr(strijk,*t1);
			// find next *t1 
			k= getLoc(strijk,*t1,k);
			// get the substring beyound j
			strijk=getSubstr(strijk,*t1);
			// find last *t1 
			l= getLoc(strijk,*t1,l);		
		while (j<n1) 
		{
			mmapIJDPKit = mmapIJDPK.find(strijkl);
			if (mmapIJDPKit!=mmapIJDPK.end())
			{
				m=mmapIJDPKit->second;
				mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
			}
			else 
			{  if (m==n1)
				{
					strijkl =  to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;
					mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
				}

			}
			strijklm =  to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1+ to_string(m);
			mit = m1.find(strijklm);
			if (mit!=m1.end()) {
				mit->second.serializetext(outfile);
				mit->second.serializetext(missfile);
				o1 = mit->second;
				ts = o1.get_tstop();
				mDP1.insert (o_Pair(strijkl,mit->second));
//			mStop.insert (tdPair(strijkl,ts));
			}
			i=j;
			j=k;
			k=l;
			l=m;
			strijkl =  to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;

		}
	}
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

/* 
	fileName(f1,"dptrc.txt",outfilename,q,-1);
	outfile.open(outfilename, ios::out  );
				{// write the dp result data to a file
					writeTextData(mmapIJDPResICost,strijklm,j,outfile,"Stij\tK");
				}

if (outfile.is_open()) {
	 outfile.close();
}
 outfile.clear();
*/
 fileName(f1,"dpPatrc.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );
	o1.serializetexthdr(outfile);
				{// write the dp result data to a file
					writeTextObjectData(mDP1,o1,strijklm,outfile,"");
				}
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 // write dp trace data of the pathof optimization
 fileName(f1,"dptrest.txt",outfilename,q,-1);
 outfile.open(outfilename, ios::out  );

				{// write the dp result data to a file
					writeTextData(mmapIJDPK,strijklm,j,outfile,"Stage\tK");
					writeTextData(mmapIJDPRes,strijklm,fcost,outfile,"KeyPath\tOptimalCost");
					writeTextData(mmapStopDPRes,j,strijklm,outfile,"Stage\tKeyPath");
					writeTextData(mmapIJDPK,strijklm,j,outfile,"i-j\tk");
					writeTextData(mmapDPredK,j,i,outfile,"i\tj");
				}
			


 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
return m1;
 }


// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename l,typename m, typename o,typename n, typename  M,typename s,typename t,typename x>
l& dpbarun3d(l& l1,m& m1, o& o1,n& n1, M& M1, s& s1,t& t1, x& x1)
{
	m::iterator mit;
	typedef map<s, o> mDP;
	mDP mDP1;
	tstop ts;
	typedef pair <s, o> o_Pair; 
	typedef pair <s,tstop> tdPair;
	mapsstop mStop;
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
	//char outfilename[ MaxStrLen +1]=""; // data output file 
    //char  ext[30]="\0"; // file extension
    //char* fext=ext; // file extension
	string strijk, strij,strjk,str1;
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
// run dp routine
	mit = m1.begin();

	for (j=n1;j>=n1-M1;j--)
	{
		strjk =  to_string(j) + *t1 + to_string(n1) ;
		mmapIJDPK.insert(strlng_Pair(strjk,n1)); // predecessor  
		mmapIJDPRes.insert(strdbl_Pair(strjk,fcostar)); 
		mmapStopDPRes.insert(lngstr_Pair(n1,strjk));
		mapDPRes.insert(lngdbl_Pair(n1,fcostar));
	}
	int q=-1;

	string strext = x1;
	strext.append("DPMissedIC.txt");
	//strcpy(fext, strext.c_str());
	//fileName(f1,fext,outfilename,q,-1);
	ofstream missfile( strext.c_str(), ios::out | ios::trunc );

	for (j=n1;j>=1;j--)
	{
		for (i=j;i>=j-M1 && i>=1 ;i--)
		{
			strij = to_string(i) + *t1 + to_string(j) ;
			for (k=j;k<=j+M1 && k<=n1;k++)
			{
					strjk = to_string(j) + *t1 + to_string(k) ;
					strijk = to_string(i) + *t1 + strjk;
						mit=m1.find(strijk);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_tstop().get_TCost();
							if (icost<0) { 
								icost=inf;
							}
							mmapKLMDPRit = mmapIJDPRes.find(strjk);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strij));
							mmapDPIC.insert(dblng_Pair(fcost,k));
							mmapIJDPResICost.insert(strlng_Pair(strij,k));
						}
						else
						{	
							if (i==j||j==k) {
								continue;
							} else {
								cout<<"key " <<strijk<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijk<<" does not have an immediate cost"<<endl;
							}
						}
			} // k - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						k =mmapDPICit->second;
						strij =mmapICDPIJRit->second;
						// strij = strijk.substr(0,strijk.find_last_of(*t1));
						mmapIJDPK.insert(strlng_Pair(strij,k)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(k,strij));
						mmapDPredK.insert(lng_Pair(j,k));
						mmapIJDPRes.insert(strdbl_Pair(strij,fcostar));

					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
		} // i loop

	} // j loop

 if (missfile.is_open()) {
	 missfile.close();
 }
 missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
	strext = x1;
	strext.append("dpobjecTrace.txt");
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
    //fileName(f1,fext,outfilename,q,-1);
    outfile.open( strext.c_str(), ios::out );
	o1.serializetexthdr(outfile);
	strext = x1;
	strext.append("dptrace.txt");
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
    //fileName(f1,fext,outfilename,q,-1);
    missfile.open(strext.c_str(), ios::out  );
	o1.serializetexthdr(missfile);

 // find the from the i,j values from k+1 to k+M1  (i.e. "1_2 to 1_4") key values and 
 //	pick the minimum of these values as the starting point then find m and recurse  
	i=1; 
	j=1;
	k=1;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
    for (k=1;k<=M1+1;k++) 
	{
		strij =(to_string(i) + *t1 + to_string(j));
		mmapKLMDPRit = mmapIJDPRes.find(strij);
		if (mmapKLMDPRit != mmapIJDPRes.end()) {
			mmapICostDPIJ.insert(dblstr_Pair(mmapKLMDPRit->second,strij));
		}
	}
	mmapICDPIJRit = mmapICostDPIJ.begin();
	strij=mmapICDPIJRit->second;

	mmapIJDPKit=mmapIJDPK.find(strij);
		while (mmapIJDPKit!=mmapIJDPK.end()) 
		{
			k=mmapIJDPKit->second;
			mmapIJDPResICost.insert(strlng_Pair(strij,k));
			strijk =  strij + *t1+ to_string(k);
			mit = m1.find(strijk);
			if (mit!=m1.end()) {
				mit->second.serializetext(outfile);
				mit->second.serializetext(missfile);
				o1 = mit->second;
				ts = o1.get_tstop();
				mDP1.insert (o_Pair(strij,mit->second));
			}
			if (strij == (to_string(n1) + *t1 + to_string(n1)))
			{ 
				break;
			}
			strij = getSubstr(strij,*t1) + *t1 + to_string(k) ;
			mmapIJDPKit = mmapIJDPK.find(strij);
		}
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

	strext = x1;
	strext.append("DPatternTrace.txt");
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
    //fileName(f1,fext,outfilename,q,-1);
	outfile.open(strext.c_str(), ios::out  );
	o1.serializetexthdr(outfile);
 // write the dp result data to a file
 writeTextObjectData(mDP1,o1,strijk,outfile,"");
 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();

 // write the path trace to a file
	strext = x1;
	strext.append("DPathTrace.out");
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
    //fileName(f1,fext,outfilename,q,-1);
	outfile.open(strext.c_str(), ios::out  );

				{// write the dp result data to a file
					writeTextData(mmapIJDPK,strijk,j,outfile,"Stage\tK");
					writeTextData(mmapIJDPRes,strijk,fcost,outfile,"KeyPath\tOptimalCost");
					writeTextData(mmapStopDPRes,j,strijk,outfile,"Stage\tKeyPath");
					writeTextData(mmapIJDPK,strijk,j,outfile,"i-j\tk");
					writeTextData(mmapDPredK,j,i,outfile,"i\tj");
				}
			


 if (outfile.is_open()) {
	 outfile.close();
 }
 outfile.clear();
  return mDP1;
 }



// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename m, typename o,typename n, typename  M,typename s,typename f,typename t>
m& dpbarun3dx(m& m1, o& o1,n& n1, M& M1, s& s1,f& f1,t& t1)
{
	m::iterator mit;
	typedef map<s, o> mDP;
	mDP mDP1;
	tstop ts;
	typedef pair <s, o> o_Pair; 
	typedef pair <s,tstop> tdPair;
	mapstop mStop;
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
	char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijk, strij,strjk,str1;
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
// run dp routine

	for (j=n1;j>=n1-M1;j--)
	{
		strjk =  to_string(j) + *t1 + to_string(n1) ;
		mmapIJDPK.insert(strlng_Pair(strjk,n1)); // predecessor  
		mmapIJDPRes.insert(strdbl_Pair(strjk,fcostar)); 
		mmapStopDPRes.insert(lngstr_Pair(n1,strjk));
		mapDPRes.insert(lngdbl_Pair(n1,fcostar));
	}
	int q=-1;
	fileName(f1,"dpmiss3dic.txt",outfilename,q,-1);
	ofstream missfile( outfilename, ios::out | ios::trunc );

	for (j=n1;j>=1;j--)
	{
		for (i=j;i>=i-M1 && j>=1 ;j--)
		{
			strij = to_string(i) + *t1 + to_string(j) ;
			for (k=j;k<=j+M1 && k<=n1;k++)
			{
					strjk = to_string(j) + *t1 + to_string(k) ;
					strijk = to_string(i) + *t1 + strjk;
						mit=m1.find(strijk);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_tstop().get_TCost();
							if (icost<0) { 
								icost=inf;
							}
							mmapKLMDPRit = mmapIJDPRes.find(strjk);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strij));
							mmapDPIC.insert(dblng_Pair(fcost,k));
							mmapIJDPResICost.insert(strlng_Pair(strij,k));
						}
						else
						{	
							if (i==j||j==k) {
								continue;
							} else {
								cout<<"key " <<strijk<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijk<<" does not have an immediate cost"<<endl;
							}
						}
			} // k - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						k =mmapDPICit->second;
						strij =mmapICDPIJRit->second;
						// strij = strijk.substr(0,strijk.find_last_of(*t1));
						mmapIJDPK.insert(strlng_Pair(strij,k)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(k,strij));
						mmapDPredK.insert(lng_Pair(j,k));
						mmapIJDPRes.insert(strdbl_Pair(strij,fcostar));

					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
		} // i loop

	} // j loop

	if (missfile.is_open()) {
	 missfile.close();
	}
	missfile.clear();


	 // trace the optimal path from the above result
	// write the resultant dp path
	 fileName(f1,"dpobjtrc3d.txt",outfilename,q,-1);
	 outfile.open( outfilename, ios::out );
	 o1.serializetexthdr(outfile);
	 fileName(f1,"dptrc3d.txt",outfilename,q,-1);
	 missfile.open(outfilename, ios::out  );
	 o1.serializetexthdr(missfile);

	 // find the from the i,j values from k+1 to k+M1  (i.e. "1_2 to 1_4") key values and 
	 //	pick the minimum of these values as the starting point then find m and recurse  
	i=1; 
	j=1;
	k=1;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
	for (k=1;k<=M1+1;k++) 
	{
		strij =(to_string(i) + *t1 + to_string(j));
		mmapKLMDPRit = mmapIJDPRes.find(strij);
		if (mmapKLMDPRit != mmapIJDPRes.end()) {
			mmapICostDPIJ.insert(dblstr_Pair(mmapKLMDPRit->second,strij));
		}
	}
	mmapICDPIJRit = mmapICostDPIJ.begin();
	strij=mmapICDPIJRit->second;

	mmapIJDPKit=mmapIJDPK.find(strij);
	while (mmapIJDPKit!=mmapIJDPK.end()) 
	{
		k=mmapIJDPKit->second;
		mmapIJDPResICost.insert(strlng_Pair(strij,k));
		strijk =  strij + *t1+ to_string(k);
		mit = m1.find(strijk);
		if (mit!=m1.end()) {
			mit->second.serializetext(outfile);
			mit->second.serializetext(missfile);
			o1 = mit->second;
			ts = o1.get_tstop();
			mDP1.insert (o_Pair(strij,mit->second));
		}
		if (strij == (to_string(n1) + *t1 + to_string(n1)))
		{ 
			break;
		}
		strij = getSubstr(strij,*t1) + *t1 + to_string(k) ;
		mmapIJDPKit = mmapIJDPK.find(strij);
	}
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

	fileName(f1,"dpPatrc3d.txt",outfilename,q,-1);
	outfile.open(outfilename, ios::out  );
	o1.serializetexthdr(outfile);
	// write the dp result data to a file
	writeTextObjectData(mDP1,o1,strijk,outfile,"");
	if (outfile.is_open()) {
	 outfile.close();
	}
	outfile.clear();

 // add the cost at the start of the route
	fileName(f1,"dptrest3d.txt",outfilename,q,-1);
	outfile.open(outfilename, ios::out  );

	{// write the dp result data to a file
		writeTextData(mmapIJDPK,strijk,j,outfile,"Stage\tK");
		writeTextData(mmapIJDPRes,strijk,fcost,outfile,"KeyPath\tOptimalCost");
		writeTextData(mmapStopDPRes,j,strijk,outfile,"Stage\tKeyPath");
		writeTextData(mmapIJDPK,strijk,j,outfile,"i-j\tk");
		writeTextData(mmapDPredK,j,i,outfile,"i\tj");
	}
			


	if (outfile.is_open()) {
	 outfile.close();
	}
	outfile.clear();
	return m1;
}



// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename l,typename m, typename o,typename n, typename  p,typename s,typename t,typename x,typename y,typename g,typename h,typename d>
l& dpbarun5d(l& l1,m& m1, o& o1,n& n1, p& M1, s& s1,t& t1, x& x1,y& y1,g& gc,h& ph,d& dbo)
{
	m::iterator mit;
	typedef map<s, o> mDP;
	mDP mDP1;
	tstop ts;
	typedef pair <s, o> o_Pair; 
	typedef pair <s,tstop> tdPair;
	mapsstop mStop;
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
//	char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijklm, strijkl,strjklm,strjkl,strklm,strijk,fname,strext;
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0,wcost=inf,rcost=inf,ocost=inf;
// run dp routine

	for (j=n1;j>=n1-M1;j--)
	{
		strjklm =  to_string(j) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
		mmapIJDPK.insert(strlng_Pair(strjklm,n1)); // predecessor  
		mmapIJDPRes.insert(strdbl_Pair(strjklm,fcostar)); 
		mmapStopDPRes.insert(lngstr_Pair(n1,strjklm));
		mapDPRes.insert(lngdbl_Pair(n1,fcostar));
	}
	int q=-1;
	strext = x1;
	strext.append("DPMissedIC.txt");
	//const char * outfilename5 = strext.c_str();
	ofstream missfile( strext.c_str(), ios::out | ios::trunc );
	// mDP.clear();

	for (k=n1;k>=1;k--)
	{
		for (j=k;j>=k-M1 && j>=1 ;j--)
		{
			for (i=j;i>=j-M1 && i>=1;i--)
			{
				for (l=k;l<=k+M1 && l<=n1;l++)
				{
					strijkl = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;
					for (m=l;m<=l+M1 && m<=n1;m++)
					{
						strjklm = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1 + to_string(m);
						strijklm = to_string(i) + *t1 + strjklm;
						mit=m1.find(strijklm);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_tstop().get_TCost();
							wcost = mit->second.get_tstop().get_WalkCost();
							rcost = mit->second.get_tstop().get_RideCost();
							ocost = mit->second.get_tstop().get_OperCost();
							if (icost<0 ||wcost<0 ||rcost<0 ||ocost<0  ) { 
								icost=inf;
							}
							//if (wcost<0) { 
							//	wcost=inf;
							//}
							//if (rcost<0 ) { 
							//	rcost=inf;
							//}
							//if (ocost<0 ) { 
							//	ocost=inf;
							//}
							mmapKLMDPRit = mmapIJDPRes.find(strjklm);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (m==n1 && l==n1 && k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							if (fcost<0) {fcost = inf;}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strijkl));
							mmapDPIC.insert(dblng_Pair(fcost,m));
							mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
						}
						else
						{	
							if (i==j||j==k||k==l||l==m) {
								continue;
							} else {
								cout<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
							}
						}
					} // i - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() && mmapDPICit != mmapDPIC.end()  )
					{
						fcostar =mmapICDPIJRit->first;
						m = mmapDPICit->second;
						strijkl =mmapICDPIJRit->second;
						strijk = strijkl.substr(0,strijkl.find_last_of(*t1));
						mmapIJDPK.insert(strlng_Pair(strijkl,m)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(m,strijkl));
						mmapDPredK.insert(lng_Pair(k,m));
						mmapIJDPRes.insert(strdbl_Pair(strijkl,fcostar));

					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
				} // j - loop
			} // k - loop
		} // l loop

	} // m loop

	 if (missfile.is_open()) {
		 missfile.close();
	 }
	 missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
	strext = x1;
	strext.append("dpobjecTrace.txt");
//	fext[0]='\0';
//	strcpy(fext, strext.c_str());
//    fileName(f1,fext,outfilename,q,-1);
	//const char * outfilename4 = strext.c_str();
	outfile.open( strext.c_str(), ios::out );
	//outfilename4=0;
	o1.serializetexthdr(outfile);
	strext = x1;
	strext.append("dptrace.txt");
//    fileName(f1,fext,outfilename,q,-1);
	//const char * outfilename3 = strext.c_str();
	missfile.open(strext.c_str(), ios::out  );
	o1.serializetexthdr(missfile);
	//outfilename3=0;
	// write the dp result data to a file
//	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
//	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
//	outfile<<txthdr<<endl;
//					writeTextObjectData(mDP1,o1,strijklm,outfile,txthdr);
//					writeTextObjectData(mStop,ts,strijklm,outfile,txthdr);

 // find the from the i,j,k & l values from k+1 to k+M1  (i.e. if M1 = 2 "1_1_1_2 to 1_1_1_4") key values and 
 //	pick the minimum of these values as the starting point then find m and recurse  
	i=1; 
	j=1;
	k=1;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
    for (l=1;l<=M1+1;l++) 
	{
		strijkl =(to_string(i) + *t1 + to_string(j)+ *t1 + to_string(k)+ *t1 + to_string(l));
		mmapKLMDPRit = mmapIJDPRes.find(strijkl);
		if (mmapKLMDPRit != mmapIJDPRes.end()) {
			mmapICostDPIJ.insert(dblstr_Pair(mmapKLMDPRit->second,strijkl));
		}
	}
	mmapICDPIJRit = mmapICostDPIJ.begin();
	strijkl=mmapICDPIJRit->second;

	mmapIJDPKit=mmapIJDPK.find(strijkl);
		while (mmapIJDPKit!=mmapIJDPK.end()) 
		{
			m=mmapIJDPKit->second;
			mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
			strijklm =  strijkl + *t1+ to_string(m);
			mit = m1.find(strijklm);
			if (mit!=m1.end()) {
				mit->second.serializetext(outfile);
				mit->second.serializetext(missfile);
				o1 = mit->second;
				ts = o1.get_tstop();
				l1.insert (o_Pair(strijkl,mit->second));
//			mStop.insert (tdPair(strijkl,ts));
			}
			if (strijkl == (to_string(m) + *t1 + to_string(m)+ *t1 + to_string(m)+ *t1 + to_string(m)))
			{ 
				break;
			}
			strijkl = getSubstr(strijkl,*t1) + *t1 + to_string(m) ;
			mmapIJDPKit = mmapIJDPK.find(strijkl);
		}
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

	strext = x1;
	strext.append("DPatternTrace.txt");
	//const char * outfname = new char[strext.length() + 1];
	//const char * outfilename2 = strext.c_str();
	outfile.open(strext.c_str(), ios::out  );
	o1.serializetexthdr(outfile);

 // write the dp result data to a file
	writeTextObjectData(l1,o1,strijklm,outfile,"");
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	//outfilename2=0;

// create result table of Stops 
	short D = gc->get_dpdimension();

	string tblDP = y1.get_tblstop() + "D" + to_string<short>(D);
	tblDP.append("t" + to_string<long> (ph.tripId()) + ph.tripPeriod());
	tblDP.append("W" + to_string<double> ((int) (gc->get_walkcost()*10)));
	tblDP.append("R" + to_string<double> ((int) (gc->get_ridecost()*10)));
	replace(tblDP.begin(),tblDP.end(),' ','_');
	ReplaceAll2(tblDP,"__","_");

	int intSRID = y1.get_srid();

	string str1 =	"( DPKey Text,StOrdr Long, StopId Long,StopName Text, Walk Integer, Ride Integer, " 
		"ScId integer, I integer, J integer, K integer, L integer, M integer,CumDist Double,"
		" CRdTm Double,undCRdTm Double,CRdTmC Double, Ons Double,Offs Double,DepVol Double,"
		"probStop Double,depDelay Double, arrDelay Double,"
		"stopDelay Double,dwlDelay Double,rideDelay Double,PVal Double,"
		"AVal Double, WkTmOns Double,WkTmOffs Double,WalkCost Double,RideCost Double,"
		"OperCost Double, TCost Double , dpCnt long  ); ";

		//strtblName = x1.get_tbltrip() + stRW + "_MPDP" ; 
	string strtblName = tblDP + "_DPTrace" ; 
	
	replace(strtblName.begin(),strtblName.end(),' ','_');
	ReplaceAll2(strtblName,"__","_");

	bool blnCreate = createSpaTbl(str1,strtblName,dbo,intSRID,outfile);

	if (!blnCreate) {
	// exit since the DP Result table could not be created 
		outfile << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
		cout << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl;
		cout << "Press enter to exit!"<<endl;
		cin>>strtblName;
		exit (0);
		
	}	

	
	string strDPSQLInserTblDef =  "(\"DPKey\",\"StopId\",\"StopName\",\"Walk\",\"Ride\",\"ScId\",\"I\",\"J\",\"K\",\"L\",\"M\","
			"\"CRdTm\",\"undCRdTm\",\"Ons\",\"Offs\",\"DepVol\",\"probStop\",\"depDelay\", "  
			"\"arrDelay\",\"dwlDelay\",\"rideDelay\",\"PVal\",\"AVal\",\"WkTmOns\",\"WkTmOffs\",\"WalkCost\","
			"\"RideCost\",\"OperCost\",\"TCost\",\"dpCnt\",\"stopDelay\" ,\"StOrdr\",\"CRdTmC\",\"CumDist\") "  
			" Values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,"
					 "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);" ;
	
	string 	strDPSQLInserTbl = "Insert into " + strtblName + "  " + strDPSQLInserTblDef ;
	
	// make the trip dpmap
	long ip = 0;//kStop1.tripId();

//	tripDPsmap = dpTripWalkRideMap(sdptsidpmmap,tripDPsmap,ip,dpstop,tDPstp,gcost);
	bool blnIns = inSpaTblDPStop(strDPSQLInserTbl,dbo,l1,o1,ip,gc,outfile);
	if (blnIns) {

	}


 // write the path trace to a file
	strext = x1;
	strext.append("DPathTrace.out");
	//const char * outfilename1 = strext.c_str();
	outfile.open(strext.c_str(), ios::out  );

	{// write the dp result data to a file
		writeTextData(mmapIJDPK,strijklm,j,outfile,"Stage\tK");
		writeTextData(mmapIJDPRes,strijklm,fcost,outfile,"KeyPath\tOptimalCost");
		writeTextData(mmapStopDPRes,j,strijklm,outfile,"Stage\tKeyPath");
		writeTextData(mmapIJDPK,strijklm,j,outfile,"i-j\tk");
		writeTextData(mmapDPredK,j,i,outfile,"i\tj");
	}
	 if (outfile.is_open()) {
		 outfile.close();
	 }
	outfile.clear();
	//outfilename1=0;
	return l1;
 }

// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename l,typename m, typename o,typename n, typename  M,typename s,typename f,typename t,typename x,typename g,typename p,typename d,typename r,typename u,typename v>
l& dpOptimalMultiPd5dbz(l& mDP1,m& m1, o& o1,n& n1, M& M1, s& tblDP,f& f1,t& t1, x& x1,g& g1,p& p1,d& outdb,r& dpstop,u& stop0,v& v1)
{
	m::iterator mit;
	typedef map<s, o> mDP;
	//mDP mDP1;
	typedef pair <s, o> o_Pair; 
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
	int M = 0, D=0, q = 0, intSRID=0,ip=0;
	bool blnCreate=false, blnIns;
	//char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijklm="", strijkl="",strjklm="",strjkl="",strklm="",strijk="",str1="",outFNm="",stRW="", strtblName="";
    //char  ext[50]="\0"; // file extension
    //char* fext=ext; // file extension
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
	M = g1->get_maxskip() + 1;

	D = g1->get_dpdimension();
	mmapIJDPK.clear();
	mmapIJDPRes.clear();
	mmapStopDPRes.clear();
	mapDPRes.clear();
	mmapICostDPIJ.clear();
	mmapDPIC.clear();
	mmapIJDPResICost.clear(); 
	mmapIJDPK.clear(); // predecessor  
	mmapStopDPRes.clear();
	mmapDPredK.clear();
	
// run dp routine
	f outFN = "";
	for (j=n1;j>=n1-M1;j--)
	{
		strjklm =  to_string(j) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
		mmapIJDPK.insert(strlng_Pair(strjklm,n1)); // predecessor  
		mmapIJDPRes.insert(strdbl_Pair(strjklm,fcostar)); 
		mmapStopDPRes.insert(lngstr_Pair(n1,strjklm));
		mapDPRes.insert(lngdbl_Pair(n1,fcostar));
	}
	stRW = "D" + to_string<short>(D);
	stRW.append("W" + to_string<double> ((int) (g1->get_walkcost()*10)));
	stRW.append("R" + to_string<double> ((int) (g1->get_ridecost()*10)));

	outFN = f1 + stRW + "MPDPMissedIC5d.txt";
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
	//fileName(f1,fext,outfilename,q,-1);
	ofstream missfile( outFN.c_str(), ios::out | ios::trunc );

	// mDP.clear();
	for (k=n1;k>=1;k--)
	{
		for (j=k;j>=k-M1 && j>=1 ;j--)
		{
			for (i=j;i>=j-M1 && i>=1;i--)
			{
				for (l=k;l<=k+M1 && l<=n1;l++)
				{
					strijkl = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;
					for (m=l;m<=l+M1 && m<=n1;m++)
					{
						strjklm = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1 + to_string(m);
						strijklm = to_string(i) + *t1 + strjklm;
						mit=m1.find(strijklm);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_TCost();
							if (icost<0) { 
								icost=inf;
							}
							mmapKLMDPRit = mmapIJDPRes.find(strjklm);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (m==n1 && l==n1 && k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strijkl));
							mmapDPIC.insert(dblng_Pair(fcost,m));
							mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
						}
						else
						{	
							if (i==j||j==k||k==l||l==m) {
								continue;
							} else {
								cout<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
							}
						}
					} // m - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						if ( mmapDPICit !=mmapDPIC.end()) {
							m =mmapDPICit->second;
						}
						strijkl =mmapICDPIJRit->second;
						strijk = strijkl.substr(0,strijkl.find_last_of(*t1));
						mmapIJDPK.insert(strlng_Pair(strijkl,m)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(m,strijkl));
						mmapDPredK.insert(lng_Pair(k,m));
						mmapIJDPRes.insert(strdbl_Pair(strijkl,fcostar));

					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
				} // l - loop
			} // i - loop
		} // j loop

	} // k loop

	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
	outFN = f1 + stRW + "MPDPResult.txt";
	outfile.open( outFN.c_str(), ios::out );
	o1.serializetext(outfile);
	outFN = f1 + stRW + "MPDPTrace.txt";
    //fileName(f1,fext,outfilename,q,-1);
	missfile.open(outFN.c_str(), ios::out  );
	o1.dpAccumHeader(missfile);
				// write the dp result data to a file
//	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
//	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
//	outfile<<txthdr<<endl;
//					writeTextObjectData(mDP1,o1,strijklm,outfile,txthdr);
//					writeTextObjectData(mStop,ts,strijklm,outfile,txthdr);

 // find the from the i,j,k & l values from k+1 to k+M1  (i.e. "1_1_1_2 to 1_1_1_4") key values and 
 //	pick the minimum of these values as the starting point then find m and recurse  
	i=1; 
	j=1;
	k=1;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
    for (l=1;l<=M1+1;l++) 
	{
		strijkl =(to_string(i) + *t1 + to_string(j)+ *t1 + to_string(k)+ *t1 + to_string(l));
		mmapKLMDPRit = mmapIJDPRes.find(strijkl);
		if (mmapKLMDPRit != mmapIJDPRes.end()) {
			mmapICostDPIJ.insert(dblstr_Pair(mmapKLMDPRit->second,strijkl));
		}
	}
	mmapICDPIJRit = mmapICostDPIJ.begin();
	strijkl=mmapICDPIJRit->second;

	mmapIJDPKit=mmapIJDPK.find(strijkl);
	while (mmapIJDPKit!=mmapIJDPK.end()) 
	{
		m=mmapIJDPKit->second;
		mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
		strijklm =  strijkl + *t1+ to_string(m);
		mit = m1.find(strijklm);
		if (mit!=m1.end()) {
			outfile<<mit->second;
			missfile<<mit->second;
			o1 = mit->second;
			mDP1.insert (o_Pair(strijkl,o1));
//			mStop.insert (tdPair(strijkl,ts));
		}
		if (strijkl == (to_string(n1) + *t1 + to_string(n1)+ *t1 + to_string(n1)+ *t1 + to_string(n1)))
		{ 
			break;
		}
		strijkl = getSubstr(strijkl,*t1) + *t1 + to_string(m) ;
		mmapIJDPKit = mmapIJDPK.find(strijkl);
	}
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

	outFN = f1 + stRW + "MPDPaTrace.txt";
	outfile.open(outFN.c_str(), ios::out  );
	outfile<<"Summary of DP Results "<<endl<<"DPKey\t";
	o1.dpAccumHeader(outfile);
// write the dp result data to a file
	writeTextData(mDP1,strijklm,o1,outfile,"");
	if (outfile.is_open()) {
		 outfile.close();
	}
	outfile.clear();

// write the dp result to a table
// create result table of Stops 
	intSRID = x1.get_srid();

	str1 =	"( DPKey Text,StOrdr Long, StopId Long,StopName Text, Walk Integer, Ride Integer, " 
		"ScId integer, I integer, J integer, K integer, L integer, M integer,"
		" CRdTm Double,undCRdTm Double,CRdTmC Double, Ons Double,Offs Double,DepVol Double,"
		"probStop Double,depDelay Double, arrDelay Double,"
		"stopDelay Double,dwlDelay Double,rideDelay Double,PVal Double,"
		"AVal Double, WkTmOns Double,WkTmOffs Double,WalkCost Double,RideCost Double,"
		"OperCost Double, TCost Double , dpCnt long  ); ";

		//strtblName = x1.get_tbltrip() + stRW + "_MPDP" ; 
		strtblName = tblDP + "_MPDPTrace" ; 
	
	replace(strtblName.begin(),strtblName.end(),' ','_');
	ReplaceAll2(strtblName,"__","_");

	blnCreate = createSpaTbl(str1,strtblName,outdb,intSRID,outfile);

	if (!blnCreate) {
	// exit since the DP Result table could not be created 
		outfile << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
		cout << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl;
		cout << "Press enter to exit!"<<endl;
		cin>>strtblName;
		exit (0);
		
	}	

	
	string strDPSQLInserTblDef =  "(\"DPKey\",\"StopId\",\"StopName\",\"Walk\",\"Ride\",\"ScId\",\"I\",\"J\",\"K\",\"L\",\"M\","
			"\"CRdTm\",\"undCRdTm\",\"Ons\",\"Offs\",\"DepVol\",\"probStop\",\"depDelay\", "  
			"\"arrDelay\",\"dwlDelay\",\"rideDelay\",\"PVal\",\"AVal\",\"WkTmOns\",\"WkTmOffs\",\"WalkCost\","
			"\"RideCost\",\"OperCost\",\"TCost\",\"dpCnt\",\"stopDelay\" ,\"StOrdr\",\"CRdTmC\") "  
			" Values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,"
					 "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);" ;
	
	string 	strDPSQLInserTbl = "Insert into " + strtblName + "  " + strDPSQLInserTblDef ;
	
	// call insert table routine to populate the stop table
	//	blnInsert = inSpaTblStop(strSQLInserTbl,pTransOutDb,*pstop,stops,tblTrip, keyFld , intSRID,logfile);
	// make the trip dpmap
	ip = 0;//kStop1.tripId();

//	tripDPsmap = dpTripWalkRideMap(sdptsidpmmap,tripDPsmap,ip,dpstop,tDPstp,gcost);
	blnIns = inSpaTblDPStop(strDPSQLInserTbl,outdb,mDP1,o1,dpstop,g1,v1);
	if (blnIns) {
		//strSQL = " select t2.RteName,t2.SchlName,t2.DirName, t2.StOrdr,t2.Stop_Id,"
		//	"t2.EdgeId,T2.PosAlong,t1.ScID, t1.DPKey, t1.i, t1.j, t1.k, t1.l, t1.m, t1.Cnt,"
		//	"t1.StopName, t1.PVal, t1.Aval, t1.Ons, t1.Offs, t1.DepVol,t1.probStop, t1.depDelay,"
		//	"t1.arrDelay, t1.dwlDelay, t1.RideDelay,t1.WkTmOns, t1.WkTmOffs, t1.WalkCost,t1.RideCost,"
		//	"t1.OperCost,t1.TCost," + " MakePoint(XC,YC," + (to_string<int>(intSRID)) + " ) Geometry "
		//	" FROM " + strtblName + " t1 , " + RTD_Stops_Rte12SBE + " t2 "
		//	" where t1.StopName = t2.StrMainCrs "
		//	" ORDER BY K ;"

		//blnCreate = execSpatialQuery(strSQL,strDPTblName,pTransOutDb,intSRID,outfile);

		//if (blnCreate) {

		//	strSQL = " SELECT recovergeometrycolumn('" + strDPTblName + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Point',2);";
		//	blnCreate = recoverSpatialGeometry(strSQL,strDPTblName,pTransOutDb,intSRID,outfile);
		//}

	}
 // write the path trace to a file
	outFN = f1 + stRW + "MPDPathTrace.out";
	outfile.open(outFN.c_str(), ios::out  );

	{// write the dp result data to a file
		writeTextData(mmapIJDPK,strijklm,j,outfile,"Stage\tK");
		writeTextData(mmapIJDPRes,strijklm,fcost,outfile,"KeyPath\tOptimalCost");
		writeTextData(mmapStopDPRes,j,strijklm,outfile,"Stage\tKeyPath");
		writeTextData(mmapIJDPK,strijklm,j,outfile,"i-j\tk");
		writeTextData(mmapDPredK,j,i,outfile,"i\tj");
	}
			
	if (outfile.is_open()) {
		 outfile.close();
	}
	outfile.clear();
	return mDP1;
 }


 // backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
// elist - DP fixed stops (not to be eliminated) in the DP algorithm  
template <typename l, typename m, typename o,typename n, typename  M,typename s,typename f,typename t,typename x,typename g,typename p,typename d,typename r,typename u,typename e,typename v>
l& dpOptimalMultiPd5dbxe(l& mDP1, m& m1, o& o1,n& n1, M& M1, s& tblDP,f& f1,t& t1, x& x1,g& g1,p& p1,d& outdb,r& dpstop,u& stop0,e& elist, v& v1)
{
	m::iterator mit;
	typedef pair <s, o> o_Pair; 
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0,i0=0,j0=0,k0=0,l0=0,m0=0;
	int M = 0, D=0, q = 0, intSRID=0,ip=0;
	bool blnCreate=false, blnIns=false, blnFix=false;
	//char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijklm="", strijkl="",strjklm="",strjkl="",strklm="",strijk="",str1="",outFNm="",stRW="", strtblName="";
    //char  ext[50]="\0"; // file extension
    //char* fext=ext; // file extension
// at stop n , f*(n) = 0
	double icost=0,opcost=0,wkcost=0,rdcost=0, fcost=0,fcostar=0;
	M = g1->get_maxskip() + 1;

	D = g1->get_dpdimension();
	mmapIJDPK.clear();
	mmapIJDPRes.clear();
	mmapStopDPRes.clear();
	mapDPRes.clear();
	mmapICostDPIJ.clear();
	mmapDPIC.clear();
	mmapIJDPResICost.clear(); 
	mmapIJDPK.clear(); // predecessor  
	mmapStopDPRes.clear();
	mmapDPredK.clear();


// run dp routine
	f outFN = "";
	for (j=n1;j>=n1-M1;j--)
	{
		blnFix = false;
		for (j0=j+1;j0<n1;j0++) {
			if ((elist.find(j0) != elist.end() ) ) {// a fixed stop is skipped in this jklm pair 
				blnFix = true;
			} 
		}
			strjklm =  to_string(j) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);

		if (blnFix) { // skip this key from being in solution
			// strjklm =  to_string(j) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
		} else {
			mmapIJDPK.insert(strlng_Pair(strjklm,n1)); // predecessor  
			mmapIJDPRes.insert(strdbl_Pair(strjklm,fcostar)); 
			mmapStopDPRes.insert(lngstr_Pair(n1,strjklm));
			mapDPRes.insert(lngdbl_Pair(n1,fcostar));
		}
	}
	stRW = "D" + to_string<short>(D);
	stRW.append("W" + to_string<double> ((int) (g1->get_walkcost()*10)));
	stRW.append("R" + to_string<double> ((int) (g1->get_ridecost()*10)));

	outFN = f1 + tblDP + "_MPDPMissIC5d.txt";
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
	//fileName(f1,fext,outfilename,q,-1);
	ofstream missfile( outFN.c_str(), ios::out | ios::trunc );

	mDP1.clear();
	for (k=n1;k>=1;k--)
	{
		for (j=k;j>=k-M1 && j>=1 ;j--)
		{
			blnFix = listmember(j,k,elist);
			for (i=j;i>=j-M1 && i>=1;i--)
			{
				for (l=k;l<=k+M1 && l<=n1;l++)
				{
					strijkl = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;
					for (m=l;m<=l+M1 && m<=n1;m++)
					{
						strjklm = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1 + to_string(m);
						strijklm = to_string(i) + *t1 + strjklm;
						blnFix = listmember(i,j,elist);
						blnFix = listmember(j,k,elist);
						blnFix = listmember(k,l,elist);
						blnFix = listmember(l,m,elist);
						if (blnFix) { 
							icost=inf;
						}
						mit=m1.find(strijklm);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = o1.get_tstop().get_TCost();
							opcost = o1.get_tstop().get_OperCost();
							rdcost = o1.get_tstop().get_RideCost();
							wkcost = o1.get_tstop().get_WalkCost();
							if (icost<0 || opcost<0 || rdcost<0 || wkcost<0 || blnFix) { 
								icost=inf;
							}
							mmapKLMDPRit = mmapIJDPRes.find(strjklm);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (m==n1 && l==n1 && k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strijkl));
							mmapDPIC.insert(dblng_Pair(fcost,m));
							mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
						}
						else
						{	
							if (i==j||j==k||k==l||l==m) {
								continue;
							} else {
								cout<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								v1<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
							}
						}
					} // i - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						if ( mmapDPICit !=mmapDPIC.end()) {
							m =mmapDPICit->second;
						}
						strijkl =mmapICDPIJRit->second;
						//strijk = strijkl.substr(0,strijkl.find_last_of(*t1));
						mmapIJDPK.insert(strlng_Pair(strijkl,m)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(m,strijkl));
						mmapDPredK.insert(lng_Pair(k,m));
						mmapIJDPRes.insert(strdbl_Pair(strijkl,fcostar));

					}
					mmapICostDPIJ.clear();
					mmapIJDPResICost.clear();
					mmapDPIC.clear();
				} // j - loop
			} // k - loop
		} // l loop

	} // m loop

	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
	outFN = f1 + tblDP + "_MPDPResult.txt";
	outfile.open( outFN.c_str(), ios::out );
	o1.serializetexthdr(outfile);
	outFN = f1 + tblDP + "_MPDPTrace.txt";
    //fileName(f1,fext,outfilename,q,-1);
	missfile.open(outFN.c_str(), ios::out  );
	missfile<<o1;
	v1<<o1;
				// write the dp result data to a file
//	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
//	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
//	outfile<<txthdr<<endl;
//					writeTextObjectData(mDP1,o1,strijklm,outfile,txthdr);
//					writeTextObjectData(mStop,ts,strijklm,outfile,txthdr);

 // find the from the i,j,k & l values from k+1 to k+M1  (i.e. "1_1_1_2 to 1_1_1_4") key values and 
 //	pick the minimum of these values as the starting point then find m and recurse  
	i=1; 
	j=1;
	k=1;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
    for (l=1;l<=M1+1;l++) 
	{
		strijkl =(to_string(i) + *t1 + to_string(j)+ *t1 + to_string(k)+ *t1 + to_string(l));
		mmapKLMDPRit = mmapIJDPRes.find(strijkl);
		if (mmapKLMDPRit != mmapIJDPRes.end()) {
			mmapICostDPIJ.insert(dblstr_Pair(mmapKLMDPRit->second,strijkl));
		}
	}
	mmapICDPIJRit = mmapICostDPIJ.begin();
	if (mmapICDPIJRit != mmapICostDPIJ.end())
	{
		strijkl=mmapICDPIJRit->second;
		mmapIJDPKit=mmapIJDPK.find(strijkl);
		while (mmapIJDPKit!=mmapIJDPK.end()) 
		{
			m=mmapIJDPKit->second;
			mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
			strijklm =  strijkl + *t1+ to_string(m);
			mit = m1.find(strijklm);
			if (mit!=m1.end()) {
				outfile<<mit->second;
				missfile<<mit->second;
				o1 = mit->second;
				mDP1.insert (o_Pair(strijkl,o1));
	//			mStop.insert (tdPair(strijkl,ts));
			}
			if (strijkl == (to_string(n1) + *t1 + to_string(n1)+ *t1 + to_string(n1)+ *t1 + to_string(n1)))
			{ 
				break;
			}
			strijkl = getSubstr(strijkl,*t1) + *t1 + to_string(m) ;
			mmapIJDPKit = mmapIJDPK.find(strijkl);
		}
	} else {
		cout<<"key " <<strijkl<<" does not have an immediate cost"<<endl;
		missfile<<"key " <<strijkl<<" does not have an immediate cost"<<endl;
		v1<<"key " <<strijkl<<" does not have an immediate cost"<<endl;
	}

	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

	outFN = f1 + tblDP + "MPDPaTrace.txt";
	outfile.open(outFN.c_str(), ios::out  );
	outfile<<"Summary of DP Results "<<endl<<"DPKey\t";
	//o1.serializetexthdr(outfile);
	o1.serializetexthdr(outfile);// write the dp result data to a file
	writeTextData(mDP1,strijklm,o1,outfile,"");
	if (outfile.is_open()) {
		 outfile.close();
	}
	outfile.clear();

// write the dp result to a table
// create result table of Stops 
	intSRID = x1.get_srid();

	str1 =	"( DPKey Text,StOrdr Long, StopId Long,StopName Text, Walk Integer, Ride Integer, " 
		"ScId integer, I integer, J integer, K integer, L integer, M integer,"
		" CRdTm Double,undCRdTm Double,CRdTmC Double, Ons Double,Offs Double,DepVol Double,"
		"probStop Double,depDelay Double, arrDelay Double,"
		"stopDelay Double,dwlDelay Double,rideDelay Double,PVal Double,"
		"AVal Double, WkTmOns Double,WkTmOffs Double,WalkCost Double,RideCost Double,"
		"OperCost Double, TCost Double , dpCnt long  ); ";

		//strtblName = x1.get_tbltrip() + stRW + "_MPDP" ; 
		strtblName = tblDP + "_DPTrace" ; 
	
	replace(strtblName.begin(),strtblName.end(),' ','_');
	ReplaceAll2(strtblName,"__","_");

	blnCreate = createSpaTbl(str1,strtblName,outdb,intSRID,outfile);

	if (!blnCreate) {
	// exit since the DP Result table could not be created 
		outfile << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
		cout << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl;
		cout << "Press enter to exit!"<<endl;
		cin>>strtblName;
		exit (0);
		
	}	

	
	string strDPSQLInserTblDef =  "(\"DPKey\",\"StopId\",\"StopName\",\"Walk\",\"Ride\",\"ScId\",\"I\",\"J\",\"K\",\"L\",\"M\","
			"\"CRdTm\",\"undCRdTm\",\"Ons\",\"Offs\",\"DepVol\",\"probStop\",\"depDelay\", "  
			"\"arrDelay\",\"dwlDelay\",\"rideDelay\",\"PVal\",\"AVal\",\"WkTmOns\",\"WkTmOffs\",\"WalkCost\","
			"\"RideCost\",\"OperCost\",\"TCost\",\"dpCnt\",\"stopDelay\" ,\"StOrdr\",\"CRdTmC\") "  
			" Values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,"
					 "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);" ;
	
	string 	strDPSQLInserTbl = "Insert into " + strtblName + "  " + strDPSQLInserTblDef ;
	
	// call insert table routine to populate the stop table
	//	blnInsert = inSpaTblStop(strSQLInserTbl,pTransOutDb,*pstop,stops,tblTrip, keyFld , intSRID,logfile);
	// make the trip dpmap
	ip = 0;//kStop1.tripId();

//	tripDPsmap = dpTripWalkRideMap(sdptsidpmmap,tripDPsmap,ip,dpstop,tDPstp,gcost);
	blnIns = inSpaTblDPStop(strDPSQLInserTbl,outdb,mDP1,o1,dpstop,g1,v1);
	if (blnIns) {

	}
 // write the path trace to a file
	outFN = f1 + tblDP + "MPDPathTrace.out";
	outfile.open(outFN.c_str(), ios::out  );

	{// write the dp result data to a file
		writeTextData(mmapIJDPK,strijklm,j,outfile,"Stage\tK");
		writeTextData(mmapIJDPRes,strijklm,fcost,outfile,"KeyPath\tOptimalCost");
		writeTextData(mmapStopDPRes,j,strijklm,outfile,"Stage\tKeyPath");
		writeTextData(mmapIJDPK,strijklm,j,outfile,"i-j\tk");
		writeTextData(mmapDPredK,j,i,outfile,"i\tj");
	}
			
	if (outfile.is_open()) {
		 outfile.close();
	}
	outfile.clear();
	return mDP1;
 }


template <typename i, typename j, typename e>
bool listmember (i& i1, j& j1, e& elist)
{
	i i0=0;
	bool blnFix=false;
	for (i0=i1+1;i0<j1;i0++) { // i & j
		if ((elist.find(i0) != elist.end() ) ) {// a fixed stop is skipped in this l-m pair 
			blnFix = true;
		} 
	}
    return blnFix;
}



// backward dp run - inputs multimap with dp objects, 
//m - multimap of dp objects, o - an initialized dp object 
// n - number of stops in the dp set, M - the maximum unit of boundary stop  
// s - input stop file name, f - output stop file name character array
template <typename l, typename m, typename o,typename n, typename  M,typename s,typename f,typename t,typename x,typename g,typename p,typename d,typename r,typename u,typename e,typename v>
l& dpOptimalMultiPd5db(l& mDP1, m& m1, o& o1,n& n1, M& M1, s& tblDP,f& f1,t& t1, x& x1,g& g1,p& p1,d& outdb,r& dpstop,u& stop0,e& elist, v& v1)
{
	m::iterator mit;
	typedef pair <s, o> o_Pair; 
	mit = m1.begin();
	long i=0,j=0,k=0,l=0,m=0,i1=0,i2=0;
	int M = 0, D=0, q = 0, intSRID=0,ip=0;
	bool blnCreate=false, blnIns;
	//char outfilename[ MaxStrLen +1]=""; // data output file 
	string strijklm="", strijkl="",strjklm="",strjkl="",strklm="",strijk="",str1="",outFNm="",stRW="", strtblName="";
    //char  ext[50]="\0"; // file extension
    //char* fext=ext; // file extension
// at stop n , f*(n) = 0
	double icost=0, fcost=0,fcostar=0;
	M = g1->get_maxskip() + 1;

	D = g1->get_dpdimension();
	mmapIJDPK.clear();
	mmapIJDPRes.clear();
	mmapStopDPRes.clear();
	mapDPRes.clear();
	mmapICostDPIJ.clear();
	mmapDPIC.clear();
	mmapIJDPResICost.clear(); 
	mmapIJDPK.clear(); // predecessor  
	mmapStopDPRes.clear();
	mmapDPredK.clear();


// run dp routine
	f outFN = "";
	for (j=n1;j>=n1-M1;j--)
	{
		strjklm =  to_string(j) + *t1 + to_string(n1) + *t1 + to_string(n1) + *t1 + to_string(n1);
		mmapIJDPK.insert(strlng_Pair(strjklm,n1)); // predecessor  
		mmapIJDPRes.insert(strdbl_Pair(strjklm,fcostar)); 
		mmapStopDPRes.insert(lngstr_Pair(n1,strjklm));
		mapDPRes.insert(lngdbl_Pair(n1,fcostar));
	}
	stRW = "D" + to_string<short>(D);
	stRW.append("W" + to_string<double> ((int) (g1->get_walkcost()*10)));
	stRW.append("R" + to_string<double> ((int) (g1->get_ridecost()*10)));

	outFN = f1 + tblDP + "_MPDPMissIC5d.txt";
	//fext[0]='\0';
	//strcpy(fext, strext.c_str());
	//fileName(f1,fext,outfilename,q,-1);
	ofstream missfile( outFN.c_str(), ios::out | ios::trunc );

	mDP1.clear();
	for (k=n1;k>=1;k--)
	{
		for (j=k;j>=k-M1 && j>=1 ;j--)
		{
			for (i=j;i>=j-M1 && i>=1;i--)
			{
				for (l=k;l<=k+M1 && l<=n1;l++)
				{
					strijkl = to_string(i) + *t1 + to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) ;
					for (m=l;m<=l+M1 && m<=n1;m++)
					{
						strjklm = to_string(j) + *t1 + to_string(k) + *t1 + to_string(l) + *t1 + to_string(m);
						strijklm = to_string(i) + *t1 + strjklm;
						mit=m1.find(strijklm);
						if (mit!=m1.end()) 
						{
							o1 = mit->second;
							icost = mit->second.get_tstop().get_TCost();
							if (icost<0) { 
								icost=inf;
							}
							mmapKLMDPRit = mmapIJDPRes.find(strjklm);
							if ( mmapKLMDPRit!=mmapIJDPRes.end() )
							{
								fcost = icost + mmapKLMDPRit->second;
							}
							else
							{
								if (m==n1 && l==n1 && k==n1) {
									fcost = icost;
								} else
								{ 
									fcost = icost + inf;
								}
							}
							mmapICostDPIJ.insert(dblstr_Pair(fcost,strijkl));
							mmapDPIC.insert(dblng_Pair(fcost,m));
							mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
						}
						else
						{	
							if (i==j||j==k||k==l||l==m) {
								continue;
							} else {
								cout<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
								missfile<<"key " <<strijklm<<" does not have an immediate cost"<<endl;
							}
						}
					} // i - loop
					mmapICDPIJRit=mmapICostDPIJ.begin();
					mmapDPICit=mmapDPIC.begin();
					if ( mmapICDPIJRit!=mmapICostDPIJ.end() )
					{
						fcostar =mmapICDPIJRit->first;
						if ( mmapDPICit !=mmapDPIC.end()) {
							m =mmapDPICit->second;
						}
						strijkl =mmapICDPIJRit->second;
						strijk = strijkl.substr(0,strijkl.find_last_of(*t1));
						mmapIJDPK.insert(strlng_Pair(strijkl,m)); // predecessor  
						mmapStopDPRes.insert(lngstr_Pair(m,strijkl));
						mmapDPredK.insert(lng_Pair(k,m));
						mmapIJDPRes.insert(strdbl_Pair(strijkl,fcostar));

					}
						mmapICostDPIJ.clear();
						mmapIJDPResICost.clear();
						mmapDPIC.clear();
				} // j - loop
			} // k - loop
		} // l loop

	} // m loop

	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();


 // trace the optimal path from the above result
// write the resultant dp path
	outFN = f1 + tblDP + "_MPDPResult.txt";
	outfile.open( outFN.c_str(), ios::out );
	o1.serializetexthdr(outfile);
	outFN = f1 + tblDP + "_MPDPTrace.txt";
    //fileName(f1,fext,outfilename,q,-1);
	missfile.open(outFN.c_str(), ios::out  );
	o1.serializetexthdr(missfile);
				// write the dp result data to a file
//	string txthdr =  "sKey|i|j|k|l|m|dpkey|CRdTm|"; 
//	txthdr.append("Ons|Offs|WalkCost|RideCost|OperCost|TCost");
//	outfile<<txthdr<<endl;
//					writeTextObjectData(mDP1,o1,strijklm,outfile,txthdr);
//					writeTextObjectData(mStop,ts,strijklm,outfile,txthdr);

 // find the from the i,j,k & l values from k+1 to k+M1  (i.e. "1_1_1_2 to 1_1_1_4") key values and 
 //	pick the minimum of these values as the starting point then find m and recurse  
	i=1; 
	j=1;
	k=1;
	size_t pos;
	mmapIJDPKit=mmapIJDPK.begin();
    for (l=1;l<=M1+1;l++) 
	{
		strijkl =(to_string(i) + *t1 + to_string(j)+ *t1 + to_string(k)+ *t1 + to_string(l));
		mmapKLMDPRit = mmapIJDPRes.find(strijkl);
		if (mmapKLMDPRit != mmapIJDPRes.end()) {
			mmapICostDPIJ.insert(dblstr_Pair(mmapKLMDPRit->second,strijkl));
		}
	}
	mmapICDPIJRit = mmapICostDPIJ.begin();
	if (mmapICDPIJRit != mmapICostDPIJ.end())
	{
		strijkl=mmapICDPIJRit->second;
		mmapIJDPKit=mmapIJDPK.find(strijkl);
		while (mmapIJDPKit!=mmapIJDPK.end()) 
		{
			m=mmapIJDPKit->second;
			mmapIJDPResICost.insert(strlng_Pair(strijkl,m));
			strijklm =  strijkl + *t1+ to_string(m);
			mit = m1.find(strijklm);
			if (mit!=m1.end()) {
				outfile<<mit->second;
				missfile<<mit->second;
				o1 = mit->second;
				mDP1.insert (o_Pair(strijkl,o1));
	//			mStop.insert (tdPair(strijkl,ts));
			}
			if (strijkl == (to_string(n1) + *t1 + to_string(n1)+ *t1 + to_string(n1)+ *t1 + to_string(n1)))
			{ 
				break;
			}
			strijkl = getSubstr(strijkl,*t1) + *t1 + to_string(m) ;
			mmapIJDPKit = mmapIJDPK.find(strijkl);
		}
	} else {
		cout<<"key " <<strijkl<<" does not have an immediate cost"<<endl;
		missfile<<"key " <<strijkl<<" does not have an immediate cost"<<endl;
	}

	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	if (missfile.is_open()) {
		 missfile.close();
	}
	missfile.clear();

	outFN = f1 + tblDP + "MPDPaTrace.txt";
	outfile.open(outFN.c_str(), ios::out  );
	outfile<<"Summary of DP Results "<<endl<<"DPKey\t";
	o1.serializetexthdr(outfile);
// write the dp result data to a file
	writeTextData(mDP1,strijklm,o1,outfile,"");
	if (outfile.is_open()) {
		 outfile.close();
	}
	outfile.clear();

// write the dp result to a table
// create result table of Stops 
	intSRID = x1.get_srid();

	str1 =	"( DPKey Text,StOrdr Long, StopId Long,StopName Text, Walk Integer, Ride Integer, " 
		"ScId integer, I integer, J integer, K integer, L integer, M integer,"
		" CRdTm Double,undCRdTm Double,CRdTmC Double, Ons Double,Offs Double,DepVol Double,"
		"probStop Double,depDelay Double, arrDelay Double,"
		"stopDelay Double,dwlDelay Double,rideDelay Double,PVal Double,"
		"AVal Double, WkTmOns Double,WkTmOffs Double,WalkCost Double,RideCost Double,"
		"OperCost Double, TCost Double , dpCnt long  ); ";

		//strtblName = x1.get_tbltrip() + stRW + "_MPDP" ; 
		strtblName = tblDP + "_DPTrace" ; 
	
	replace(strtblName.begin(),strtblName.end(),' ','_');
	ReplaceAll2(strtblName,"__","_");

	blnCreate = createSpaTbl(str1,strtblName,outdb,intSRID,outfile);

	if (!blnCreate) {
	// exit since the DP Result table could not be created 
		outfile << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
		cout << "Could not create the DP result table" <<strtblName<<" . Exiting !"<<endl;
		cout << "Press enter to exit!"<<endl;
		cin>>strtblName;
		exit (0);
		
	}	

	
	string strDPSQLInserTblDef =  "(\"DPKey\",\"StopId\",\"StopName\",\"Walk\",\"Ride\",\"ScId\",\"I\",\"J\",\"K\",\"L\",\"M\","
			"\"CRdTm\",\"undCRdTm\",\"Ons\",\"Offs\",\"DepVol\",\"probStop\",\"depDelay\", "  
			"\"arrDelay\",\"dwlDelay\",\"rideDelay\",\"PVal\",\"AVal\",\"WkTmOns\",\"WkTmOffs\",\"WalkCost\","
			"\"RideCost\",\"OperCost\",\"TCost\",\"dpCnt\",\"stopDelay\" ,\"StOrdr\",\"CRdTmC\") "  
			" Values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,"
					 "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);" ;
	
	string 	strDPSQLInserTbl = "Insert into " + strtblName + "  " + strDPSQLInserTblDef ;
	
	// call insert table routine to populate the stop table
	//	blnInsert = inSpaTblStop(strSQLInserTbl,pTransOutDb,*pstop,stops,tblTrip, keyFld , intSRID,logfile);
	// make the trip dpmap
	ip = 0;//kStop1.tripId();

//	tripDPsmap = dpTripWalkRideMap(sdptsidpmmap,tripDPsmap,ip,dpstop,tDPstp,gcost);
	blnIns = inSpaTblDPStop(strDPSQLInserTbl,outdb,mDP1,o1,dpstop,g1,v1);
	if (blnIns) {

	}
 // write the path trace to a file
	outFN = f1 + tblDP + "MPDPathTrace.out";
	outfile.open(outFN.c_str(), ios::out  );

	{// write the dp result data to a file
		writeTextData(mmapIJDPK,strijklm,j,outfile,"Stage\tK");
		writeTextData(mmapIJDPRes,strijklm,fcost,outfile,"KeyPath\tOptimalCost");
		writeTextData(mmapStopDPRes,j,strijklm,outfile,"Stage\tKeyPath");
		writeTextData(mmapIJDPK,strijklm,j,outfile,"i-j\tk");
		writeTextData(mmapDPredK,j,i,outfile,"i\tj");
	}
			
	if (outfile.is_open()) {
		 outfile.close();
	}
	outfile.clear();
	return mDP1;
 }



#endif // DP_OPTIMIZATION_H

#ifndef NETWORKPROCESSING_H
#define NETWORKPROCESSING_H
template <typename m, typename  o,  typename  f >
m& matchStop2Edge(m& m1, o& o1,f& f1 )
{
long id=0,eid=0;
double vcost=0,xcost=0;
	m::interator mit;
	mit=m1.end();
	while (mit!=m1.end) 
	{
		eid= mit->first;
		o1= mit->second;
		
	 // check if this edge is the link for the voronoi-generator (bus stop)
		id = o1->get_orig();
		if (id > 0 ) 
		{ // insert the begining and ending vertices of this edge into the heap set
		
		  // get the stop from the stop list and update the edge and vertex properties
			tsidmit = tsidmap.find(id);
			pstop = &(tsidmit->second);
			pstop->set_Edgeid(eid);
			mmapStopEdge.insert(lng_Pair(id,eid));
			mmapEdgeStop.insert(lng_Pair(eid,id));
			vcost = pstop->get_CRdTm();
			xcost = pstop->get_CRdTmE();
			// check to see if the direction of the edge is forward or reverse using the start cost (scost) & end cost (tcost)
            forwd = o1->get_scost()<=o1->get_tcost();
			
		//	double d =  (rand() / ((double)RAND_MAX));
		//	o1->set_palong(d);
			if (forwd) {
				o1->set_scost(vcost+d*o1->get_cost());
				o1->set_tcost(vcost+(1-d)*o1->get_cost());
	   			vx1->set_origin_vals(-1,o1->get_lbl(),o1->get_orig(),o1->get_scost(),o1->get_scost());
			}
			else
			{
				o1->set_tcost(vcost+d*o1->get_cost());
				o1->set_scost(vcost+(1-d)*o1->get_cost());
	   			vx1->set_origin_vals(-1,o1->get_lbl(),o1->get_orig(),o1->get_tcost(),o1->get_tcost());
			}
			vcost = vx1->get_cost(); // sort multi-map using cost for the vertex
			mmapSver.insert(dblng_Pair(vcost,v1));
//			mmapSver.insert(dblng_Pair(vcost,eid));
			mmapVx0.insert(dblng_Pair((vcost+d*o1->get_cost()),v1));
			mmapVx1.insert(dblng_Pair((xcost+d*o1->get_cost()),v1));
			mmapVx0R.insert(lngdbl_Pair(v1,(vcost+d*o1->get_cost())));
			mmapVx1R.insert(lngdbl_Pair(v1,(xcost+d*o1->get_cost())));
			mmapVx0.insert(dblng_Pair((vcost+(1-d)*o1->get_cost()),v2));
			mmapVx1.insert(dblng_Pair((xcost+(1-d)*o1->get_cost()),v2));
			mmapVx0R.insert(lngdbl_Pair(v2,(vcost+(1-d)*o1->get_cost())));
			mmapVx1R.insert(lngdbl_Pair(v2,(xcost+(1-d)*o1->get_cost())));

			mmapStopVx.insert(lng_Pair(id,v1));
			mmapStopVx.insert(lng_Pair(id,v2));
			mmapVxStop.insert(lng_Pair(v1,id));
			mmapVxStop.insert(lng_Pair(v2,id));
			if ( forwd ) { //(roundupx(vx2->get_cost(),3)) == (roundupx((vcost+o1->get_cost()),3))
	   			vx2->set_origin_vals(-1,o1->get_lbl(),o1->get_orig(),o1->get_tcost(),o1->get_tcost());
				mmapSver.insert(dblng_Pair((vx2->get_cost()),v2));
			}
			else  {
	   			vx2->set_origin_vals(-1,o1->get_lbl(),o1->get_orig(),o1->get_scost(),o1->get_scost());
				mmapSver.insert(dblng_Pair((vx2->get_cost()),v2));
			}
		}
	}
}

//m& funCoStopVxEdge(m& m1,n& n1, o& o1, s& s1, p& p1,x& mapvert,v& vx1, f& f1, mmaplng& pmmStopEdge, mmaplng& pmmEdgeStop, mmapdblng& pmmSver, 
//mmapdblng& pmmVx0, mmapdblng& pmmVx1, mmapdblng& pmmVx0R, mmapdblng& pmmVx1R, mmaplng& pmmStopVx,mmaplng& pmmVxStop)

// Assign Edge Costs using stops including multi-part (such as polylines) edges sum the cost of the edge 

template < typename m, typename n, typename  o,  typename  s, typename  p,typename  x, 
typename  v, typename q, typename r, typename  f>
m& funCoStopVxEdge(m& m1,n& n1, o& o1,s& s1, p& p1,x& mapvert,v& vx, n& mmapVxStop,n& mmapStopVx,
	  n& mmapStopEdge,n& mmapEdgeStop,q& mmapSver,q& mmapVx0,q& mmapVx1,r& mmapVx0R,r& mmapVx1R,
	  x& mapvertor0,x& mapvertor1,f& f1)
{
	long id=0,eid=0,eoid,v1=0,v2=0;
	double vcost=0,xcost=0,posalong=0,xdist=0,ecost=0,epalong=0;
	bool bnEdge=true;
	m::iterator mit; // edge iterator
	m::iterator mit2; 
	m m2;  // edge map set sorted by ending vertex id
	m m3;  // edge map set
	v Vx1a; // vertex
	v Vx2a;
	v Vx1b;
	v Vx2b;
	v* pVx1;
	v* pVx2;

	typedef n::iterator nit;
	nit nit1;
	s::iterator sit;
	o* pO;
	pO = &o1; // edge
	o o2; // edge
	typedef pair<long,o> obj_Pair;
	typedef pair<long,v> ovx_pair;

	mmapVxStop.clear();
	mmapStopVx.clear();
	mmapStopEdge.clear();
	mmapEdgeStop.clear();
	mmapVx0.clear();
	mmapVx1.clear();
	mmapVx0R.clear();
	mmapVx1R.clear();
	mapvertor0.clear(); 
	mapvertor1.clear();
	mmapSver.clear();

	sit=s1.begin();
	for(sit = s1.begin(); sit != s1.end(); ++sit) 
	{
		bnEdge=false;
		id= sit->first;
		p1= sit->second;
		eid = p1.get_Edgeid();
		posalong = p1.get_posalong();
	    pair<nit, nit> objrange = n1.equal_range(eid);
		size_t j = distance(objrange.first,objrange.second);
		if (j>1) { // polyline edge having multiple links
			m2.clear();
			for (nit1=objrange.first; nit1!=objrange.second;nit1++)
			{
			// for multi-part (such as polyline edges) edges sum the cost of the edge pieces to determine the two vertices where the stop 
			// is located. If j = 1 then there is no need to do this
				eoid = nit1->second;
				mit = m1.find(eoid);
				if (mit !=m1.end()) 
				{ // insert the edges by using the ending vertex id of each piece of edge into a set
						pO= &(mit->second);
						v1 = pO->get_toid();
						m2.insert(obj_Pair(v1,*pO));
				}
			}
			for (mit=m2.begin();mit!=m2.end();++mit)
			{// find the edge piece without a from id key, that will be the beginning of the edges
				mit2 = m2.find(mit->second.get_frid());
				if ( mit2== m2.end())
				{ // if the edge from id is not found among the set keys then it must be the beginning edge  
					o2= mit->second;
					break;
				}
			}
			m3.clear();
			for (mit=m2.begin();mit!=m2.end();++mit)
			{// use the beginning vertex as a key to create a set of edges emanating out of the vertex
				m3.insert(obj_Pair(mit->second.get_frid(),mit->second));
			}
			m2.clear();
			m2 = sortEdgeObj(m3,o2,m2,id);
			xcost = 0;
			for (mit=m2.begin();mit!=m2.end();++mit)
			{
				pO= &(mit->second);
				xcost = xcost + pO->get_cost();
			}
			mit = m2.begin();
			epalong = 0;
			while (mit!=m2.end()) 
			{
				pO= &(mit->second);
				epalong = epalong + pO->get_cost()/xcost;
				pO->set_palong( epalong );
				if (posalong<epalong) { // this edge is the stop link
					pO->set_toStop(true);
					epalong = (pO->get_cost()-(epalong-posalong)*xcost)/pO->get_cost();
					posalong = epalong;
					bnEdge=true;
					break;
				}
				mit++;
			}
		} else if (j==1) { // if a unique edge is found in the edge map
			nit1= objrange.first;
			eoid = nit1->second;			
			mit = m1.find(eoid);
			if (mit !=m1.end()) 
			{ 
				id= mit->first;
				pO = &(mit->second);
				bnEdge=true;
			} // if an edge is found in the edge map
		} else {
			bnEdge=false;
			cout << "There is no edge id - serial id multi-map key with " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
			f1 << "There is no edge id - serial id multi-map key with " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
		}

		if (bnEdge) 
		{

			// get the vertices from the vert list to append it to the current edge
			v1 = pO->get_frid();
			v2 = pO->get_toid();
			vcost = p1.get_CRdTm();
			xcost = p1.get_CRdTmE();
			// get the stop from the stop list and update the edge and vertex properties
			id = p1.get_id();
			pO->set_orig(id);
			pO->set_lbl(-1);
			pO->set_dirn(1);
			mmapStopEdge.insert(lng_Pair(id,eoid));
			mmapEdgeStop.insert(lng_Pair(eoid,id));
			pO->set_palong(posalong);
			vxmap_Iter = mapvert.find(v1);
			xdist = posalong * pO->get_cost();
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				pVx1 =  &(vxmap_Iter->second);
				pO->set_scost(vcost+xdist);
				pO->set_tcost((vcost+pO->get_cost()-xdist));
	 			mmapSver.insert(dblng_Pair(vcost,v1));
				mmapVx0.insert(dblng_Pair((vcost+xdist),v1));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v1));
				mmapVx0R.insert(lngdbl_Pair(v1,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v1,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v1));
				mmapVxStop.insert(lng_Pair(v1,id));
				pVx1->set_origin_vals(-1,-1,id,vcost+xdist,vcost+xdist);
				Vx1a = *pVx1;
				mapvertor0.insert(ovx_pair(v1,Vx1a));
				pVx1->set_origin_vals(-1,-1,id,xcost+xdist,xcost+xdist);
				Vx1b = *pVx1;
				mapvertor1.insert(ovx_pair(v1,Vx1b));
			} else		
			{
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<< endl;
			}  // 4i

			vxmap_Iter = mapvert.find(v2);
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				pVx2 =  &(vxmap_Iter->second);
				xdist = (1-posalong) * pO->get_cost();
				mmapVx0.insert(dblng_Pair((vcost+xdist),v2));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v2));
				mmapVx0R.insert(lngdbl_Pair(v2,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v2,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v2));
				mmapVxStop.insert(lng_Pair(v2,id));
				pVx2->set_origin_vals(-1,-1,id,vcost+xdist,vcost+xdist);
				Vx2a = *pVx2;
				mapvertor0.insert(ovx_pair(v2,Vx2a));
				pVx2->set_origin_vals(-1,-1,id,xcost+xdist,xcost+xdist);
				Vx2b = *pVx2;
				mmapSver.insert(dblng_Pair(pVx2->get_cost(),v2));
				mapvertor1.insert(ovx_pair(v2,Vx2b));
			} else {
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<< endl;
			}  // 4i
		}// if the edge is in this link 
		else 
		{
			cout << "There is no edge with a key of " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
			f1 << "There is no edge with a key of " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
		}  // if an edge is found in the edge map
	} // loop through the stops 4i if the edge id is in the edge-object map
	return m1;
} // end routine


// Using the stop set, fix the vertex and edge costs that are the entry points into the network.  
template < typename m, typename n, typename  o,  typename  s, typename  p,typename  x, 
typename  v, typename q, typename r, typename  f>
m& fixVxEdgeCoStop(m& m1, n& n1, o& o1,s& s1, p& p1,x& mapvert,v& vx, n& mmapVxStop,n& mmapStopVx,
	  n& mmapStopEdge,n& mmapEdgeStop,q& mmapSver,q& mmapVx0,q& mmapVx1,r& mmapVx0R,r& mmapVx1R,
	  x& mapvertor0,x& mapvertor1,f& f1)
{
	long id=0,eid=0,eoid,v1=0,v2=0;
	double vcost=0,xcost=0,posalong=0,xdist=0,ecost=0,epalong=0;
	bool bnEdge=true;
	m::iterator mit; // edge iterator
	m::iterator mit2; 
	m m2;
	m m3;
	v Vx1a; // vertex
	v Vx2a;
	v Vx1b;
	v Vx2b;
	v* pVx1;
	v* pVx2;
	typedef n::iterator nit;
	nit nit1;

	typedef s::iterator sit;
	sit sit1;
	o* pO;
	pO = &o1; // edge
	o o2; // edge
	typedef pair<long,o> obj_Pair;
	typedef pair<long,v> ovx_pair;

	mmapVxStop.clear();
	mmapStopVx.clear();
	mmapStopEdge.clear();
	mmapEdgeStop.clear();
	mmapVx0.clear();
	mmapVx1.clear();
	mmapVx0R.clear();
	mmapVx1R.clear();
	mapvertor0.clear();
	mapvertor1.clear();
	mmapSver.clear();

	sit1=s1.begin();
	for(sit1 = s1.begin(); sit1 != s1.end(); ++sit1) 
	{
		bnEdge=false;
		id= sit1->first;
		p1= sit1->second;
		eid = p1.get_Edgeid();
		posalong = p1.get_posalong();
	    pair<nit, nit> objrange = n1.equal_range(eid);
		size_t j = distance(objrange.first,objrange.second);
		if (j>=1) { 
			nit1= objrange.first;
			eoid = nit1->second;			
			mit = m1.find(eoid);
			if (mit !=m1.end()) 
			{ 
				pO = &(mit->second);
				bnEdge=true;
			} else {
				bnEdge=false;
				cout << "There is no edge with a key of " << eid<< 
					" for the stop " << id << " in this set. "<<endl;
				f1 << "There is no edge with a key of " << eid<< 
					" for the stop " << id << " in this set. "<<endl;
			}  // if an edge is found in the edge map
		} else {// if an edge is missing in the edge map
			bnEdge=false;
			cout << "There is no edge id - serial id multi-map key with " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
			f1 << "There is no edge id - serial id multi-map key with " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
		}  // if an edge is found in the edge map

		if (bnEdge) 
		{
			// get the vertices from the vert list to append it to the current edge
			v1 = pO->get_frid();
			v2 = pO->get_toid();
			vcost = p1.get_CRdTm();
			xcost = p1.get_CRdTmE();
			// get the stop from the stop list and update the edge and vertex properties
			id = p1.get_id();
			pO->set_orig(id);
			pO->set_lbl(-1);
			pO->set_dirn(1);
			mmapStopEdge.insert(lng_Pair(id,eid));
			mmapEdgeStop.insert(lng_Pair(eid,id));
			pO->set_palong(posalong);
			vxmap_Iter = mapvert.find(v1);
			xdist = posalong * pO->get_cost();
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				pVx1 =  &(vxmap_Iter->second);
				pO->set_scost(vcost+xdist);
				pO->set_tcost((vcost+pO->get_cost()-xdist));
	 			mmapSver.insert(dblng_Pair(vcost,v1));
				mmapVx0.insert(dblng_Pair((vcost+xdist),v1));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v1));
				mmapVx0R.insert(lngdbl_Pair(v1,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v1,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v1));
				mmapVxStop.insert(lng_Pair(v1,id));
				pVx1->set_origin_vals(-1,-1,id,vcost+xdist,vcost+xdist);
				Vx1a = *pVx1;
				mapvertor0.insert(ovx_pair(v1,Vx1a));
				pVx1->set_origin_vals(-1,-1,id,xcost+xdist,xcost+xdist);
				Vx1b = *pVx1;
				mapvertor1.insert(ovx_pair(v1,Vx1b));
			} else		
			{
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<< endl;
			}  // 4i

			vxmap_Iter = mapvert.find(v2);
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				pVx2 =  &(vxmap_Iter->second);
				xdist = (1-posalong) * pO->get_cost();
				mmapVx0.insert(dblng_Pair((vcost+xdist),v2));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v2));
				mmapVx0R.insert(lngdbl_Pair(v2,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v2,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v2));
				mmapVxStop.insert(lng_Pair(v2,id));
				pVx2->set_origin_vals(-1,-1,id,vcost+xdist,vcost+xdist);
				Vx2a = *pVx2;
				mapvertor0.insert(ovx_pair(v2,Vx2a));
				pVx2->set_origin_vals(-1,-1,id,xcost+xdist,xcost+xdist);
				Vx2b = *pVx2;
				mmapSver.insert(dblng_Pair(pVx2->get_cost(),v2));
				mapvertor1.insert(ovx_pair(v2,Vx2b));
			} else {
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<< endl;
			}  // 4i
		}// if the edge is in this link 
	} // loop through the stops 4i if the edge id is in the edge-object map
	return m1;
} // end routine

// maped,mapEidOid,ev,tsidmap,*pstop,mapvert,*pVx1,mmapVxStop,mmapStopVx,mmapStopEdge,
//	mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,mmapVx0R,mmapVx1R,mapvertor0,mapvertor1,outedgefile

template < typename m, typename n, typename  o,  typename  s, typename  p,typename  x, 
typename  v, typename q, typename r, typename  f>
m& funCoStopVxEdge2(m& m1,n& n1, o& o1,s& s1, p& p1,x& mapvert,v& vx, n& mmapVxStop,n& mmapStopVx,
	  n& mmapStopEdge,n& mmapEdgeStop,q& mmapSver,q& mmapVx0,q& mmapVx1,r& mmapVx0R,r& mmapVx1R,
	  x& mapvertor0,x& mapvertor1,f& f1)
{
	long id=0,eid=0,eoid,v1=0,v2=0;
	double vcost=0,xcost=0,posalong=0,xdist=0,ecost=0,epalong=0;
	bool bnEdge=true;
	m::iterator mit; // edge iterator
	m::iterator mit2; // edge iterator
	m m2;
	m m3;
	v Vx1a; // vertex
	v Vx2a;
	v Vx1b;
	v Vx2b;
	v* pVx1;
	v* pVx2;

	typedef n::iterator nit; // edge id oid iterator
	nit nit1;  // edge id oid iterator
	s::iterator sit; // stop iterator
	o* pO;
	pO = &o1; // edge
	o o2; // edge
	typedef pair<long,o> obj_Pair;
	typedef pair<long,v> ovx_pair;

	mmapVxStop.clear();
	mmapStopVx.clear();
	mmapStopEdge.clear();
	mmapEdgeStop.clear();
	mmapVx0.clear();
	mmapVx1.clear();
	mmapVx0R.clear();
	mmapVx1R.clear();
	mapvertor0.clear();
	mapvertor1.clear();
	mmapSver.clear();

	sit=s1.begin();
	for(sit = s1.begin(); sit != s1.end(); ++sit) 
	{
		bnEdge=false;
		id= sit->first;
		p1= sit->second;
		eid = p1.get_Edgeid();
		posalong = p1.get_posalong();
	    pair<nit, nit> objrange = n1.equal_range(eid);
		size_t j = distance(objrange.first,objrange.second);
		if (j>1) { 
			m2.clear();
			for (nit1=objrange.first; nit1!=objrange.second;nit1++)
			{
			// for multi-part (such as polylines) edges sum the cost of the edge pieces to determine where the stop 
			// is located and find out the two vertices that it belongs. If j = 1 then there is no need to do this
				eoid = nit1->second;
				mit = m1.find(eoid);
				if (mit !=m1.end()) 
				{ // insert the edges by using the ending vertex id of each piece of edge into a set
						pO= &(mit->second);
						v1 = pO->get_toid();
						m2.insert(obj_Pair(v1,*pO));
				}
			}
			for (mit=m2.begin();mit!=m2.end();++mit)
			{// find the edge piece with out a key, that will be the beginning of the edges
				mit2 = m2.find(mit->second.get_frid());
				if ( mit2== m2.end())
				{ // if the edge from id is not found among the set keys then it must be the beginning edge  
					o2= mit->second;
					break;
				}
			}
			m3.clear();
			for (mit=m2.begin();mit!=m2.end();++mit)
			{// use the beginning vertex as a key to create a set out of the one above
				m3.insert(obj_Pair(mit->second.get_frid(),mit->second));
			}
			m2.clear();
			m2 = sortEdgeObj(m3,o2,m2,id);
			xcost = 0;
			for (mit=m2.begin();mit!=m2.end();++mit)
			{
				pO= &(mit->second);
				xcost = xcost + pO->get_cost();
			}
			mit = m2.begin();
			epalong = 0;
			while (mit!=m2.end()) 
			{
				pO= &(mit->second);
				epalong = epalong + pO->get_cost()/xcost;
				pO->set_palong( epalong );
				if (posalong<epalong) { // this edge is the stop link
					pO->set_toStop(true);
					epalong = (pO->get_cost()-(epalong-posalong)*xcost)/pO->get_cost();
					posalong = epalong;
					bnEdge=true;
					break;
				}
				mit++;
			}
		} else if (j==1) {
			nit1= objrange.first;
			eoid = nit1->second;			
			mit = m1.find(eoid);
			if (mit !=m1.end()) 
			{ 
				id= mit->first;
				pO = &(mit->second);
				bnEdge=true;
			}
		} else {
			bnEdge=false;
		}

		if (bnEdge) 
		{

			// get the vertices from the vert list to append it to the current edge
			v1 = pO->get_frid();
			v2 = pO->get_toid();
			vcost = p1.get_CRdTm();
			xcost = p1.get_CRdTmE();
			// get the stop from the stop list and update the edge and vertex properties
			id = p1.get_id();
			pO->set_orig(id);
			pO->set_lbl(-1);
			pO->set_dirn(1);
			mmapStopEdge.insert(lng_Pair(id,eoid));
			mmapEdgeStop.insert(lng_Pair(eoid,id));
			pO->set_palong(posalong);
			vxmap_Iter = mapvert.find(v1);
			xdist = posalong * pO->get_cost();
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				pVx1 =  &(vxmap_Iter->second);
				pO->set_scost(vcost+xdist);
				pO->set_tcost((vcost+pO->get_cost()-xdist));
	 			mmapSver.insert(dblng_Pair(vcost,v1));
				mmapVx0.insert(dblng_Pair((vcost+xdist),v1));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v1));
				mmapVx0R.insert(lngdbl_Pair(v1,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v1,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v1));
				mmapVxStop.insert(lng_Pair(v1,id));
				pVx1->set_origin_vals(-1,-1,id,vcost+xdist,vcost+xdist);
				Vx1a = *pVx1;
				mapvertor0.insert(ovx_pair(v1,Vx1a));
				pVx1->set_origin_vals(-1,-1,id,xcost+xdist,xcost+xdist);
				Vx1b = *pVx1;
				mapvertor1.insert(ovx_pair(v1,Vx1b));
			} else		
			{
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<< endl;
			}  // 4i

			vxmap_Iter = mapvert.find(v2);
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				pVx2 =  &(vxmap_Iter->second);
				xdist = (1-posalong) * pO->get_cost();
				mmapVx0.insert(dblng_Pair((vcost+xdist),v2));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v2));
				mmapVx0R.insert(lngdbl_Pair(v2,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v2,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v2));
				mmapVxStop.insert(lng_Pair(v2,id));
				pVx2->set_origin_vals(-1,-1,id,vcost+xdist,vcost+xdist);
				Vx2a = *pVx2;
				mapvertor0.insert(ovx_pair(v2,Vx2a));
				pVx2->set_origin_vals(-1,-1,id,xcost+xdist,xcost+xdist);
				Vx2b = *pVx2;
				mmapSver.insert(dblng_Pair(pVx2->get_cost(),v2));
				mapvertor1.insert(ovx_pair(v2,Vx2b));
			} else {
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<< endl;
			}  // 4i
		}// if the edge is in this link 
		else 
		{
			cout << "There is no edge with a key of " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
			f1 << "There is no edge with a key of " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
		}  // if an edge is found in the edge map
	} // loop through the stops 4i if the edge id is in the edge-object map
	return m1;
} // end routine

// m - stop map, n - dpstop map with long id, o - dpstop map with string id, a - stop, b -dpstop
// q - iteration number , r key for n, s - key for o
template < typename t, typename n,   typename  o, typename a,typename b,typename q,typename r,typename s>
o& insertDPStops(t& t1,n& n1,o& o1,a& a1,b& b1,q& q1,r& r1, s& s1)
{
	double vcost=0;
	t::iterator tit; // stop map iterator
	typedef n::iterator nit; // dpstop map iterator
	nit nit1;
	typedef pair<r,b> r1o1Pair; // r - long key, b - dpstop
	typedef pair<s,b> r2o1Pair; // s - string key, b - dpstop
	long i=0,j=0,k=0,l=0,m=0;
	static long dpi=0;
	a a2; // stop
	a* pA1;
//	a* pA2;
	b b2; // dpstop
//	b* pb2;

	tit=t1.begin();
	i = tit->first;
	while (tit!=t1.end())
		{// store the i-j-k-l-m data in the dp store
			a1 = tit->second;
			pA1 = &a1;
			k = a1.get_id();
				if (tit == t1.begin()) {
					 i=j=k;
					} else {
						tit--;
	 					a2 = tit->second;
						j=a2.get_id();
						if (tit == t1.begin()) {
							i=j;
						} else {
							tit--;
	 						a2 = tit->second;
							i=a2.get_id();
							tit++;
						}
							tit++;
					}
					if (++tit == t1.end()) {
					 l=m=k;
					} else {
	 					a2 = tit->second;
						l=a2.get_id();
						if (++tit == t1.end()) {
							m=l;
						} else {
	 						a2 = tit->second;
							m=a2.get_id();
						}
							tit--;
					}
							tit--;
					b b2((*pA1,i,j,k,l,m,q1));
					b2.set_tstop(a1);
					b2.set_i(i);
					b2.set_j(j);
					b2.set_k(k);
					b2.set_l(l);
					b2.set_m(m);
					b2.set_q(q1);
					b2.set_dpkey(b2.makey(b2));

					n1.insert(r1o1Pair(++dpi,b2));
					o1.insert(r2o1Pair(b2.get_dpkey(),b2));
					tit++;
				}
return o1;
}


// m - stop map, n - dpstop map with long id, o - dpstop map with string id, a - stop, b -dpstop
// q - iteration number , r key for n, s - key for o , x a string to be appended to output filename
template < typename t, typename n,   typename  o, typename a,typename b,typename q,typename r,typename s>
o& insert5DPStops(t& t1,n& n1,o& o1,a& a1,b& b1,q& q1,r& r1, s& s1)
{
	double vcost=0;
	t::iterator tit; // stop map iterator
	typedef n::iterator nit; // dpstop map iterator
	nit nit1;
	typedef pair<r,b> r1o1Pair; // r - long key, b - dpstop
	typedef pair<s,b> r2o1Pair; // s - string key, b - dpstop
	long i=0,j=0,k=0,l=0,m=0;
	static long dpi=0;
	a a2; // stop
	a* pA1;
//	a* pA2;
	b b2; // dpstop
//	b* pb2;

	tit=t1.begin();
	i = tit->first;
	while (tit!=t1.end())
		{// store the i-j-k-l-m data in the dp store
			a1 = tit->second;
			pA1 = &a1;
			k = a1.get_id();
				if (tit == t1.begin()) {
					 i=j=k;
					} else {
						tit--;
	 					a2 = tit->second;
						j=a2.get_id();
						if (tit == t1.begin()) {
							i=j;
						} else {
							tit--;
	 						a2 = tit->second;
							i=a2.get_id();
							tit++;
						}
							tit++;
					}
					if (++tit == t1.end()) {
					 l=m=k;
					} else {
	 					a2 = tit->second;
						l=a2.get_id();
						if (++tit == t1.end()) {
							m=l;
						} else {
	 						a2 = tit->second;
							m=a2.get_id();
						}
							tit--;
					}
							tit--;
					b b2((*pA1,i,j,k,q1));
					b2.set_tstop(a1);
					b2.set_i(i);
					b2.set_j(j);
					b2.set_k(k);
					b2.set_l(l);
					b2.set_m(m);
					b2.set_q(q1);
					b2.set_dpkey(b2.makey(b2));

					n1.insert(r1o1Pair(++dpi,b2));
					o1.insert(r2o1Pair(b2.get_dpkey(),b2));
					tit++;
				}
return o1;
}


// m - stop map, n - dpstop map with long id, o - dpstop map with string id, a - stop, b -dpstop
// q - iteration number , r key for n, s - key for o
template < typename t, typename n,   typename  o, typename a,typename b,typename q,typename r,typename s>
o& insert3DPStops(t& t1,n& n1,o& o1,a& a1,b& b1,q& q1,r& r1, s& s1)
{
	double vcost=0;
	t::iterator tit; // stop map iterator
	typedef n::iterator nit; // dpstop map iterator
	nit nit1;
	typedef pair<r,b> r1o1Pair; // r - long key, b - dpstop
	typedef pair<s,b> r2o1Pair; // s - string key, b - dpstop
	long i=0,j=0,k=0,l=0,m=0;
	static long dpi=0;
	a a2; // stop
	a* pA1;
//	a* pA2;
	b b2; // dpstop
//	b* pb2;

	tit=t1.begin();
	i = tit->first;
	while (tit!=t1.end())
		{// store the i-j-k-l-m data in the dp store
			a1 = tit->second;
			pA1 = &a1;
			j = a1.get_id();
				if (tit == t1.begin()) {
					 i=j;
					} else {
						tit--;
	 					a2 = tit->second;
						i=a2.get_id();
							tit++;
					}
					if (++tit == t1.end()) {
					 k=j;
					} else {
	 					a2 = tit->second;
						k=a2.get_id();
						}

					b b2((*pA1,i,j,k,q1));
					b2.set_tstop(a1);
					b2.set_i(i);
					b2.set_j(j);
					b2.set_k(k);
					b2.set_q(q1);
					b2.set_dpkey(b2.makey3d(b2));

					n1.insert(r1o1Pair(++dpi,b2));
					o1.insert(r2o1Pair(b2.get_dpkey(),b2));
				}
return o1;
}

// stops - aStop - stopVx - VxOrigins - Vertex Set
template < typename m,  typename  s, typename n,   typename  v, typename  x,typename  o,typename  d,typename  e,typename  g>
m& updateVxOrigins2(m& m1,s& s1,n& n1,v& v1,v& v2,x& vx,o& o1, d& d1,e &e1,g& g1)
{
	double vcost=0;
	o o2,o3; // key value
	s s2; // stop object
	m::iterator mit; // stop iterator
	m::iterator mit2; // stop iterator
	e::iterator eit; // edge iterator
	typedef n::iterator nit;
	nit nit1;
	typedef v::iterator vit;
	vit vit1;
	vit vit2;
	typedef pair<o,x> ovxPair;

	x Vx1; // vertex
	x Vx2;
	x* pVx1;
	x* pVx2;
	mit=m1.begin();
	for(mit = m1.begin(); mit != m1.end(); ++mit) 
	{
		o1= mit->first;
		s1= mit->second; //stop obejct
		eit = e1.find(s1.get_Edgeid());
		g1 = eit->second;
		if (eit!=e1.end()) { // edge is in the set
			if (d1==0) {
				vcost = s1.get_CRdTm() + g1.get_cost()*s1.get_posalong();
			} else {
				vcost = s2.get_CRdTmE() + g1.get_cost()*s1.get_posalong();
			}
			pair<nit,nit> vertSet = n1.equal_range(o1); // find the vertices that belong to this stop
	   		size_t j = distance(vertSet.first,vertSet.second);
			if (j>=1) {
				for (nit1=vertSet.first; nit1!=vertSet.second;nit1++)
				{
    				o2 = nit1->second;
					vit1 = v1.find(o2);  // is there  a Vx in the origin set
					if (vit1 != v1.end()) {
						Vx1 = vit1->second;
						if (o1 != Vx1.get_orig()) {
							mit2 = m1.find(	Vx1.get_orig());
							if (mit2!=m1.end()) {
								s2 = mit2->second;
								if (vcost<Vx1.get_cost()) {
									Vx1.set_cost(vcost);
									Vx1.set_orig(o1);
								}
							} // stop is part of the set
							else { // Vx1's stop is not in the stop set switch it to the stop
								Vx1.set_orig(o1);
								Vx1.set_cost(vcost);
							}
						}
						pVx1=&Vx1;
						vit2 = v2.find(o2);
						Vx2 = vit2->second;
						pVx2 = &Vx2;
						Vx2 = Vx1;
					} else { // Vx1 is not in the origin vertex set
						cout << " Origin Vertex "<< o2<< " Not found for stop id "<<o1<< " !"<<endl;
						vit2 = v2.find(o2);
						Vx2 = vit2->second;
						Vx2.set_orig(o1);
						Vx2.set_cost(vcost);
						
					}
					v2.erase(vit2);
					v2.insert(ovxPair(o2,Vx2));
				} // loop over the vertices for this stop 
			} else { 
				cout << " No vertices found for stop id "<<o2<< " !"<<endl;
			}
		} else {
			cout<< " Edge id "<<s1.get_Edgeid()<<" is not not found for stop "<<s1.get_id()<<endl;
		}

	} // loop through the vertex origins map
	return m1;
} // end routine
// m - StopVx, n - VxCost, v - CostVx,  
template < typename m, typename n,   typename  t,typename  v, typename  x,typename  o>
m& updateVxOrigins0(m& m1,n& n1,t& t1,v& v1,v& v2,x& vx,o& o1)
{
	double y1=0,y2=0;
	m::iterator mit; // edge iterator
	o o2,o3;
	typedef n::iterator nit;
	nit nit1;
	typedef v::iterator vit;
	vit vit1;
	vit vit2;
	typedef t::iterator tit;
	tit tit1;
	tit tit2;
	typedef pair<o,x> ovxPair;
	x Vx1; // vertex
	x Vx2;
	x* pVx1;
	x* pVx2;
	mit=m1.begin();
	for(mit = m1.begin(); mit != m1.end(); ++mit) 
	{
		o1 = mit->first; 
		o2 = mit->second; 
		pair<nit,nit> vertSet = n1.equal_range(o2); // find the vertices that belong to this stop
	   	size_t j = distance(vertSet.first,vertSet.second);
		if (j>0) {
			y1 = 1E6;
			for (nit1=vertSet.first; nit1!=vertSet.second;nit1++)
			{
    				y2 = nit1->second;
					if (y2<y1) {
						tit1 = t1.find(y1);
						y1=y2;
						if (tit1!=t1.end()) { // erase this key
							t1.erase(tit1);
						}
						vit1 = v1.find(o2);
						if (vit1 != v1.end()) {
							Vx1 = vit1->second;
							pVx1=&Vx1;
							pVx1->set_orig(o1);
							pVx1->set_cost(y1);
							v1.erase(vit1);
							v1.insert(ovxPair(o2,Vx1));
						}
					}
						vit2 = v2.find(o2);
						if (vit2 != v2.end()) {
							Vx2 = vit2->second;
							pVx2 = &Vx2;
							Vx2 = Vx1;
							v2.erase(vit2);
							v2.insert(ovxPair(o2,Vx2));
						}
			}
		} else {
			cout << " Vertex " << o2 <<" not found for stop "<<o1<<endl;
		}
	} // loop through the vertex origins map
	return m1;
} // end routine


template < typename m, typename n,   typename  v, typename  x,typename  o>
m& updateVxOrigins(m& m1,n& n1,v& v1,v& v2,x& vx,o& o1)
{
	double vcost=0;
	m::iterator mit; // edge iterator
	typedef n::iterator nit;
	nit nit1;
	typedef v::iterator vit;
	vit vit1;
	vit vit2;
	typedef pair<o,x> ovxPair;

	x Vx1; // vertex
	x Vx2;
	x* pVx1;
	x* pVx2;
	mit=m1.begin();
	for(mit = m1.begin(); mit != m1.end(); ++mit) 
	{
		vcost= mit->first;
		o1= mit->second;
		vit1 = v1.find(o1);
		Vx1 = vit1->second;
		pVx1=&Vx1;
		vit2 = v2.find(o1);
		Vx2 = vit2->second;
		pVx2 = &Vx2;
		Vx2 = Vx1;
		v2.erase(vit2);
		v2.insert(ovxPair(o1,Vx2));
		vit2 = v2.find(o1);
	} // loop through the vertex origins map
	return m1;
} // end routine


template < typename m, typename n, typename  o,  typename  s, typename  p,typename  x, typename  v,typename  u,  typename  f>
m& upCoStopVxEdge(m& m1,n& n1, o& o1, s& s1, p& p1,x& mapvert,v& vx1,u& u1, f& f1)
{
	long id=0,eid=0,eoid,v1=0,v2=0;
	double vcost=0,xcost=0,posalong=0,xdist=0,ecost=0,epalong=0;
	bool bnEdge=true;
	m::iterator mit;
	m::iterator mit2;
	m m2;
	m m3;
//	v vx;
	v vx2;
	typedef n::iterator nit;
	nit nit1;
	s::iterator sit;
	o* pO;
	pO = &o1;
	o o2;
	typedef pair<long,o> obj_Pair;
	sit=s1.begin();
	for(sit = s1.begin(); sit != s1.end(); ++sit) 
	{
		bnEdge=false;
		id= sit->first;
		p1= sit->second;
		eid = p1.get_Edgeid();
		posalong = p1.get_posalong();
	    pair<nit, nit> objrange = n1.equal_range(eid);
		size_t j = distance(objrange.first,objrange.second);
		if (j>1) { 
			for (nit1=objrange.first; nit1!=objrange.second;nit1++)
			{
			// sum the cost of the broken edge to determine where the stop is located 
			// and find out the two vertices that it belongs. If j = 1 then there is no need to do this
				eoid = nit1->second;			
				mit = m1.find(eoid);
				if (mit !=m1.end()) 
				{ // insert the begining and ending vertices of this edge into the heap set
						pO= &(mit->second);
						v1 = pO->get_toid();
						m2.insert(obj_Pair(v1,*pO));
				}
			}
			for (mit=m2.begin();mit!=m2.end();++mit)
			{
				mit2 = m2.find(mit->second.get_frid());
				if ( mit2== m2.end())
				{
					o2= mit->second;
				}
			}
			for (mit=m2.begin();mit!=m2.end();++mit)
			{
				m3.insert(obj_Pair(mit->second.get_frid(),mit->second));
			}
			m2.clear();
			m2 = sortEdgeObj(m3,o2,m2,id);
			xcost = 0;
			for (mit=m2.begin();mit!=m2.end();++mit)
			{
				pO= &(mit->second);
				xcost = xcost + pO->get_cost();
			}
			mit = m2.begin();
			while (mit!=m2.end()) 
			{
				pO= &(mit->second);
				epalong = epalong + pO->get_cost()/xcost;
				pO->set_palong( epalong );
				if (posalong<epalong) { // this edge is the stop link
					pO->set_toStop(true);
					epalong = (pO->get_cost()-(epalong-posalong)*xcost)/pO->get_cost();
					posalong = epalong;
					bnEdge=true;
					break;
				}
				mit++;
			}
		} else if (j=1) {
			nit1= objrange.first;
			eoid = nit1->second;			
			mit = m1.find(eoid);
			if (mit !=m1.end()) 
			{ // insert the begining and ending vertices of this edge into the heap set
				id= mit->first;
				pO = &(mit->second);
				bnEdge=true;
			}
		} else {
			bnEdge=false;
		}

		if (bnEdge) 
		{

			// get the vertices from the vert list to append it to the current edge
			v1 = pO->get_frid();
			v2 = pO->get_toid();
			vcost = p1.get_CRdTm();
			xcost = p1.get_CRdTmE();
			// get the stop from the stop list and update the edge and vertex properties
			id = p1.get_id();
			pO->set_orig(id);
			pO->set_lbl(-1);
			pO->set_dirn(1);
			pO->set_palong(posalong);
			vxmap_Iter = mapvert.find(v1);
			xdist = posalong * pO->get_cost();
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				vx1 =  &(vxmap_Iter->second);
				pO->set_scost(vcost+xdist);
				pO->set_tcost((vcost+pO->get_cost()-xdist));
				vx1->set_origin_vals(-1,pO->get_lbl(),pO->get_orig(),pO->get_scost(),pO->get_scost());

			} else		
			{
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v1 << " for Edge "<<eid<< endl;
			}  // 4i

			vxmap_Iter = mapvert.find(v2);
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				vx2 =  &(vxmap_Iter->second);
				xdist = (1-posalong) * pO->get_cost();
				vx2->set_origin_vals(-1,pO->get_lbl(),pO->get_orig(),pO->get_tcost(),pO->get_tcost());
			} else		
			{
				cout << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<<endl;
				f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<< endl;
			}  // 4i
		}// if the edge is in this link 
		else 
		{
			cout << "There is no edge with a key of " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
			f1 << "There is no edge with a key of " << eid<< 
				" for the stop " << id << " in this set. "<<endl;
		}  // if an edge is found in the edge map
	} // loop through the stops 4i if the edge id is in the edge-object map
	return m1;
} // end routine

// maplngvx -m  m1 - mapvert , mmaplnged n, - n1 - mmaped, mmaplng - r, r1 - mmapved,
//  mmaplng - s, s1 - mmapEdgeStop,short - t, t1 - onoff
template <typename m,typename n,typename  r,typename  s,typename  v,typename  e,typename  t,typename  q,typename f>
n& edgeVoronoi (m& m1, n& n1,r& r1,s& s1,v& v1,e& e1, t& t1,q& q1, f& f1)
{
    v* pVx1,*pVx2;
    e* pEg;
    long o1=0,o2=0, eid=0,esid=0; //  eid - edge id, eoid - edge object id
    int j=0,ne=0,nv=0,ne1=0,ne2=0,ne3=0,ne4=0; // ne - no of edges, nv - no of vertices,# fixed, #loop, # bound, # free
    string str1;
	double vcost=0,xcost=0,ecost=0,vcost1=0;
	bool fixed,bound,free,free1,free2,forwd,bacwd,loop;
	m::iterator mit;
	n::iterator nit;
// process the edge voronoi routine now

nit = n1.begin();
while(nit != n1.end()) 
{ 
	 // get the edge
	pEg = &(nit->second);
    esid =(*nit).first;
	
	mit = m1.find(pEg->get_frid());
	if (mit!=m1.end()) {
		pVx1 = &(mit->second);
		o1 = pVx1->get_orig();
	} else {
		f1<<" Start Vertex id "<<o1<< " for edge id "<<pEg->get_id()<< " missing from Vertex Set!"<<endl;
		continue;
	}
	mit = m1.find(pEg->get_toid());
	if (mit!=m1.end()) {
		pVx2 = &(mit->second);
		o2 = pVx2->get_orig();
	} else {
		f1<<" End Vertex id "<<o2<< " for edge id "<<pEg->get_id()<< " missing from Vertex Set!"<<endl;
		continue;
	}
	fixed = o1>0 && o2>0 && o1==o2;
	bound = o1>0 && o2>0 && o1!=o2;
	free1 = o1<=0 && o2>0; 
	free2 = o2<=0 && o1>0;
	free = o1<=0 && o2<=0;
	if (fixed || bound) 
	{ // both ends are fixed including boundary edges
		forwd = (pVx1->get_id() == pEg->get_frid() && pVx2->get_id()==pEg->get_toid());
		bacwd = pVx1->get_idp() == pVx2->get_id();
		loop = pVx1->get_id() == pVx2->get_id();
		//ceil( ( num * pow( 10,x ) ) - 0.5 ) / pow( 10,x );
		vcost1=	ceil((abs(pVx1->get_cost()-pVx2->get_cost())*pow(10.0,3.0))-0.5)/pow(10.0,3.0);
		ecost = ceil((pEg->get_cost()*pow(10.0,3.0))-0.5)/pow(10.0,3.0);
					ne++;
		if (fixed)  
		{ // if it is a forward labelled arc 
			pEg->set_orig(pVx1->get_orig());
			pEg = upEdgeCostDirn(pEg,pVx1,pVx2,forwd);
		}
			if (bound) // boundary edge - beg. & end vertices are assigned to different stops
			{ // find the indifference point between the vertices 
				pEg->set_orig(-2);
				xcost = (abs(pVx1->get_cost()-pVx2->get_cost())+ pEg->get_cost())/(2*pEg->get_cost());
				pEg->set_palong(xcost);
				forwd = true;
				pEg = upEdgeCostDirn(pEg,pVx1,pVx2,forwd);
				str1.erase();
				str1= "Bdry " + to_string(pVx1->get_orig()) + " - " + to_string(pVx2->get_orig());
				pEg->set_enote( str1 );
//				pEg->show_edge(outedgebdyfile);
				ne2++;
			}
		}
		else if (free1) // fixed at v2
		{
				pEg->set_tcost(pVx2->get_cost() + pEg->get_cost());
				pEg->set_orig(pVx2->get_orig());
				pVx1->set_idp(pVx2->get_id());
				pVx1->set_orig(pVx2->get_orig());
				pVx1->set_cost(pEg->get_tcost());
				pVx1->set_tcost(pEg->get_tcost());
				pEg->set_enote( "free1" );
				ne3++;
		}
		else if (free2) // fixed at v2
		{
				pEg->set_tcost(pVx1->get_cost() + pEg->get_cost());
				pEg->set_orig(pVx1->get_orig());
				pVx2->set_idp(pVx1->get_id());
				pVx2->set_orig(pVx1->get_orig());
				pVx2->set_cost(pVx1->get_cost()+pEg->get_cost());
				pVx2->set_tcost(pVx1->get_cost());
				pEg->set_enote( "free2" );
				ne4++;
		}
		else if (free)
		{ // it spans unfixed vertices
				pEg->set_enote( "free" );
				pEg->set_orig(0);
//				pEg->show_edge(cout);
				nv++;
		}
//		pEg->show_edge(outfile);
/*		if ((pEg->get_orig()>0)||(pEg->get_orig()==-2))
		{
			pEg->show_edge(outedgefile);
     //		pEg->show_edge(cout);
		}
	*/
		nit++;
} // while map origin is present
		f1<<" End of Edge Voronoi Process "<<endl;
		f1<<" # of Fixed Edges  = "<< ne<<endl;
		f1<<" # of Loop Edges =  "<< ne1<<endl;
		f1<<" # of Boundary Edges = "<< ne2<<endl;
		f1<<" # of Free Vertex 1 Edges = "<< ne3<<endl;
		f1<<" # of Free Vertex 2 Edges = "<< ne4<<endl;
		f1<<" # of Free Edges = "<< nv<<endl;
		f1<<" Total # Edges read = "<< (nv+ne+ne1+ne3+ne4)<<endl;

		return n1;
} // end edge voronoi


//maplngvx& vxvtarjan (mmapdblng m& mmapSver,maplngvx& v& mapvert,mmaplnged& v& mapvert2,
//	mmaplng& p& mmapved, mmaplng& p& mmapVxStop,short i& onoff) 

template <typename l,typename m,typename n,typename v,typename e,typename i,typename g>
v& vxvtarx (l& l1,m& m1, n& n1,v& v1,v& v2,e& e1,i& i1,g& gc)  
{
	long ip=0;
    vertexp vx;
    vertexp* pVx1;
    pVx1 = &vx;
	edgev ev;
	edgev* pev=&ev; 
	typedef l::iterator lIter;
	lIter lIt;
	typedef m::iterator mIter;
	mIter mIt;
	m m2;
	typedef e::iterator eIter;
	eIter eIt;
	typedef v::iterator vIter;
	vIter vIt;
	typedef pair <long,vertexp> vxPair;
	typedef pair <double,long> costPair;
// prepare the starting voronoi generators in the m2 set also save them in v2  
	for(mIt = m1.begin(); mIt != m1.end(); ++mIt) {
		m2.insert(costPair(mIt->first,mIt->second));
		vIt = v1.find(mIt->second);
		if (vIt != v1.end()) {
			pVx1=&(vIt->second);
    		pVx1->set_lbl(-1);
			v2.insert(vxPair(mIt->second,*pVx1));
		}
	}
// run a Vertex Voronoi		
	mIt = m2.begin();
    while (mIt!=m2.end()) {
		vIt = v1.find(mIt->second);
		if (vIt != v1.end()) {
			pVx1=&(vIt->second);
    		pVx1->set_lbl(-1);
			pCalcVxEdges(l1,m2,v1,v2,e1,pVx1,gc);
			m2.erase(mIt);
		}
		mIt = m2.begin();
	}

//	mmapvert1 = remapVertId2Orig(v1,mmapvert1,*pVx1,ip);
return v2;
}

template <typename p,typename s,typename n,typename v,typename x,typename i,typename g,typename h>
v& vxveuclid (p& p1,s& s1, n& n1,v& v1,x& x1,i& i1,g& gc, h& blnHist )  
{
	long ip=0, origOn=0, origOff=0;
    vertexp vx;
    vertexp* pVx1 = &vx;
    tstop ts;
    tstop* pts = &ts;
	typedef s::iterator sIter;
	sIter sIt;
	s s2;
    parcel par1;
    parcel* pp = &par1;
	typedef p::iterator pIter;
	pIter pIt;
	typedef n::iterator nIter; // long-long multimap
	nIter nIt1,nIt2;
	n n2;  // long-long multimap
	typedef v::iterator vIter;  // long-vertex multimap
	vIter vIt;
	typedef x::iterator xIter;  // long-vertex multimap
	xIter xIt;
	typedef i::iterator iIter; // double-long multimap
	iIter iIt1,iIt2;
	i i2;  // double-long multimap
	typedef pair <long,vertexp> vxPair;
	typedef pair <double,long> costPair;
	double rdTm=0, rdTm2E=0,ecost=0, wkCostOn=0, wkCostOff=0;
// run a Vertex Voronoi		
// use the set l - stops as the starting voronoi generators in the m2 set also save them in v2  
	for(pIt = p1.begin(); pIt != p1.end(); ++pIt) {
		// m is the parcel set l is the stops set
		// loop over all the stops and get the minimum cost asignment for ons and offs
		pp = &pIt->second;
		for(sIt = s1.begin(); sIt != s1.end(); ++sIt) { // loop over the stops
			pts = &sIt->second;
			rdTm = pts->get_CRdTm();
			rdTm2E = pts->get_CRdTmE();
			ecost = edist(*pp,*pts)/(gc.get_walkspd()*60);
			i1.insert(costPair((rdTm+ecost),pts->get_id()));
			i2.insert(costPair((rdTm2E+ecost),pts->get_id()));
		}
// get the off minimum
		iIt1 = i1.begin();
		origOff = iIt1->second;
		sIt = s1.find(origOff);
		pts = &sIt->second;
		ecost = edist(*pp,*pts)/(gc.get_walkspd()*60);
		wkCostOff = iIt1->first - pts->get_CRdTm();
		// make the resulting vertex data
		vertexp vx1(pp->get_id(),pts->get_id(),-1,wkCostOff,origOff,origOff,-1);
		v1.insert(vxPair(pp->get_id(),vx1));
// get the on minimum
		iIt2 = i2.begin();
		origOn = iIt2->second;
		sIt = s1.find(origOn);
		pts = &sIt->second;
		ecost = edist(*pp,*pts)/(gc.get_walkspd()*60);
		wkCostOn = iIt2->first - pts->get_CRdTmE();
		// make the resulting vertex data
		vertexp vx2(pp->get_id(),pts->get_id(),-1,wkCostOn,origOn,origOn,-1);
		x1.insert(vxPair(pp->get_id(),vx2));
		// update the parcel's on and off cost as well as the on/off stop id
		pIt->second.set_origoff(origOff);
		pIt->second.set_origon(origOn);
		if (blnHist) {
			pIt->second.set_hwkoffcost(wkCostOff);
			pIt->second.set_hwkoncost(wkCostOn);
		} else {
			pIt->second.set_wkoffcost(wkCostOff);
			pIt->second.set_wkoncost(wkCostOn);
		}
		i1.clear();
		i2.clear();
	}			
	return v1;
}


template <typename m, typename v>
m& vxvor (v& pVx1,m& m1) 
{
	long ip=0;
// run a Vertex Voronoi		
	for(mmapSverit = mmapSver.begin(); mmapSverit != mmapSver.end(); ++mmapSverit) {
			mmapSver2.insert(dblng_Pair(mmapSverit->first,mmapSverit->second));
		mapverit = mapvert.find(mmapSverit->second);
		if (mapverit != mapvert.end()) {
			pVx1=&(mapverit->second);
    		pVx1->set_lbl(-1);
			mapvert2.insert(vx_pair(mmapSverit->second,*pVx1));
		}
	}
	mmapSverit = mmapSver2.begin();
    while (mmapSverit!=mmapSver2.end()) {
		mapverit = mapvert.find(mmapSverit->second);
		if (mapverit != mapvert.end()) {
			pVx1=&(mapverit->second);
    		pVx1->set_lbl(-1);
			tarjanrx(pVx1);
			mmapSver2.erase(mmapSverit);
		}
		mmapSverit = mmapSver2.begin();
	}
	mmapvert1.clear();
	mmapvert1 = remapVertId2Orig(mapvert,mmapvert1,*pVx1,ip);
return m1;
}

// map Stop - Vertex container cleanup
template <typename u, typename v, typename x, typename y>
u& copyMap(u& u1, v& v1, x& x1, y& y1)
{  
u::iterator uit;
v::iterator vit;
typedef pair <x,y> obj_pair;

	for(uit=u1.begin();uit!=u1.end();uit++)
	{
		x1 = uit->first;
		y1 = uit->second;
				v1.insert(obj_pair(x1,y1));
	}
  return v1;
}

// adjust the vertex cost by removing the ride time for offs (ride time from start to alighting stop) 
// and ons (ride time from boarding stop to the end) 
template < typename m, typename  o, typename  s, typename  p, typename  u, typename  f>
m& reduceVxToTime2WalkTime(m& m1,o& o1, s& s1, p& p1,u& u1, f& f1)
{
	long sid=0,v1=0;
	double vcost=0,xcost=0;
	typedef m::iterator mit;
	mit mit1;
	m m2;
	s::iterator sit;
	o* pO;
	pO = &o1;
	typedef pair<long,o> obj_Pair;
	sit=s1.begin();
	for(sit = s1.begin(); sit != s1.end(); ++sit) 
	{
			// get the stop from the stop list and update the edge and vertex properties
		sid= sit->first;
		p1= sit->second;
		if (u1) {
			// u1 = 1 - ons ride time to the end
			xcost = p1.get_CRdTmE();
		} else {
			// u1 = 0 offs - ride time to current stop 
			xcost = p1.get_CRdTm();
		}
	    pair<mit, mit> objrange = m1.equal_range(sid);
		size_t j = distance(objrange.first,objrange.second);
		if (j>=1) { 
			for (mit1=objrange.first; mit1!=objrange.second;mit1++)
			{
			// adjust the vertex cost by subtracting the ride time for offs and ons using the stop fields  
						pO= &(mit1->second);
						vcost = pO->get_cost();
						if (vcost >= xcost) {
							// subtract out the ride cost
							vcost -=xcost;
						} else {
							cout << "The vertex " << pO->get_id() << " cost "<<vcost << 
								" is less than " << xcost << " the cost of stop id " << sid << 
								" for "<< u1 << " direction !"<<endl;
							f1 << endl<< "The vertex " << pO->get_id() << " cost "<< vcost << 
								" is less than " << xcost << " the cost of stop id " << sid << 
								" for "<< u1 << " direction !"<<endl;
							if (vcost < 0) {
								vcost=0;
							}
						}
							pO->set_cost(vcost);
							m2.insert(obj_Pair(sid,*pO));
			}
		} else	{
				cout << "The stop id " << sid << " doesn't have any vertices assigned to it!"<<endl;
				f1 << endl<< "The stop id " << sid << " doesn't have any vertices assigned to it!"<<endl;
			}  // 4i
	} // loop through the stops 4i
	return m1;
} // end routine


#endif NETWORKPROCESSING_H

#ifndef NETWORKVERIFICATION_H
#define NETWORKVERIFICATION_H
// v - multimap vertices  
// e - multimap edges 
template <typename v, typename e> 	
void tarjan_SCC (v& v1,e& e1)
{
//Input: Graph G = (V, E), Start node v0
// index = 0                       // DFS node number counter 
// S = empty                       // An empty stack of nodes
// tarjan(v0)                      // Start a DFS at the start node

	maplngvx S0;
	vertexp n0;
	v1::iterator vit;
	n0 = vit->second;
	S0 = tarjan(n0,v1,e1,S0);

}

template <typename n,typename v, typename e>
v& tarjan (n& n0, v& v1, e& e1, v& s0)
{
	v::iterator mumvIter;
	mumvIter mumvedit;
	typedef pair <long,n> obj_pair;
	long idx=0, vid1=0,vid2=0;
	n0.set_index(idx);
	n0.set_lowlink(idx);
	idx=idx+1;
	vid1 = n0.get_id();
	s0.insert(vid1,n0);
		  // search through the edges coming out of this vertex - vid
	pair<mumvIter, mumvIter> vedrng = e1.equal_range(vid1);
	  //size_t j = distance(vedrange.first,vedrange.second);
	  for (mumvedit = vedrng.first; mumvedit!=vedrng.second;mumvedit++)
	  {
		vid2 = (*mumvedit).second;
		// find the vertex for this edge
		mumvedit = v1.find(vid2);
		if ((*mumvedit).second.get_index<=0) // vertex has been visited
		{
			tarjan(n& n0, v& v1, e& e1, v& s0);
		}
		else
		{

		}
		//size_t j = distance(edrange.first,edrange.second);
	  }
/*
//Input: Graph G = (V, E), Start node v0

procedure tarjan(v)
  v.index = index               // Set the depth index for v
  v.lowlink = index
  index = index + 1
  S.push(v)                     // Push v on the stack
  forall (v, v') in E do        // Consider successors of v 
    if (v'.index is undefined)  // Was successor v' visited? 
      tarjan(v')                // Recurse
      v.lowlink = min(v.lowlink, v'.lowlink)
    elseif (v' in S)            // Is v' on the stack?
      v.lowlink = min(v.lowlink, v'.lowlink)
  if (v.lowlink == v.index)     // Is v the root of an SCC?
    print "SCC:"
    repeat
      v' = S.pop
      print v'
    until (v' == v)
*/

}

template <typename n>
n& tarjanr (n& n0)
{
	long  vid1=0,vid2=0;
	static	long idx=0;
	typedef mmaplng::iterator mumvIter;
	mumvIter  mumvedit;
	typedef pair <long,vertexp> obj_pair;
    vertexp vx;
    vertexp* n1;
    n1 = &vx;
	n0->set_index(idx);
	n0->set_lowlink(idx);
	idx=idx+1;
	vid1 = n0->get_id();
	scc.insert(obj_pair(vid1,*n0));
    // search through the edges emanating from this vertex - vid
	  pair<mumvIter, mumvIter> vedrng = mmapV1V2.equal_range(vid1);
	  size_t j = distance(vedrng.first,vedrng.second);
/*	  while (j==0 )
	  {
        mumVSit= (vedrng.first);
		
		n0= &(mapvert.find(mumVSit->first))->second;
		n0->set_index(idx);
		n0->set_lowlink(idx);
		vid1 = n0->get_id();
		scc.clear();
		scc.insert(obj_pair(vid1,*n0));
		vedrng = mmapV1V2.equal_range(vid1);
		j = distance(vedrng.first,vedrng.second);
	  }
	  */
	  for (mumvedit = vedrng.first; mumvedit!=vedrng.second;mumvedit++)
	  {
		vid2 = (*mumvedit).second;
		// find the vertex for this edge
		vxmap_Iter = mapvert.find(vid2);
		n1 = &(vxmap_Iter->second);
		if ((n1->get_index())==inf) // vertex has not been visited
		{
			tarjanr(n1);
			n0->set_lowlink(min(n0->get_lowlink(),n1->get_lowlink()));
		}
		else
		{
			vomap_Iter =scc.find(n1->get_id()); 
			if (vomap_Iter!=scc.end()) // n1 is on the stack
			{
				n0->set_lowlink(min(n0->get_lowlink(),n1->get_lowlink()));
			}		
		}
	  }
	  if (n0->get_lowlink() == n0->get_index()) {   // Is v the root of an SCC?
			    cout <<"SCC:";
				vomap_Iter = scc.begin();
				while (vomap_Iter !=scc.end()) 
				{
					n1 = &(vomap_Iter->second);
					if (n1!=n0) 
					{
						cout <<n1->get_id()<<" - ";
					}
					else
					{
						break;
					}

					scc.erase(vomap_Iter);
					vomap_Iter=scc.begin(); // vomap_Iter++;
				}
	  }
return n0;

/*
//Input: Graph G = (V, E), Start node v0

procedure tarjan(v)
  v.index = index               // Set the depth index for v
  v.lowlink = index
  index = index + 1
  S.push(v)                     // Push v on the stack
  forall (v, v') in E do        // Consider successors of v 
    if (v'.index is undefined)  // Was successor v' visited? 
      tarjan(v')                // Recurse
      v.lowlink = min(v.lowlink, v'.lowlink)
    elseif (v' in S)            // Is v' on the stack?
      v.lowlink = min(v.lowlink, v'.lowlink)
  if (v.lowlink == v.index)     // Is v the root of an SCC?
    print "SCC:"
    repeat
      v' = S.pop
      print v'
    until (v' == v)
*/
}


template <typename n>
n& tarjanrx (n& n0)
{
	long  vid1=0,vid2=0, esid=0;
	double vcost=0,ecost=0;
	static	long idx=0, ldx=0;
	bool bnForwd=true,bnDups=true;
	typedef mmaplng::iterator mumvIter;
	mumvIter  mumvedit;
	typedef pair <long,vertexp> vxPair;
	typedef pair <long,double> costPair;
    vertexp vx;
    vertexp* n1;
	edgev ev;
	edgev* pev=&ev; 
    n1 = &vx;
	n0->set_index(idx);
	n0->set_lowlink(ldx);
	vid1 = n0->get_id();
    // search through the edges emanating from this vertex - vid
	pair<mumvIter, mumvIter> vedrng = mmapV1V2EgId.equal_range(vid1);
	size_t j = distance(vedrng.first,vedrng.second);
	for (mumvedit = vedrng.first; mumvedit!=vedrng.second;mumvedit++)
	{
		esid = (*mumvedit).second;
			// find the edge for this vertex
		mapedi=maped.find(esid);
		if (mapedi!=maped.end()) {
			// get the edge for this vertex
			pev = &mapedi->second;
			bnForwd = pev->get_frid()==vid1;
			ecost = pev->get_cost();
			if (bnForwd) {
				vid2 =  pev->get_toid();
			} else {
				vid2 =  pev->get_frid();
			}
			bnDups = vid2==vid1;
			if (!bnDups) {
			// find the vertex for this edge
				vxmap_Iter = mapvert.find(vid2);
				if (vxmap_Iter != mapvert.end()) {
					n1 = &(vxmap_Iter->second);
					if (n1->get_lbl()!=-1)
					{
						vcost = n0->get_cost() + ecost;
						if ((n1->get_cost())==inf) // vertex has not been visited
						{
							n1->set_cost(vcost);
							n1->set_idp(vid1);
							n1->set_orig(n0->get_orig());
							mmapSver2.insert(dblng_Pair(vcost,vid2));
							n1->set_index(n0->get_index()+1);
							n1->set_lowlink(n0->get_lowlink());
							idx = n1->get_index();
							ldx = n1->get_lowlink();
							mapvert2.insert(vxPair(vid2,*n1));
							//tarjanrx(n1);
						} else if (vcost<n1->get_cost()) {
							n1->set_cost(vcost);
							n1->set_idp(vid1);
							n1->set_orig(n0->get_orig());
							n1->set_index(n0->get_index()+1);
							n1->set_lowlink(n0->get_lowlink());
							vomap_Iter=mapvert2.find(vid2);
							mapvert2.erase(vomap_Iter);
							mapvert2.insert(vxPair(vid2,*n1));
							idx = n1->get_index();
							ldx = n1->get_lowlink();
						} else if ((n1->get_cost())<inf && vcost == n1->get_cost()) {
				// tie breaker may not occur in a key that is double precision
						}
					} // if n1 is labelled skip it
				} // if vertex2 is not found 
			} // if vertex is the same as the fist in case of a loop 
		} // if edge is not found in the edge set
	  }
    return n0;
}
// m - mmapdblng map - mmapSver2, v - maplngvx - mapvert & mapvert2, 
// maplnged - e - maped, l - mmaplng (mmapV1V2EOId)  
template <typename l,typename m,typename v,typename e,typename x,typename g>
x& pCalcVxEdges (l& l1,m& m1, v& v1,v& v2,e& e1,x& x1,g& gc)
{
	long  vid1=0,vid2=0, esid=0;
	double vcost=0,ecost=0;
	static	long idx=0, ldx=0;
	bool bnForwd=true,bnLoop=false;
	typedef l::iterator lIter;
	lIter lIt;
	typedef e::iterator eIter;
	eIter eIt;
	typedef v::iterator vIter;
	vIter vIt;
	typedef pair <long,vertexp> vxPair;
	typedef pair <long,double> costPair;
    vertexp vx;
    x x2;
	edgev ev;
	edgev* pEg=&ev; 
    x2 = &vx;
	x1->set_index(idx);
	x1->set_lowlink(ldx);
	vid1 = x1->get_id();
    // search through the edges emanating from this vertex - vid
	pair<lIter, lIter> vedrng = l1.equal_range(vid1);
	size_t j = distance(vedrng.first,vedrng.second);
	for (lIt = vedrng.first; lIt!=vedrng.second;lIt++)
	{
		esid = (*lIt).second;
			// find the edge for this vertex
		eIt=e1.find(esid);
		if (eIt!=e1.end()) {
			// get the edge for this vertex
			pEg = &eIt->second;
			bnForwd = pEg->get_frid()==vid1;
			ecost = pEg->get_cost();
			if (bnForwd) {
				vid2 =  pEg->get_toid();
			} else {
				vid2 =  pEg->get_frid();
			}
			bnLoop = vid2==vid1;
			if (!bnLoop) {
			// find the vertex for this edge
				vxmap_Iter = v1.find(vid2);
				if (vxmap_Iter != v1.end()) {
					x2 = &(vxmap_Iter->second);
					if (x2->get_lbl()!=-1)
					{
						vcost = x1->get_cost() + ecost;
						if ((x2->get_cost())==inf) // vertex has not been visited
						{
							x2->set_cost(vcost);
							x2->set_idp(vid1);
							x2->set_orig(x1->get_orig());
							m1.insert(dblng_Pair(vcost,vid2));
							x2->set_index(x1->get_index()+1);
							x2->set_lowlink(x1->get_lowlink());
							idx = x2->get_index();
							ldx = x2->get_lowlink();
							v2.insert(vxPair(vid2,*x2));
							pEg->set_orig(x1->get_orig());
							upEdgeCostDirn(pEg, x1,x2,bnForwd);
							//tarjanrx(x2);
						} else if (vcost<x2->get_cost()) {
							x2->set_cost(vcost);
							x2->set_idp(vid1);
							x2->set_orig(x1->get_orig());
							pEg->set_orig(x1->get_orig());
							upEdgeCostDirn(pEg,x1,x2,bnForwd);
							x2->set_index(x1->get_index()+1);
							x2->set_lowlink(x1->get_lowlink());
							vIt=v2.find(vid2);
							v2.erase(vIt);
							v2.insert(vxPair(vid2,*x2));
							idx = x2->get_index();
							ldx = x2->get_lowlink();
						} else if ((x2->get_cost())<inf && vcost == x2->get_cost()) {
				// tie breaker may not occur in a key that is double precision
						}
					} else if (pEg->get_orig()==-1) {
						if (x1->get_orig()==x2->get_orig()) {
							pEg->set_orig(x1->get_orig());	
						} // else {pEg->set_orig(-2)}
					}// if x2 is labelled skip it after checking it the origin is not yet set
				} // if vertex2 is not found 
			} else { // if vertex is the same as the fist in case of a loop 
				pEg->set_orig(x1->get_orig());
			}
		} // if edge is not found in the edge set
	  }
    return x1;
}

template <typename e,typename v, typename b>
e& upEdgeCostDirn (e& pEg,v& pVx1, v& pVx2, b& forwd)
{

		if (forwd ) {
			pEg->set_scost(pVx1->get_cost());
			pEg->set_tcost(pVx2->get_cost());
			pEg->set_dirn(1);	
		} else if (!forwd ) {
			pEg->set_dirn(-1);	
			pEg->set_scost(pVx2->get_cost());
			pEg->set_tcost(pVx1->get_cost());
		}
		return pEg;
}
template <typename e,typename v,  typename x>
e& calcEdgeBndPos (e& pEg,v& pVx1, v& pVx2, x& xPos)
{

				xPos = (abs(pVx1->get_cost()-pVx2->get_cost())+ pEg->get_cost())/(2*pEg->get_cost());
				pEg->set_palong(xPos);
		return pEg;
}

#endif // NETWORKVERIFICATION_H ///:~

#ifndef  MISCFUNCTIONS_H ///:~
#define  MISCFUNCTIONS_H ///:~

// rounding functions

template < typename o1, typename o2>
double edist(o1& o3,o2& o4)
{
	return sqrt(pow((o3.get_xc() - o4.get_xc()),2) + pow((o3.get_yc() - o4.get_yc()),2));
}

template < typename o1,typename o2,typename o3,typename o4>
double edist(o1& o5,o2& o6,o3& o7,o4& o8)
{
	return sqrt(pow((o5 - o6),2) + pow((o7 - o8),2));
}

// rounding functions
double roundupx(double num,int x)
{
	return ceil( ( num * pow( 10.0,x ) ) - 0.5 ) / pow( 10.0,x );
}	

double roundownx(double num,int x)
{
	return floor( ( num * pow( 10.0,x ) ) - 0.5 ) / pow( 10.0,x );
}	

template <typename m,typename n,typename h, typename f>
void writeObjData(m& m1, n& n1,h& h1,f& f1)
{
	m::iterator it;
// write out the obejct container m - container, n - object, h - header text, f - output stream
// "Summary of Vertex attributes by id " <<endl
//    << "Id" << "\t"<< "Stop" << "\t" << " id " << "\t" << " cost "<< "\t" << "TotCost "
//	 << "\t"<< "ParCnt"<< "\t"<<"index "<< "\t"<<"LowLink "
//			f1<<h1<< endl;
			for (it=m1.begin();it!=m1.end();++it)
			{
			    f1 << *it;
			}

}

 string datetimeStamp(ofstream& pts)
 {
         char dateStr [9];
         char timeStr [9];
		 string datetime="";
         _strdate( dateStr);
         _strtime( timeStr );
         pts<<(" Date %s  ", dateStr);
         pts<<(" Time %s  ", timeStr)<<endl;
		datetime.append("Date : " );
		datetime.append(dateStr);
		datetime.append(" Time : " );
		datetime.append(timeStr);
		return datetime;
 }

template <typename s, typename t,typename i>
i& getLoc(s& s1, t& t1,i& i1)
{
	typedef i i0;
	long p1;
	p1=s1.find(t1);
	i1 = from_String<i>(s1.substr(0,p1));
	return i1;
}

template <typename s, typename t>
s& getSubstr(s& s1, t& t1)
{
	long p1;
	s& s2=s1;
	p1=s1.find(t1);
	s2 = s1.substr(p1+1,s1.length()-p1-1);
	return s2;
}


#endif // MISCFUNCTIONS_H ///:~

/*
import struct, datetime, decimal, itertools

def dbfreader(f):
    """Returns an iterator over records in a Xbase DBF file.

    The first row returned contains the field names.
    The second row contains field specs: (type, size, decimal places).
    Subsequent rows contain the data records.
    If a record is marked as deleted, it is skipped.

    File should be opened for binary reads.

    """
    # See DBF format spec at:
    #     http://www.pgts.com.au/download/public/xbase.htm#DBF_STRUCT

void dbfread (ifstream dbfile) {
  char vertinfname[ MaxStrLen +1]="C:\\NEU\\bustops\\Boston\\VertexFC0b.dbf\0"; // vertex data input file
  ifstream in(vertinfname, ios::in | ios::binary);
  assure(in, vertinfname);
  in.read(readData[0], STR_LEN);
  assert(strcmp(readData[0], "Hickory dickory dus. . .")

int fd = open ("filename", O_RDONLY);
rec r;
while (read(fd, &r, sizeof(r)))
do_something_with_data(rec);
close(fd);

}
*/
/*
    numrec, lenheader = struct.unpack('<xxxxLH22x', f.read(32))    
    numfields = (lenheader - 33) / 32

    fields = []
    for fieldno in xrange(numfields):
        name, typ, size, deci = struct.unpack('<11sc4xBB14x', f.read(32))
        name = name.replace('\0', '')       # eliminate NULs from string   
        fields.append((name, typ, size, deci))
    yield [field[0] for field in fields]
    yield [tuple(field[1:]) for field in fields]

    terminator = f.read(1)
    assert terminator == '\r'

    fields.insert(0, ('DeletionFlag', 'C', 1, 0))
    fmt = ''.join(['%ds' % fieldinfo[2] for fieldinfo in fields])
    fmtsiz = struct.calcsize(fmt)
    for i in xrange(numrec):
        record = struct.unpack(fmt, f.read(fmtsiz))
        if record[0] != ' ':
            continue                        # deleted record
        result = []
        for (name, typ, size, deci), value in itertools.izip(fields, record):
            if name == 'DeletionFlag':
                continue
            if typ == "N":
                value = value.replace('\0', '').lstrip()
                if value == '':
                    value = 0
                elif deci:
                    value = decimal.Decimal(value)
                else:
                    value = int(value)
            elif typ == 'D':
                y, m, d = int(value[:4]), int(value[4:6]), int(value[6:8])
                value = datetime.date(y, m, d)
            elif typ == 'L':
                value = (value in 'YyTt' and 'T') or (value in 'NnFf' and 'F') or '?'
            result.append(value)
        yield result

*/

/*
void alternative_pattern()
{
// maximum number of stops to be removed l = 4
		int i=0,j=0,k=0,l=4,m=0,n=0,o=0,p=0;

	// pattern making 
		 for (i = 0;i<=l;i++)
		{
			for (j = i;j<=l;j++)
			{
				if (i==j) 
				{ k=i+1;}
				else 
				{ k=i+j+2;} // k - pattern repetition for a specific i,j combination
			// make stop pattern for this set keep1=i,skip2=j starting at i from the start
				for (m=1;m<=k;m++) 
				{
					n = tsidpmap.size();
					o=m;
					while (o<n)
					{
						o=o+i; // i - keep equal number of stops
						for (p=o+1;p<(o+i+j);p++) 
						{ //j - remove equal number of stops
							tsidmit = tsidpmap.find(p);
							tsidpmap.erase(tsidmit);
						}
						o=p;
					}
				}
			}
		 }

}
*/

/*
    DFS(G)
    {
    make a new vertex x with edges x->v for all v
    initialize a counter N to zero
    initialize list L to empty
    build directed tree T, initially a single vertex {x}
    visit(x)
    }

    visit(p)
    {
    add p to L
    dfsnum(p) = N
    increment N
    low(p) = dfsnum(p)
    for each edge p->q
        if q is not already in T
        {
        add p->q to T
        visit(q)
        low(p) = min(low(p), low(q))
        } else low(p) = min(low(p), dfsnum(q))

    if low(p)=dfsnum(p)
    {
        output "component:"
        repeat
        remove last element v from L
        output v
        remove v from G
        until v=p
    }
    }
*/

/* pointers in classes how to copy data into an embeded object
class Dog {
  string nm;
public:
  Dog(const string& name) : nm(name) {
    cout << "Creating Dog: " << *this << endl;
  }
  // Synthesized copy-constructor & operator= 
  // are correct.
  // Create a Dog from a Dog pointer:
  Dog(const Dog* dp, const string& msg) 
    : nm(dp->nm + msg) {
    cout << "Copied dog " << *this << " from "
         << *dp << endl;
  }
  ~Dog() { 
    cout << "Deleting Dog: " << *this << endl;
  }
  void rename(const string& newName) {
    nm = newName;
    cout << "Dog renamed to: " << *this << endl;
  }
  friend ostream&
  operator<<(ostream& os, const Dog& d) {
    return os << "[" << d.nm << "]";
  }
};

class DogHouse {
  Dog* p;
  string houseName;
public:
  DogHouse(Dog* dog, const string& house)
   : p(dog), houseName(house) {}
  DogHouse(const DogHouse& dh)
    : p(new Dog(dh.p, " copy-constructed")),
      houseName(dh.houseName 
        + " copy-constructed") {}
  DogHouse& operator=(const DogHouse& dh) {
    // Check for self-assignment:
    if(&dh != this) {
      p = new Dog(dh.p, " assigned");
      houseName = dh.houseName + " assigned";
    }
    return *this;
  }
  void renameHouse(const string& newName) {
    houseName = newName;
  }
  Dog* getDog() const { return p; }
  ~DogHouse() { delete p; }
  friend ostream&
  operator<<(ostream& os, const DogHouse& dh) {
    return os << "[" << dh.houseName 
      << "] contains " << *dh.p;
  }
}; 
*/

#ifndef SS_Point_H
#define SS_Point_H

class Vector;

//==================================================================
//  Point Class Definition
//==================================================================

class Point {
friend class Vector;
friend class Line;
friend class Linep;
friend class vertexp;
friend class ShapePt;
protected:
	int dimn;            // # coords (1, 2, or 3 max here)
	Error err;           // error indicator
public:
	double x, y, z;      // z=0 for 2D, y=z=0 for 1D
	
	//----------------------------------------------------------
	// Lots of Constructors (add more as needed)
	Point() { dimn=3; x=y=z=0; err=Enot; }
	// 1D Point
	Point( int a) { dimn=1; x=a; y=z=0; err=Enot; }
	Point( double a) { dimn=1; x=a; y=z=0; err=Enot; }
	// 2D Point
	Point( int a, int b) { dimn=2; x=a; y=b; z=0; err=Enot; }
	Point( double a, double b) { dimn=2; x=a; y=b; z=0; err=Enot; }
	// 3D Point
	Point( int a, int b, int c) { dimn=3; x=a; y=b; z=c; err=Enot; }
	Point( double a, double b, double c) { dimn=3; x=a; y=b; z=c; err=Enot; }
	// n-dim Point
//	Point( int n, int a[]);
//	Point( int n, double a[]);
	// Destructor
	~Point() {};

	int dim() { return dimn; }      // get dimension
//	int setdim( int);               // set new dimension

	Vector operator-( Point);       // Vector difference


	//----------------------------------------------------------
	// Point Scalar Operations (convenient but often illegal)
	// using any type of scalar (int, float, or double)
	//    are not valid for points in general,
	//    unless they are 'affine' as coeffs of 
	//    a sum in which all the coeffs add to 1,
	//    such as: the sum (a*P + b*Q) with (a+b == 1).
	//    The programmer must enforce this (if they want to).

	// Scalar Multiplication
	friend Point operator*( int, Point);
	friend Point operator*( double, Point);
	friend Point operator*( Point, int);
	friend Point operator*( Point, double);
	// Scalar Division
	friend Point operator/( Point, int);
	friend Point operator/( Point, double);
	friend Point operator+( Point, Point);     // add points

	//----------------------------------------------------------
	// Point and Vector Operations (always valid) 
	Point  operator+( Vector);      // +translate
	Point  operator-( Vector);      // -translate
	Point& operator+=( Vector);     // inc translate
	Point& operator-=( Vector);     // dec translate


/*
	//----------------------------------------------------------
	// Input/Output streams
	friend istream& operator>>( istream&, Point&);
	friend ostream& operator<<( ostream&, Point);

	//----------------------------------------------------------
	// Assignment "=": use the default to copy all members

	//----------------------------------------------------------
	// Comparison (dimension must match, or not)
	int operator==( Point);
	int operator!=( Point);

	//----------------------------------------------------------
	// Point Addition (also convenient but often illegal)
	//    is not valid unless part of an affine sum.
	//    The programmer must enforce this (if they want to).

	// Affine Sum
	// Returns weighted sum, even when not affine, but...
	// Tests if coeffs add to 1.  If not, sets: err = Esum.
	friend Point asum( int, int[], Point[]);
	friend Point asum( int, double[], Point[]);

	//----------------------------------------------------------
	// Point Relations
	friend double d( Point, Point);         // Distance
	friend double d2( Point, Point);        // Distance^2
	double isLeft( Point, Point);           // 2D only
	double Area( Point, Point); 		// any dim for triangle PPP
	void  set_x( double x1); 		// x coordinate of point
	void  set_y( double y1); 		// y coordinate of point
	void  set_z( double z1); 		// z coordinate of point
	double get_x(); 		// x coordinate of point
	double get_y(); 		// y coordinate of point
	double get_z(); 		// z coordinate of point

	// Collinearity Conditions (any dim n)
	boolean isOnLine( Point, Point, char);  // is On line (char= flag)
	boolean isOnLine( Point, Point);        // is On line (flag= all)
	boolean isBefore( Point, Point);        // is On line (flag= before)
	boolean isBetween( Point, Point);       // is On line (flag= between)
	boolean isAfter( Point, Point);         // is On line (flag= after)
	boolean isOnRay( Point, Point);         // is On line (flag= between|after)

	//----------------------------------------------------------
	// Error Handling
	void  clerr() { err = Enot;}            // clear error
	int   geterr() { return err;}           // get error
	char* errstr();                         // return error string
	void Point::serialize(ofstream& pVx); // write binary output
	void Point::deserialize(ifstream& pVx); // read binary input
	void Point::deserializetext(ifstream& pVx); // read binary input

*/
	
//==================================================================
// Point Class Methods
//==================================================================

//------------------------------------------------------------------
// Constructors (add more as needed)
//------------------------------------------------------------------

// n-dim Point
Point::Point( int n, int a[]) {
	x = y = z = 0;
	err = Enot;
	switch (dimn = n) {
	case 3: z = a[2];
	case 2: y = a[1];
	case 1: x = a[0];
		break;
	default:
		err=Edim;
	}
}

Point::Point( int n, double a[]) {
	x = y = z = 0.0;
	err = Enot;
	switch (dimn = n) {
	case 3: z = a[2];
	case 2: y = a[1];
	case 1: x = a[0];
		break;
	default:
		err=Edim;
	}
};


//------------------------------------------------------------------
// Assign (set) dimension
//------------------------------------------------------------------

int Point::setdim( int n) {
	switch (n) {
	case 1: y = 0;
	case 2: z = 0;
	case 3:
		return dimn = n;
	default:                      // out of range value
		err = Edim;           // just flag the error
		return ERROR;
	}
}

//------------------------------------------------------------------
// Comparison (note: dimension must compare)
//------------------------------------------------------------------

int Point::operator==( Point Q)
{
	if (dimn != Q.dim()) return FALSE;
	switch (dimn) {
	case 1:
		return (x==Q.x);
	case 2:
		return (x==Q.x && y==Q.y);
	case 3:
	default:
		return (x==Q.x && y==Q.y && z==Q.z);
	}
}

int Point::operator!=( Point Q)
{
	if (dimn != Q.dim()) return TRUE;
	switch (dimn) {
	case 1:
		return (x!=Q.x);
	case 2:
		return (x!=Q.x || y!=Q.y);
	case 3:
	default:
		return (x!=Q.x || y!=Q.y || z!=Q.z);
	}
}

//------------------------------------------------------------------
// Affine Sums
// Returns weighted sum, even when not affine, but...
// Tests if coeffs add to 1.  If not, sets: err = Esum.
//------------------------------------------------------------------

Point asum( int n, int c[], Point Q[])
{
	int        maxd = 0;
	int        cs = 0;
	Point      P;

	for (int i=0; i<n; i++) {
		cs += c[i];
		if (Q[i].dim() > maxd)
			maxd = Q[i].dim();
	}
	if (cs != 1)                 // not an affine sum
		P.err = Esum;        // flag error, but compute sum anyway

	for (int i=0; i<n; i++) {
		P.x += c[i] * Q[i].x;
		P.y += c[i] * Q[i].y;
		P.z += c[i] * Q[i].z;
	}
	P.dimn = maxd;
	return P;
}

Point asum( int n, double c[], Point Q[])
{
	int        maxd = 0;
	double     cs = 0.0;
	Point      P;

	for (int i=0; i<n; i++) {
		cs += c[i];
		if (Q[i].dim() > maxd)
			maxd = Q[i].dim();
	}
	if (cs != 1)                 // not an affine sum
		P.err = Esum;        // flag error, but compute sum anyway

	for (int i=0; i<n; i++) {
		P.x += c[i] * Q[i].x;
		P.y += c[i] * Q[i].y;
		P.z += c[i] * Q[i].z;
	}
	P.dimn = maxd;
	return P;
}

//------------------------------------------------------------------
// Distance between Points
//------------------------------------------------------------------

double ds( Point P, Point Q) {      // Euclidean distance
	double dx = P.x - Q.x;
	double dy = P.y - Q.y;
	double dz = P.z - Q.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
};

double d2( Point P, Point Q) {     // squared distance (more efficient)
	double dx = P.x - Q.x;
	double dy = P.y - Q.y;
	double dz = P.z - Q.z;
	return (dx*dx + dy*dy + dz*dz);
};

//------------------------------------------------------------------
// Sidedness of a Point wrt a directed line P1->P2
//        - makes sense in 2D only
//------------------------------------------------------------------

double Point::isLeft( Point P1, Point P2) {
	if (dimn != 2 || P1.dim() != 2 || P2.dim() != 2) {
		err = Edim;        // flag error, but compute anyway
	}
	return ((P1.x - x) * (P2.y - y) - (P2.x - x) * (P1.y - y));
};

//------------------------------------------------------------------
// Error Routines
//------------------------------------------------------------------

char* Point::errstr() {            // return error string
	switch (err) {
	case Enot:
		return "no error";
	case Edim:
		return "error: invalid dimension for operation";
	case Esum:
		return "error: Point sum is not affine";
	default:
		return "error: unknown err value";
	}
};

//------------------------------------------------------------------
// Binary output/input Routines
//------------------------------------------------------------------

void Point::serialize(ofstream& pVx)
{
 pVx.write(reinterpret_cast<char *>(&x), sizeof(x));
 pVx.write(reinterpret_cast<char *>(&y), sizeof(y));
 pVx.write(reinterpret_cast<char *>(&z), sizeof(z));
 pVx.write(reinterpret_cast<char *>(&dimn), sizeof(dimn));
 pVx.write(reinterpret_cast<char *>(&err), sizeof(err));
}

void Point::deserialize(ifstream& pVx)
{
 pVx.read(reinterpret_cast<char *>(&x), sizeof(x));
 pVx.read(reinterpret_cast<char *>(&y), sizeof(y));
 pVx.read(reinterpret_cast<char *>(&z), sizeof(z));
 pVx.read(reinterpret_cast<char *>(&dimn), sizeof(dimn));
 pVx.read(reinterpret_cast<char *>(&err), sizeof(err));
}

void Point::set_x(double x1)
{
	x = x1;
}

void Point::set_y(double y1)
{
	y = y1;
}

void Point::set_z(double z1)
{
	z = z1;
}

double Point::get_x()
{
	return x;
}

double Point::get_y()
{
	return y;
}

double Point::get_z()
{
	return z;
}


};


//------------------------------------------------------------------
// IO streams
//------------------------------------------------------------------

// Read input Point format: "(%f)", "(%f, %f)", or "(%f, %f, %f)"
istream& operator>>( istream& input, Point& P) {
	char c;
	input >> c;                // skip '('
	input >> P.x;
	input >> c;                
	if (c == ')') {
		P.setdim(1);       // 1D coord
		return input;
	}
	// else                    // skip ','
	input >> P.y;
	input >> c;
	if (c == ')') {
		P.setdim(2);       // 2D coord
		return input;
	}
	// else                    // skip ','
	input >> P.z;
	P.setdim(3);               // 3D coord
	input >> c;                // skip ')'
	return input;
}

// Write output Point in format: "(%f)", "(%f, %f)", or "(%f, %f, %f)"
ostream& operator<<( ostream& output, Point P) {
	switch (P.dim()) {
	case 1:
		output << "(" << P.x << ")";
		break;
	case 2:
		output << "(" << P.x << ", " << P.y << ")";
		break;
	case 3:
		output << "(" << P.x << ", " << P.y << ", " << P.z << ")";
		break;
	default:
		output << "Error: P.dim = " << P.dim();
	}
	return output;
}

//------------------------------------------------------------------
// Point Scalar Operations (convenient but often illegal)
//        are not valid for points in general,
//        unless they are 'affine' as coeffs of 
//        a sum in which all the coeffs add to 1,
//        such as: the sum (a*P + b*Q) with (a+b == 1).
//        The programmer must enforce this (if they want to).
//------------------------------------------------------------------

Point operator*( int c, Point Q ) {
	Point P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	P.dimn = Q.dim();
	return P;
}

Point operator*( double c, Point Q) {
	Point P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	P.dimn = Q.dim();
	return P;
}

Point operator*( Point Q, int c) {
	Point P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	P.dimn = Q.dim();
	return P;
}

Point operator*( Point Q, double c) {
	Point P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	P.dimn = Q.dim();
	return P;
}

Point operator/( Point Q, int c) {
	Point P;
	P.x = Q.x / c;
	P.y = Q.y / c;
	P.z = Q.z / c;
	P.dimn = Q.dim();
	return P;
}

Point operator/( Point Q, double c) {
	Point P;
	P.x = Q.x / c;
	P.y = Q.y / c;
	P.z = Q.z / c;
	P.dimn = Q.dim();
	return P;
}

//------------------------------------------------------------------
// Point Addition (also convenient but often illegal)
//    is not valid unless part of an affine sum.
//    The programmer must enforce this (if they want to).
//------------------------------------------------------------------

Point operator+( Point Q, Point R)
{
	Point P;
	P.x = Q.x + R.x;
	P.y = Q.y + R.y;
	P.z = Q.z + R.z;
	P.dimn = max( Q.dim(), R.dim());
	return P;
}



#endif SS_Point_H

#ifndef SS_Vector_H
#define SS_Vector_H

//==================================================================
//  Vector Class Definition
//==================================================================

class Vector : public Point {
public:
	// Constructors same as Point class
	Vector() : Point() {};
	Vector( int a) : Point(a) {};
	Vector( double a) : Point(a) {};
	Vector( int a, int b) : Point(a,b) {};
	Vector( double a, double b) : Point(a,b) {};
	Vector( int a, int b, int c) : Point(a,b,c) {};
	Vector( double a, double b, double c) : Point(a,b,c) {};
	Vector( int n, int a[]) : Point(n,a) {};
	Vector( int n, double a[]) : Point(n,a) {};
	~Vector() {};
	//----------------------------------------------------------
	// IO streams and Comparisons: inherit from Point class

	//----------------------------------------------------------
	// Vector Unary Operations
	Vector operator-();                // unary minus
	Vector operator~();                // unary 2D perp operator

	//----------------------------------------------------------
	// Scalar Multiplication
	friend Vector operator*( int, Vector);
	friend Vector operator*( double, Vector);
	friend Vector operator*( Vector, int);
	friend Vector operator*( Vector, double);
	// Scalar Division
	friend Vector operator/( Vector, int);
	friend Vector operator/( Vector, double);

	//----------------------------------------------------------
	// Vector Arithmetic Operations
	Vector operator+( Vector);        // vector add
	Vector operator-( Vector);        // vector subtract
	double operator*( Vector);        // inner dot product
	double operator|( Vector);        // 2D exterior perp product
	Vector operator^( Vector);        // 3D exterior cross product

	Vector& operator*=( double);      // vector scalar mult
	Vector& operator/=( double);      // vector scalar div
	Vector& operator+=( Vector);      // vector increment
	Vector& operator-=( Vector);      // vector decrement
	Vector& operator^=( Vector);      // 3D exterior cross product
	//----------------------------------------------------------
	// Vector Properties
	double len() {                    // vector length
		return sqrt(x*x + y*y + z*z);
	}
	double len2() {                   // vector length squared (faster)
		return (x*x + y*y + z*z);
	}
/*
	//----------------------------------------------------------
	// Special Operations
	void normalize();                 // convert vector to unit length
	friend Vector sum( int, int[], Vector[]);     // vector sum
	friend Vector sum( int, double[], Vector[]);  // vector sum
*/

//------------------------------------------------------------------
//  Special Operations
//------------------------------------------------------------------

void Vector::normalize() {                      // convert to unit length
	double ln = sqrt( x*x + y*y + z*z );
	if (ln == 0) return;                    // do nothing for nothing
	x /= ln;
	y /= ln;
	z /= ln;
}

Vector sum( int n, int c[], Vector w[] ) {     // vector sum
	int     maxd = 0;
	Vector  v;

	for (int i=0; i<n; i++) {
		if (w[i].dim() > maxd)
			maxd = w[i].dim();
	}
	v.dimn = maxd;

	for (int i=0; i<n; i++) {
		v.x += c[i] * w[i].x;
		v.y += c[i] * w[i].y;
		v.z += c[i] * w[i].z;
	}
	return v;
}

Vector sum( int n, double c[], Vector w[] ) {  // vector sum
	int     maxd = 0;
	Vector  v;

	for (int i=0; i<n; i++) {
		if (w[i].dim() > maxd)
			maxd = w[i].dim();
	}
	v.dimn = maxd;

	for (int i=0; i<n; i++) {
		v.x += c[i] * w[i].x;
		v.y += c[i] * w[i].y;
		v.z += c[i] * w[i].z;
	}
	return v;
}

};

//------------------------------------------------------------------
// Point Vector Operations
//------------------------------------------------------------------

Vector Point::operator-( Point Q)        // Vector diff of Points
{
	Vector v;
	v.x = x - Q.x;
	v.y = y - Q.y;
	v.z = z - Q.z;
	v.dimn = max( dimn, Q.dim());
	return v;
}

Point Point::operator+( Vector v)        // +ve translation
{
	Point P;
	P.x = x + v.x;
	P.y = y + v.y;
	P.z = z + v.z;
	P.dimn = max( dimn, v.dim());
	return P;
}

Point Point::operator-( Vector v)        // -ve translation
{
	Point P;
	P.x = x - v.x;
	P.y = y - v.y;
	P.z = z - v.z;
	P.dimn = max( dimn, v.dim());
	return P;
}

Point& Point::operator+=( Vector v)        // +ve translation
{
	x += v.x;
	y += v.y;
	z += v.z;
	dimn = max( dimn, v.dim());
	return *this;
}

Point& Point::operator-=( Vector v)        // -ve translation
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	dimn = max( dimn, v.dim());
	return *this;
}


//==================================================================
// Vector Class Methods
//==================================================================

//------------------------------------------------------------------
//  Unary Ops
//------------------------------------------------------------------

// Unary minus
Vector Vector::operator-() {
	Vector v;
	v.x = -x; v.y = -y; v.z = -z;
	v.dimn = dimn;
	return v;
};

// Unary 2D perp operator
Vector Vector::operator~() {
	if (dimn != 2) err = Edim;   // and continue anyway
	Vector v;
	v.x = -y; v.y = x; v.z = z;
	v.dimn = dimn;
	return v;
};

//------------------------------------------------------------------
//  Scalar Ops
//------------------------------------------------------------------

// Scalar multiplication
Vector operator*( int c, Vector w ) {
	Vector v;
	v.x = c * w.x;
	v.y = c * w.y;
	v.z = c * w.z;
	v.dimn = w.dim();
	return v;
};

Vector operator*( double c, Vector w ) {
	Vector v;
	v.x = c * w.x;
	v.y = c * w.y;
	v.z = c * w.z;
	v.dimn = w.dim();
	return v;
};

Vector operator*( Vector w, int c ) {
	Vector v;
	v.x = c * w.x;
	v.y = c * w.y;
	v.z = c * w.z;
	v.dimn = w.dim();
	return v;
};

Vector operator*( Vector w, double c ) {
	Vector v;
	v.x = c * w.x;
	v.y = c * w.y;
	v.z = c * w.z;
	v.dimn = w.dim();
	return v;
};

// Scalar division
Vector operator/( Vector w, int c ) {
	Vector v;
	v.x = w.x / c;
	v.y = w.y / c;
	v.z = w.z / c;
	v.dimn = w.dim();
	return v;
}

Vector operator/( Vector w, double c ) {
	Vector v;
	v.x = w.x / c;
	v.y = w.y / c;
	v.z = w.z / c;
	v.dimn = w.dim();
	return v;
}

//------------------------------------------------------------------
//  Arithmetic Ops
//------------------------------------------------------------------

Vector Vector::operator+( Vector w ) {
	Vector v;
	v.x = x + w.x;
	v.y = y + w.y;
	v.z = z + w.z;
	v.dimn = max( dimn, w.dim());
	return v;
}

Vector Vector::operator-( Vector w ) {
	Vector v;
	v.x = x - w.x;
	v.y = y - w.y;
	v.z = z - w.z;
	v.dimn = max( dimn, w.dim());
	return v;
}

//------------------------------------------------------------------
//  Products
//------------------------------------------------------------------

// Inner Dot Product
double Vector::operator*( Vector w ) {
	return (x * w.x + y * w.y + z * w.z);
}

// 2D Exterior Perp Product
double Vector::operator|( Vector w ) {
	if (dimn != 2) err = Edim;    // and continue anyway
	return (x * w.y - y * w.x);
}

// 3D Exterior Cross Product
Vector Vector::operator^( Vector w ) {
	Vector v;
	v.x = y * w.z - z * w.y;
	v.y = z * w.x - x * w.z;
	v.z = x * w.y - y * w.x;
	v.dimn = 3;
	return v;
}

//------------------------------------------------------------------
//  Shorthand Ops
//------------------------------------------------------------------

Vector& Vector::operator*=( double c ) {        // vector scalar mult
	x *= c;
	y *= c;
	z *= c;
	return *this;
}

Vector& Vector::operator/=( double c ) {        // vector scalar div
	x /= c;
	y /= c;
	z /= c;
	return *this;
}

Vector& Vector::operator+=( Vector w ) {        // vector increment
	x += w.x;
	y += w.y;
	z += w.z;
	dimn = max(dimn, w.dim());
	return *this;
}

Vector& Vector::operator-=( Vector w ) {        // vector decrement
	x -= w.x;
	y -= w.y;
	z -= w.z;
	dimn = max(dimn, w.dim());
	return *this;
}

Vector& Vector::operator^=( Vector w ) {        // 3D exterior cross product
	double ox=x, oy=y, oz=z;
	x = oy * w.z - oz * w.y;
	y = oz * w.x - ox * w.z;
	z = ox * w.y - oy * w.x;
	dimn = 3;
	return *this;
}




#endif SS_Vector_H

#ifndef SS_Line_H
#define SS_Line_H

class Linep
{
    private:
        Point *P0;
        Point *P1;
    public:
        Linep( Point &Pb,  Point &Pe) : P0(&Pb), P1(&Pe) {}

        void setPoints(Point &Pb, Point &Pe)
        {
            this->P0 = &Pb;
            this->P1 = &Pe;
        }
		double getX ()
		{
			return this->P0->x;
		}
		double getY ()
		{
			return this->P0->y;
		}
};

class Line
{
    private:
        Point P0;
        Point P1;
    public:
        Line(const Point & pb, const Point & pe) : P0(pb), P1(pe) {}

        void setPoints(const Point & pb, const Point & pe)
        {
            this->P0 = pb;
            this->P1 = pe;
        }
		const Point &startPoint() const
        {
            return P0;
        }

        Point &startPoint()
        {
            return P0;
        }

        const Point &endPoint() const
        {
            return P1;
        }

        Point &endPoint()
        {
            return P1;
        }

		double get_length() 
		{
			double dx = (P1.get_x()-P0.get_x());
			double dy = (P1.get_y()-P0.get_y());
			return (sqrt(dy*dy + dx*dx));
		}
		};


// closest2D_Point_to_Line(): finds the closest 2D Point to a Line
//    Input:  an array P[] of n points, and a Line L
//    Return: the index i of the Point P[i] closest to L

template <typename m, typename o, typename p, typename  k, typename  x >
k& closestdistpoint(m& m0, o& o0, p& p0, k& k0,x& x0)
{
	o* o1;
	o1=&o0;
	p p1;
	p p2;
	p* pp1=&p1;
	p* pp2=&p2;
    typedef pair <k,o> kpair;
    // initialize min index and distance to first point
	double x,y,min=0,dist=0;
	m::iterator mit;
	x::iterator xit;
	k k1;
		min = (1e+10);

// use the ordered collection of points using int(x*10) to get the minimum distance
    // loop through Point array testing for min distance to L
	for(mit=m0.begin();mit!=m0.end();mit++)
	{
		o0 = mit->second;
		// get the point for this shape
		p1 = o0.get_pt();
		// calculate the dist
		dist = ds(p1,p0);
	    if (dist < 0) dist = -dist;    // absolute value
		    if (dist < min) {    // this point is closer
			    k1 = o0.get_edge();          // so have a new minimum
				min = dist;
			}
	}

    return k1;    // the index of the closest Point P[mi]
}
//===================================================================

// dist_Point_to_Line(): get the distance of a point to a line.
//    Input:  a Point P and a Line L (in any dimension)
//    Return: the shortest distance from P to L

float dist_Point_to_Line( Point P, Line L)
{
    Vector v = L.endPoint() - L.startPoint();
    Vector w = P - L.startPoint();

    double c1 = dot(w,v);
    double c2 = dot(v,v);
    double b = c1 / c2;

    Point Pb = L.startPoint() + b * v;
	c1 = d(P,Pb);
    return c1;
};
//===================================================================

// dist_Point_to_Segment(): get the distance of a point to a segment.
//    Input:  a Point P and a Segment S (in any dimension)
//    Return: the shortest distance from P to S
float dist_Point_to_Segment( Point P, Line S)
{
    Vector v = S.endPoint() - S.startPoint();
    Vector w = P - S.endPoint();

    double c1 = dot(w,v);
    if ( c1 <= 0 )
        return d(P, S.startPoint());

    double c2 = dot(v,v);
    if ( c2 <= c1 )
        return d(P, S.endPoint());

    double b = c1 / c2;
    Point Pb = S.startPoint() + b * v;
    return d(P, Pb);
};

// closest2D_Point_to_Point(): finds the closest 2D Point to a Point
//    Input:  an array P[] of n points, and another 
//    Return: the index i of the collection to P

int closest2D_Point_to_Line ( Point P[], int n, Line L)
{
    // Get coefficients of the implicit line equation.
    // Do NOT normalize since scaling by a constant
    // is irrelevant for just comparing distances.
    float a = L.startPoint().y - L.endPoint().y;
    float b = L.endPoint().x - L.startPoint().x;
    float c = L.startPoint().x * L.endPoint().y - L.endPoint().x * L.startPoint().y;

    // initialize min index and distance to P[0]
    int mi = 0, i=0;
    float min = a * P[0].x + b * P[0].y + c;
    if (min < 0) min = -min;    // absolute value

    // loop through Point array testing for min distance to L
    for (i=1; i<n; i++) {
        // just use dist squared (sqrt not needed for comparison)
        float dist = a * P[i].x + b * P[i].y + c;
        if (dist < 0) dist = -dist;   // absolute value
        if (dist < min) {    // this point is closer
            mi = i;          // so have a new minimum
            min = dist;
        }
    }
    return mi;    // the index of the closest Point P[mi]
}


#endif SS_Line_H

#if !defined( _VERTEXP_H_ )
#define _VERTEXP_H_ 

// define the vertex object 
class Ptxy;
class Point;

class vertexp  {

	friend class ShapePt;
	friend class Point;

public:
//	vertexp();
	vertexp(const vertexp& vx )  {
		pid = vx.pid;
		pidp = vx.pidp;
		lbl = vx.lbl;
		cost= vx.cost;
		fid = vx.fid;
		orig = vx.orig;
		tcost= vx.tcost;
		toStop= vx.toStop;
		index= vx.index;
		lowlink= vx.lowlink;
		x = vx.x;
		y = vx.y;
		pt = vx.pt;
	//  cout << "Vx[" << pid << "]" << endl;
	//    ++copycons;
  }
  vertexp& operator=(const vertexp& vx) {
    //cout << "(" << pid << ")=[" << vx.pid << "]" << endl;
		pid = vx.pid;
		pidp = vx.pidp;
		lbl = vx.lbl;
		cost= vx.cost;
		fid = vx.fid;
		orig = vx.orig;
		tcost= vx.tcost;
		toStop= vx.toStop;
		index= vx.index;
		lowlink= vx.lowlink;
		x = vx.x;
		y = vx.y;
		pt = vx.pt;
    return *this;
  }

	bool operator==(vertexp vx);

	vertexp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1,short toStop1, long index1,long lowlink1);

	vertexp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1,short toStop1, long index1,long lowlink1, Point pt,double x=-1, double y=-1);

	vertexp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1, double tcost1,short toStop1, long index1,long lowlink1);

	vertexp(Point pt, long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1, double tcost1,short toStop1, long index1,long lowlink1);

	vertexp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1, double tcost1,short toStop1, long index1,long lowlink1, Point pt);

	void vxp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1,short toStop1, long index1,long lowlink1);

	void vxp(Point pt1, long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1,short toStop1, long index1,long lowlink1);

	void set_id (long pid1) {pid=pid1;}

	long get_id () const {return pid;}

	void set_idp (long pidp1) {pidp=pidp1;}

	long get_idp () const {return pidp;}  // predecessor

	void set_lbl (short lbl1) {lbl=lbl1;}

	short get_lbl () const {return lbl;}

	void set_toStop (short toStop1) {toStop=toStop1;}

	short get_toStop () const {return toStop;}

	void set_x (double x1) {x = x1;}

	double get_x () const {return x;}

	void set_y (double y1) {y=y1;}

	double get_y () const {return y;}

	void set_z (double z1) {z=z1;}

	double get_z () const {return z;}

	void set_cost (double cost1) {cost=cost1;}

	double get_cost () const {return cost;}

	void set_tcost (double tcost1) {tcost=tcost1;}

	double get_tcost () const {return tcost;}

	void set_fid (long fid1) {fid=fid1;}

	long get_fid () const {return fid;}

	void set_orig (long orig1) {orig=orig1;}

	long get_orig () const {return orig;}

	void set_lowlink (long lowlink1) {lowlink=lowlink1;}

	long get_lowlink () const {return lowlink;}

	void set_index (long index1) {index=index1;}

	long get_index () const {return index;}

	void set_pt (Point pt1) {pt=pt1;}

	const Point get_pt () const {return pt;}


	void vertexp::serialize(ofstream& pfVx);
	void vertexp::deserialize(ifstream& pfVx);
	void vertexp::serializetext(ofstream& pVx);
	void vertexp::serializetexthdr(ofstream& pVx);

    void show_vertof(ostream& out)
    { 
		out << get_id() << "\t" << get_idp() << "\t" << get_orig() << "\t" 
		  << get_cost() << "\t" << get_tcost() << "\t" << get_lbl()<< "\t" 
		  <<get_fid() << "\t" <<pt<< endl;
	};

    void show_verthdr(ostream& out)
    { 
		out  << "vid" << "\t" << "pid" << "\t" << "orig" << "\t" 
		 << "vcost" << "\t" << "tcost" << "\t" << "lbl"<<"\t"<<"fid" << endl;
	};

	void show_vertex(void)
    { 
		cout << "Vertex Id : " << get_id() << endl;
		cout << "Pred Id: " << get_idp() << endl;
		cout << "Vertex Feature Id: " << get_fid() << endl;
		cout << "Vertex Origin : " << get_orig() << endl;
		cout << "Vertex Cost : " << get_cost() << endl;
		cout << "Vertex Total Cost : " << get_tcost() << endl;
		cout << "Vertex lbl : " << get_lbl() << endl;
	};

void set_origin_vals(long pidp1, short lbl1, long orig1, double cost1, double tcost1,short toStop1=0)
    { 
		vertexp::pidp = pidp1;
		vertexp::lbl = lbl1;
		vertexp::orig=orig1;
		vertexp::cost=cost1;
		vertexp::tcost=tcost1;
	};
	
	void vxp(long pid1=0, long pidp1=0, short lbl1=0, double cost1=inf,  long fid1=0, long orig1=0,
		double tcost1=inf,short toStop1=0, long index1=inf, long lowlink1=inf) 
	{
		vertexp::pid = pid1;
		vertexp::pidp = pidp1;
		vertexp::lbl = lbl1;
		vertexp::cost= cost1;
		vertexp::fid = fid1;
		vertexp::orig = orig1;
		vertexp::tcost= tcost1;
		vertexp::toStop= toStop1;
		vertexp::index= index1;
		vertexp::lowlink= lowlink1;
	}

	void vxp(Point pt1 , long pid1=0, long pidp1=0, short lbl1=0, double cost1=inf,  long fid1=0, long orig1=0,
		double tcost1=inf,short toStop1=0, long index1=inf, long lowlink1=inf) 
	{
		vertexp::pid = pid1;
		vertexp::pidp = pidp1;
		vertexp::lbl = lbl1;
		vertexp::cost= cost1;
		vertexp::fid = fid1;
		vertexp::orig = orig1;
		vertexp::tcost= tcost1;
		vertexp::toStop= toStop1;
		vertexp::index= index1;
		vertexp::lowlink= lowlink1;
		vertexp::pt= pt1;
	}

	virtual ~vertexp(){
		//	cout << "Vertex Object is deleted! "<<endl;
			//delete[] pv;
//		    plist.clear();
	}

protected:
	unsigned long pid;  // vertex id
	         long pidp; // predecessor
	         short lbl; // label
         	 double cost;  // cost at vertex
         	 double tcost; // total cost
			 long fid;  //  
			 long orig;  // Origin for vertex least cost path
			 short toStop;  // travel direction toStop=1 (on) or from stop=0 (off) 
			 long index;  // index of vertex in the graph of the Strongly Connected Component (SCC)
			 long lowlink;  // indicator of the root index of the SCC, if index=lowlink, this is root
			 double x; // x  data 
			 double y; // y  data 
			 double z; // z  data 
			 Point pt; // x,y,z point data 
};

void vertexp::serialize(ofstream& pVx)
{
 pVx.write(reinterpret_cast<char *>(&pid), sizeof(pid));
 pVx.write(reinterpret_cast<char *>(&pidp), sizeof(pidp));
 pVx.write(reinterpret_cast<char *>(&lbl), sizeof(lbl));
 pVx.write(reinterpret_cast<char *>(&cost), sizeof(cost));
 pVx.write(reinterpret_cast<char *>(&tcost), sizeof(tcost));
 pVx.write(reinterpret_cast<char *>(&fid), sizeof(fid));
 pVx.write(reinterpret_cast<char *>(&orig),sizeof(orig));
 pVx.write(reinterpret_cast<char *>(&toStop), sizeof(toStop));
 pVx.write(reinterpret_cast<char *>(&index), sizeof(index));
 pVx.write(reinterpret_cast<char *>(&lowlink), sizeof(lowlink));
 pVx.write(reinterpret_cast<char *>(&x), sizeof(x));
 pVx.write(reinterpret_cast<char *>(&y), sizeof(y));
 pVx.write(reinterpret_cast<char *>(&z), sizeof(z));
 pt.serialize(pVx);
}

void vertexp::deserialize(ifstream& pVx)
{
 pVx.read(reinterpret_cast<char *>(&pid), sizeof(pid));
 pVx.read(reinterpret_cast<char *>(&pidp), sizeof(pidp));
 pVx.read(reinterpret_cast<char *>(&lbl), sizeof(lbl));
 pVx.read(reinterpret_cast<char *>(&cost), sizeof(cost));
 pVx.read(reinterpret_cast<char *>(&tcost), sizeof(tcost));
 pVx.read(reinterpret_cast<char *>(&fid), sizeof(fid));
 pVx.read(reinterpret_cast<char *>(&orig),sizeof(orig));
 pVx.read(reinterpret_cast<char *>(&toStop), sizeof(toStop));
 pVx.read(reinterpret_cast<char *>(&index), sizeof(index));
 pVx.read(reinterpret_cast<char *>(&lowlink), sizeof(lowlink));
 pVx.read(reinterpret_cast<char *>(&x), sizeof(x));
 pVx.read(reinterpret_cast<char *>(&y), sizeof(y));
 pVx.read(reinterpret_cast<char *>(&z), sizeof(z));
 pt.deserialize(pVx);
}

vertexp::vertexp(long pid1=0, long pidp1=0, short lbl1=0, double cost1=inf,  long fid1=0, long orig1=0,
				 double tcost1=inf,short toStop1=0, long index1=inf, long lowlink1=inf) 
{
	vertexp::pid = pid1;
	vertexp::pidp = pidp1;
	vertexp::lbl = lbl1;
	vertexp::cost= cost1;
	vertexp::fid = fid1;
	vertexp::orig = orig1;
	vertexp::tcost= tcost1;
	vertexp::toStop= toStop1;
	vertexp::index= index1;
	vertexp::lowlink= lowlink1;
}

vertexp::vertexp(Point pt1, long pid1=0, long pidp1=0, short lbl1=0, double cost1=inf,  long fid1=0, long orig1=0,
				 double tcost1=inf,short toStop1=0, long index1=inf, long lowlink1=inf) 
{
	vertexp::pid = pid1;
	vertexp::pidp = pidp1;
	vertexp::lbl = lbl1;
	vertexp::cost= cost1;
	vertexp::fid = fid1;
	vertexp::orig = orig1;
	vertexp::tcost= tcost1;
	vertexp::toStop= toStop1;
	vertexp::index= index1;
	vertexp::lowlink= lowlink1;
	vertexp::pt = pt1;
}

void vertexp::serializetext(ofstream& pVx)
{
 pVx <<pid<<"\t"<<pidp<<"\t"<<lbl<<"\t"<<cost<<"\t"<<tcost<<"\t"<<fid<<"\t"<<orig<<
	 "\t"<<toStop<<"\t"<<index<<"\t"<<lowlink<<"\t"<<x<<"\t"<<"\t"<<y<<"\t"<<z<<pt<<endl;
}
void vertexp::serializetexthdr(ofstream& pVx)
{
 pVx <<"pid"<<"\t"<<"pidp"<<"\t"<<"lbl"<<"\t"<<"cost"<<"\t"<<"tcost"<<"\t"<<"fid"<<"\t"
	 <<"orig"<<"\t"<<"toStop"<<"\t"<<"index"<<"\t"<<"lowlink"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"z"<<"\t"<<"Point details"<<endl;
}

bool vertexp::operator==(vertexp vx)
 {
   if(index!=vx.index)
      return false;
   if(lowlink!=vx.lowlink)
      return false;
   if(orig!=vx.orig)
      return false;
   if(pid!=vx.pid)
      return false;
   if(pidp!=vx.pidp)
      return false;
   if(lbl!=vx.lbl)
      return false;
   if(cost!=vx.cost)
      return false;
   if(toStop!=vx.toStop)
      return false;
   if(tcost!=vx.tcost)
      return false;
   if(fid!=vx.fid)
      return false;
   return true;
 }


#endif // _VERTEXP_H_ 

#ifndef ShapePt_H //Walk, Ride, Operating cost
#define ShapePt_H

class ShapePt : public vertexp {

	friend class vertexp;
	friend class Point;

public:
 
	ShapePt() : vertexp() {};

	void set_edge (long edge1) {edge=edge1;}

	long get_edge () {return edge;}

	~ShapePt() {};

private:
    long edge;    

};

#endif // ShapePt_H_ 



class VxVorAccum {
	long id;
	unsigned long pid;  // vertex id
	         long pidp; // predecessor
	         short lbl; // label
         	 double cost;  // cost at vertex
         	 double tcost; // total cost
			 long fid;  //  
			 long orig;  // Origin for vertex least cost path
			 short toStop;  // travel direction toStop=1 (on) or from stop=0 (off) 
			 long index;  // index of vertex in the graph of the Strongly Connected Component (SCC)
			 long lowlink;  // indicator of the root index of the SCC, if index=lowlink, this is root
			long Vxcnt; // count 
public:
  VxVorAccum() : id(0), pid(0),pidp(0), lbl(0), cost(0),tcost(0),fid(0),orig(0), 
	  toStop(0),index(0),lowlink(inf),Vxcnt(0) {}

  void operator()(const vertexp& Vx) {
	id = Vx.get_id();
	orig = Vx.get_orig();
	index = Vx.get_index();
	lowlink = min(lowlink,Vx.get_lowlink());
    cost += Vx.get_cost();
    tcost += Vx.get_tcost();
	Vxcnt += 1;
  }

  friend ostream&  operator<<(ostream& os, const VxVorAccum& Va) {
	// << "Summary of Vertex attributes by id " <<endl
   // << "Id" << "\t"<< "Stop" << "\t" << " id " << "\t" << " cost "<< "\t" << "TotCost "
	// << "\t"<< "ParCnt"<< "\t"<<"index "<< "\t"<<"LowLink "<< endl
    return os << Va.id<< "\t" << Va.orig << "\t" << Va.cost << "\t" << Va.tcost <<  "\t" 
		<< Va.Vxcnt<<  "\t" << Va.index <<  "\t" << Va.lowlink<< endl;
  }
	long get_id () {return id;}
	long get_orig () {return orig;}
	double get_cost() {return cost;}
	double get_tcost() {return tcost;}
	double get_index() {return index;}
	double get_lowlink() {return lowlink;}
	long get_Vxcnt() {return Vxcnt;} 

};




class VertAccum {
	long id;
	unsigned long pid;  // vertex id
	         long pidp; // predecessor
	         short lbl; // label
         	 double cost;  // cost at vertex
         	 double tcost; // total cost
			 long fid;  //  
			 long orig;  // Origin for vertex least cost path
			 short toStop;  // travel direction toStop=1 (on) or from stop=0 (off) 
			 long index;  // index of vertex in the graph of the Strongly Connected Component (SCC)
			 long lowlink;  // indicator of the root index of the SCC, if index=lowlink, this is root
			long Vxcnt; // count 
public:
  VertAccum() : id(0), pid(0),pidp(0), lbl(0), cost(0),tcost(0),fid(0),orig(0), 
	  toStop(0),index(0),lowlink(inf),Vxcnt(0) {}
  void operator()(const vertexp& Vx) {
	id = Vx.get_id();
	orig = Vx.get_orig();
	index = max(index,Vx.get_index());
	lowlink = min(lowlink,Vx.get_lowlink());
    cost += Vx.get_cost();
    tcost += Vx.get_tcost();
	Vxcnt += 1;
  }

  friend ostream&  operator<<(ostream& os, const VertAccum& Va) {
	/*  << "Summary of Vertex attributes by id " <<endl
    << "Id" << "\t"<< "Stop" << "\t" << " id " << "\t" << " cost "<< "\t" << "TotCost "
	<< "\t"<< "ParCnt"<< "\t"<<"index "<< "\t"<<"LowLink "<< endl */
    return os << Va.id<< "\t" << Va.orig << "\t" << Va.cost << "\t" << Va.tcost <<  "\t" 
		<< Va.Vxcnt<< "\t" << Va.index <<  "\t" << Va.lowlink<< endl;
  }
	long get_id () {return id;}
	long get_orig () {return orig;}
	double get_cost() {return cost;}
	double get_tcost() {return tcost;}
	double get_index() {return index;}
	double get_lowlink() {return lowlink;}
	long get_Vxcnt() {return Vxcnt;} 

};



template <typename o, typename m >
 o& readvertex1(string& rec1, o& o1, m& maphdrit, char *seps) // *seps= "\t" 
 {
	int i=0; 
	int ib=10;
	unsigned long vid=0;
	         int vidp;
			 int vp;
	         short vlbl;
         	 double vcost;
         	 double tcost;
         	 Point pt;
   static int rno = 0; //record number

// VertOId,VertFC,PVertOId,PVertFC,OriginOId,OriginFC,OriginName,VertCost,TCost,Labeled,VertNote,OBJECTID

   const int bufsz = 1000; // Buffer size;
//	const char *strx;
   string  f1, fldval, fldhdr;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
//   map< int,  string>::iterator maphdrit;
	rno++;

  mapflds1 = recoread(rec1);

  mapfldsit = mapflds1.begin();

while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;

	fldval = mapfldsit->second;
//    strx = fldval.c_str(); //strtok(stredge,seps);
	fldhdr = maphdrit->second;

   f1 = (stringUpper<string>(fldhdr)); 
   if (f1 ==  "VERTOID" || f1 == "VERTID" || f1 == "OBJECTID")
   {
      vid = fromString<long>(fldval); //strtol(strx,&stop1,ib);
	  if (vid>0) 
	  {
		  o1.set_id(vid);
	  }
	  else
	  {
		  return o1;
	  }
	}
	else if (f1 == "PVERTOID" || f1 == "IDP" )  //predecssor node id
	{
		vidp = fromString<long>(fldval); //strtol(strx,&stop1,ib);
	    o1.set_idp(vidp);
	}
	else if (f1 == "ORIGINOID")  //voronoi attractor/source/facility id
	{
		vp = fromString<long>(fldval); //strtol(strx,&stop1,ib);
			o1.set_idp((int) vp);
			o1.set_orig(vp);
	}
	else if (f1 == "X" || f1 == "XC" || f1 == "XCOORD")  // X-Coordinate
	{
		  vcost = fromString<double>(fldval); //strtod(strx,&stop1);
//		  o1.pt.set_x(vcost);
	      o1.set_x(vcost);
	} 
	else if (f1 == "Y" || f1 == "YC" || f1 == "YCOORD")  // X-Coordinate
	{
		  vcost = fromString<double>(fldval); //strtod(strx,&stop1);
//		  o1.pt.set_y(vcost);
	      o1.set_y(vcost);
	} 
	else if (f1 == "TCOST") // total cost at node (the least cost path from the source to the node)
	{
		tcost = fromString<double>(fldval); //strtod(strx,&stop1);
		o1.set_tcost(tcost);
	}
	else if (f1 == "LABELED") // label 
	{
		if (vid>0) {vlbl = fromString<short>(fldval); //strtol(strx,&stop1,ib); 
	  o1.set_lbl((short) vlbl);}
	}
             
   // get the next fld data
	mapfldsit++;
	maphdrit++;
} // close while
   if (o1.get_orig()>0)
   {
	cout<<"Stop point id "<<vid<<", origin : "<<o1.get_orig()<<", Cost "<<vcost<<endl;
   }
   return o1;
}

	template <class T>
	string  makeyTrip( T& kt)
	{
		return (kt.route() + 'x' + (kt.trip()) + 'x' +  to_string<short> (kt.dirn())) ;
	}

	template<typename U>  
	U from_string(const std::string& s) {
		std::istringstream is(s);
		U u;
		is >> u;
		return u;
	}

	class cInTime {
	private: 
		int mId;
		double mTime;
	public:
		cInTime();
		cInTime (int id, double ctime);  
	};

	class cDblTime {
	private: 
		int mId;
		double mTime;
	public:
		cDblTime();
		cDblTime (double ctime , int id);  
	};

	class HMS {
	private: 
		int mId;
		int hTime;
		int mTime;
		int sTime;
	public:
		HMS();
		HMS (int id, int ht , int mt , int st);  
	};

	template <typename a, typename b,typename o,  typename c, typename d> 
	c& stopTableData(a& netdb ,b& sql,o& o1, c& mapData,d& logFile)
	{
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows;
		int columns;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		o* pStop;
		typedef pair <int,o> sopair;
		string svertid="", sstopid="", sFldVal="";
		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( ret != SQLITE_OK )
			return mapData;
		if ( rows < 1 )
			;
		else
		{
			// query the stop detail data from the stop table
			mapData.clear();
			fprintf(logFile," Sql Query %s rows %d cols %d !",sql.c_str(),rows,columns);
			for ( i = 1; i <= rows; i++ )
			{
				pStop = &o1;
				for (j = 0; j <= columns-1; j++ )
				{
					if (j==0) // stopseq
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						pStop->set_StOrdr(iFldVal);

					} else if (j==1) // agency id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setAgencyId(sFldVal);
					} else if (j==2) // route name 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setRouteName(sFldVal);
					} else if (j==3) // route_id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setRouteId(sFldVal);
					} else if (j==4) // route 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setRoute(sFldVal);
					} else if (j==5) // Route Length
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setRouteLength(dFldVal);
					} else if (j==6) // Trip Id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setTripId(sFldVal);
					} else if (j==7) // stop Id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						pStop->setid(iFldVal);
					} else if (j==8) // stop Name
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setStopName(sFldVal);
					} else if (j==9) // stop Code
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setStopCode(sFldVal);
					} else if (j==10) // Arrival Time
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setArrTime(sFldVal);
					} else if (j==11) // Stop Latitude 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setStopLat(sFldVal);
					} else if (j==12) // Stop Longitude 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setStopLong(sFldVal);
					}  else if (j==13) // Parent Station 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setParentStn(sFldVal);
					}  else if (j==14) // Direction Id 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						pStop->setDirn(iFldVal);
					}  else if (j==15) // Trip Block Id 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setTripBlkId(sFldVal);
					}  else if (j==16) // Trip Shape Id 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setTripShpId(sFldVal);
					}   else if (j==17) // Trip Head Sign 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setTripHdSign(sFldVal);
					}   else if (j==18) // Percent Distance  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->posAlong(dFldVal);
					}   else if (j==19) // Cumulative distance 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setCumDist(dFldVal);
					}   else if (j==20) // X - Coordinate
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setxc(dFldVal);
					}    else if (j==21) // Y - Coordinate
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setyc(dFldVal);
					} 

				} // loop over cols
				if (pStop->stOrdr())
				{
					mapData.insert(sopair(pStop->stOrdr(),o1));
				}
			} // loop over rows
			sqlite3_free_table( results );
		} // if rows are found
		return mapData;
	}

	template <typename a, typename b,typename o,  typename c,typename k,typename g,typename p,typename t,typename r, typename d> 
	c& stopTableRTD1(a& netdb ,b& sql,o& o1,c& m1,k& kstop,g& gc,p& ph,t& tp,r& srid, d& logFile)
	{
		int ret=0 , i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		sqlite3_stmt *stmt2 = NULL;
		int rows=0, columns=0 , col=0;
		char *zErrMsg = NULL;
		string sstopi="",strCDim = "xy", strGType="Point";
		bool blnRet = false;
		int pkCount = 0 , fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0, iGDim=2, iMeas=0, iSrid=srid;
		long lFldVal=0;
		double dFldVal=0,mcumDist=0,mpctDist=0,dschlTime=0, dDwell=0,dpDwell=0, dcumDwell=0, drunTimeC=0,dpundRdTm=0;
		double mRteLen=0,dpprobStoph=0,dpDwellTm=0,dpArrDelay=0,dpDepDelay=0,dpArrTm=0,dpDepTm=0,dpsegRdTm=0,dpCRdTm=0;
		// dwell - dwell delay from boarding and alighting cum (DepTime - ArrTime) minutes
		// dcumDwell - cumulative dwell delay from begining of route to the current stop minutes
		// dsegRunTm = dArrTm(current Stop) - dpDepTm(Prev Stop) - Dwell(Prev Stop)
		o* pStop;
		b tblStop="",tbLine="", tblBuf="", sqlQ="",tblBufx2="",stokens="", sstopid="", sFldVal="";
		typedef pair <int,o> sopair;
		typedef c::iterator InpIt; 
		InpIt stIt;

		vector<string> svecHMS, svecHAMPM;	
		vector<int> ivecHMS;	
		double dschlTm=0, dacArrTm =0,dwellTm =0,drunTm =0, dacDepTm=0,dundelTm=0, dschTmPt = 0, dbegArrTm=0, dbegDepTm=0,dpDepVol=0; 
		if (kstop.hOnSum() < kstop.hOffSum()) {
			dpDepVol = kstop.hOffSum() - kstop.hOnSum();
		}
		// drop the old table named without the direction indicator
		tblStop = "s" + kstop.route()+ kstop.schlName()+ kstop.dir()+ to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
    	ReplaceAll2(tblStop,"__","_");

		sqlQ="Drop Table " + tblStop + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}

		sqlQ="Drop Table " + tblStop +"_Buf" + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		
		tblStop = "s" + kstop.route()+ kstop.schlName()+ kstop.dir()+ to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
    	ReplaceAll2(tblStop,"__","_");

		//replace(tblStop.begin(),tblStop.end(),'__','_');
		// Create a stop table for this run if doesn't exist already if it does drop and recreate it
		sqlQ = " as " + sql;
		blnRet = createSpaTbl (sqlQ,tblStop,netdb,srid,logFile);
		if ( blnRet )
		{  // recover geometry
			// first get teh geometry type and dimension details for the created table
			sqlQ = "SELECT Count(*) NumObj, GeometryType(\"Geometry\") GType, Srid(\"Geometry\") Srid, CoordDimension(\"Geometry\") CoorDim, ST_NDims(\"Geometry\") Dim , St_IsMeasured(\"Geometry\") Measured" 
				" from " + tblStop + " " ;
			if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
			{
					//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
					logFile<<"Table "<<tblStop<< " Geometric information could not be read !"<<endl;// assume generic
			} else { 

				// read tabel geoemtric information if it is available information
				rows=0;

				while ( sqlite3_step( stmt ) == SQLITE_ROW )
				{
					rows++;
					// query the stop geometry detail data from the stop table
					col=0; // count of geometries
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
					}
					col = 1;  // Geometric type "Point", "LineString", "PolyGon" , etc
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strGType = sFldVal;
					}
					col=2; // SRID 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iSrid = iFldVal;
					}
					col=3; // Coordinate Dimension as a string XY, XYM, XYZ, XYZM
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strCDim = sFldVal;
					}
					col=4; // Coordinate Dimension as an integer 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iGDim = iFldVal;
					}
					col=5; // If this is measured data 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iMeas = iFldVal;
					}
				}
			}
			if (rows) {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
					"'" + strGType +"' , "  + to_string<int>(iGDim) + " );";
			} else  {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " ,'Point',2);";
			}
			blnRet = recoverSpatialGeometry(sqlQ,tblStop,netdb,srid,logFile);
			// create a line shape using the stop tbl just created to compute the cumulative distance to each stop 
			tbLine =  tblStop + "_RteGeom";
			sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
				" makeline(casttomultipoint( Geometry))  Geometry from " ;
			sqlQ.append(tblStop) ;
			sqlQ.append (" group by SchlName, RteName, DirName, TimePeriod, TripId;" );
			blnRet = createSpaTbl (sqlQ,tbLine,netdb,srid,logFile);

			if (blnRet) {// recover spatial geometry
				if (rows) {
					sqlQ = " SELECT recovergeometrycolumn('" +  tbLine  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
						" 'Linestring', "  + to_string<int>(iGDim) + " );";
				} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tbLine + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Linestring',2);";
				}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
			// create a buffer polygon using the stop tbl just created to select the demand distribution area 
				tblBuf =  tblStop + "_Buf";
				sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer(( Geometry) , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tbLine + "  "
					" group by SchlName, RteName, DirName, TimePeriod, TripId;" ;
				blnRet = createSpaTbl (sqlQ,tblBuf,netdb,srid,logFile);
				if ( blnRet )
				{ // recover spatial geometry
					if (rows) {
						sqlQ = " SELECT recovergeometrycolumn('" +  tblBuf  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
							" 'Polygon', "  + to_string<int>(iGDim) + " );";
					} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tblBuf + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
					}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
					// create a second buffer polygon using the polygon just created to select the edge distribution area 
					tblBufx2 = tblBuf + "x2" ;
					// Create a stop table for this run if doesn't exist already if it does drop an recreate it
					sqlQ  = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer( Geometry , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tblBuf + " " ;
					blnRet = createSpaTbl (sqlQ,tblBufx2,netdb,srid,logFile);
					if (blnRet) {
						if (rows) {
							sqlQ = " SELECT recovergeometrycolumn('" +  tblBufx2  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
								" 'Polygon', "  + to_string<int>(iGDim) + " );";
						} else  {
							sqlQ = " SELECT recovergeometrycolumn('" + tblBufx2 + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
						}
						blnRet = recoverSpatialGeometry(sqlQ,tblBufx2,netdb,srid,logFile);
					}
				} else
				{
				  // some error occurred
				  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				  logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
				}
			} // if make line succeeds then create the buffer polygons 
			else {
				  // some error occurred
				fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
			}
		} else
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
		  //fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"\nSQLite error: "<<sqlQ <<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		}

		// query the stop table for this run and store it in the stop table

//		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"SQLite prepare error: \tSQL error: "<< sql<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		  return m1;
		}
		m1.clear();

		// query the stop detail data from the stop table
		rows=0;
//		if(sqlite3_data_count(stmt) > 0 ) { 

			while ( sqlite3_step( stmt ) == SQLITE_ROW )
			{
				rows++;
				pStop = &o1;
				// query the stop detail data from the stop table
				col=0;
				iFldVal = sqlite3_column_int(stmt, col);
				if (iFldVal>0) {  // id
					pStop->set_id(rows);
					pStop->set_StOrdr(iFldVal);
				}
				// j==1 Schl Name
				// j==2 route name 
				// j==3 DirName
				// j==4 Time period 
				// j==5 set trip start time
				// j==6 Trip Key or Id
				col = 6;
				iFldVal = sqlite3_column_int(stmt, col);
				pStop->settripId(iFldVal);
				// j==7 Trip Number
				// j==8 schduled time  
				col = 8 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dschlTime = triphrs(svecHMS,iFldVal,dschlTime); 
					if (svecHAMPM.at(1) == "PM" && dschlTime < 12.0) {
						dschlTime += 12 ; 
					}
				}			//pStop->setschlTime(sFldVal);
				// j==9 Actual Arrival Time / Cumulative Ride time to stop
				col = 9 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dacArrTm = from_string<double>(sFldVal); 
					if (rows==1 )  {
						dbegArrTm = dacArrTm;
						dpArrTm = dacArrTm;
					}
						drunTm = (dacArrTm - dbegArrTm) ;  
				}			//pStop->setschlTime(sFldVal);
				// j==10 get Actual Departure time = actual arrival time + dwell time
				col = 10 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dundelTm = from_string<double>(sFldVal); 
					if (rows==1 )  {
						dbegDepTm = dundelTm;
						dpDepTm = dundelTm;
					}

				}			//pStop->setschlTime(sFldVal);
				// j==11 CumDist  
				col = 11 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dFldVal= from_string<double> (sFldVal);
					pStop->set_CumDist(dFldVal);
				} 
				// j==12 Stop Id
				col = 12;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopLbl(sFldVal);
					// calculate the cumulative distance from the begining of the route using the line geometry  
					//sqlQ = " select t1.Stop_Id , "
					//" st_length(st_splitLeft(makeline(casttomultipoint(t1.Geometry)),st_PointN(makeline(casttomultipoint(t1.Geometry)),t1.Id)))/5280 SplitLen, "
					//" ST_length(t2.Geometry) RteLength "
					sqlQ = " select t1.Stop_Id , ST_Line_Locate_Point(t2.Geometry, t1.Geometry) PctDist, st_length(t2.Geometry)/5280 RteLength "
						", st_length(t2.Geometry) * ST_Line_Locate_Point(t2.Geometry, t1.Geometry) CumDist "  
							" from " + tblStop + " t1 , " + tbLine + " t2  where t1.Stop_Id like '" + sFldVal + "' ; " ;
					// query the data 
					if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt2, NULL ) != SQLITE_OK )
					{
					  // some error occurred
					  fprintf(stderr, "\nStop Distance Query error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
					  logFile<<"Stop Distance Query: \tSQL error: "<< sqlQ<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
					  return m1;
					}
					// j==1 PctDist 
					while ( sqlite3_step( stmt2 ) == SQLITE_ROW )
					{
						col = 1;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
						mpctDist = sqlite3_column_double(stmt2, col);
						}
						// j==2 Rte Length  
						col = 2;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mRteLen = sqlite3_column_double(stmt2, col);
						}  
						// j==3 Cum Length  
						col = 3;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mcumDist = sqlite3_column_double(stmt2, col);
						}  
					}  
					mcumDist = mRteLen * mpctDist;
					
						//st_PointN(makeline(casttomultipoint(t1.Geometry)),3)) Point3,
						//st_length(makeline(casttomultipoint(t1.Geometry))) TotalLeng
						//pStop->set_CumDist(mcumDist);

				}
				// j==13 Main and Cross Street Intersection 
				col = 13;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopName(sFldVal);
				}
				// j==14 Trip Boardings  
				col = 14;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOns(dFldVal);
					pStop->set_Ons(dFldVal);
				}  
				// j==15 Trip alightings 
				col = 15;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOffs(dFldVal);
					pStop->set_Offs(dFldVal);
				}
				// j==16 location latitude , xc
				col = 16;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_xc(dFldVal);
				}   
				// j==17 location longitude , yc 
				col = 17;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_yc(dFldVal);
				} 
				// j==18 EdgeID 
				col = 18;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					lFldVal= sqlite3_column_int(stmt, col);
					pStop->set_Edgeid(lFldVal);
				} 
				// j==19 Pos Along  Edge  
				col = 19;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_posalong(dFldVal);
				} 
				// j==20 arrDelay  
				col = 20;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_arrDelay(dFldVal); // this is omitted as Arrival Ride time data includes the deceleration delay
					//pStop->set_arrDelay(0);
				} 
				col = 21; // departure delay
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_depDelay(dFldVal); // omitted due to double counting as it is included in the segment time for RTD data
					//pStop->set_depDelay(0);
				} 
				if (sqlite3_column_count(stmt) > 22) {
					col = 22; // schl run time
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dFldVal= sqlite3_column_double(stmt, col);
						pStop->set_nearDist(dFldVal);
					} 
				}
				if (sqlite3_column_count(stmt) > 23) {
					col = 23; // Historic indicator
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnHist(iFldVal); // historic stop?
					} 
				}
				if (sqlite3_column_count(stmt) > 24) {
					col = 24; // Dwell Delay
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dwellTm= sqlite3_column_double(stmt, col);
						pStop->set_dwellDelay(dFldVal);
					} 
				}
				if (sqlite3_column_count(stmt) > 25) {
					col = 25; // segment ride time
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dpsegRdTm= sqlite3_column_double(stmt, col);
						pStop->set_rideDelay(dpsegRdTm);
					} 
				}
				if (pStop->get_id())
				{
					if(rows>1) { // compute the ride time using the ons and offs
						pStop->cdwlDelayh(gc->get_unitontm(),gc->get_unitofftm(),ph.get_hdway());
						dacDepTm = dacArrTm + dwellTm;
						dpsegRdTm = (dpsegRdTm)*60 ; // run time in minutes between prev stop and current stop
						if (dpsegRdTm<=0) { // wrong data
							fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dpDepTm );
							logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Prev. Stop Dep. Tm. " << dpDepTm <<endl;
							//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",);
						}
					}
					if (dacArrTm > dacDepTm) { // wrong data
						fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Dep. Tm. " << dacDepTm <<endl;
						dacDepTm = dacArrTm + pStop->get_dwellDelay()/3600 ;
					}
					if (tp) {
						// compute the dwell time and historic run time using ons and offs
						//pStop->set_dwellDelay(max(pStop->get_dwellDelay(),(dacDepTm-dacArrTm)*60));
						pStop->cprobstoph(ph.get_hdway());
						if(rows>1) { // compute the ride time using the ons and offs
							pStop->set_CRdTm( dpCRdTm + dpsegRdTm+ dpDwellTm +(dpDepDelay/60)*dpprobStoph + (pStop->get_arrDelay()/60)*pStop->get_probStoph()); // drunTm*60);
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
							pStop->set_dwellDelay(0);
						}
					} else {
						if(rows>1) { // compute the ride time using the ons and offs
						pStop->set_CRdTm(drunTm*60); // actual arrival time from beg of route in mins 
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
						}
					}
					drunTimeC = pStop->get_CRdTm();
					if (rows==1) {
						pStop->set_dwellDelay(0);
						dpDwellTm=0;
						pStop->set_undCRdTm(pStop->get_CRdTm());
						dpArrDelay = 0;
						dpDepTm = dacDepTm;
					} else {
						pStop->set_undCRdTm(dpundRdTm+dpsegRdTm);  // add the segment ride time including the unaccounted delay at prev stop
						dDwell = pStop->get_dwellDelay()/60; // (max<double>((dacDepTm - dacArrTm),pStop->get_dwellDelay()/60)*60);  // no trips / hr = (60/ph.get_hdway()
						dpArrDelay = pStop->get_arrDelay();
					}
					dpDepVol += (pStop->get_HistOns() - pStop->get_HistOffs());
					pStop->set_HistDepVol(dpDepVol);
					m1.insert(sopair(pStop->get_id(),o1));
					dpDwellTm = pStop->get_dwellDelay()*60 ; // minutes
					dpDwell = (dacDepTm-dacArrTm)*60;
					dDwell = max<double>(0.0,(dpDwell-dpDwellTm));  // dwell delay that is not accounted by demand
					dpsegRdTm = (dacArrTm - dpDepTm)*60 + dDwell ;  // segement ride time
					dcumDwell += dpDwellTm + dDwell;
					dpDepDelay = pStop->get_depDelay();
					dpprobStoph = pStop->get_probStoph();
					dpArrTm = dacArrTm;
					dpDepTm = dacDepTm;
					dpundRdTm = pStop->get_undCRdTm();
					dpCRdTm = pStop->get_CRdTm();
				}
			} // loop over all rows 
			// update the last stop depDelay (0) and probStopping (1)  
			if (m1.size()>1) {
				InpIt begin= m1.begin();
				InpIt end = m1.end();
				end--;
				end = m1.find(end->first);
				if(end!=m1.end()) { 
					o1 = end->second;
					o1.set_probStop(1.0);
					o1.set_depDelay(0.0);
					m1.erase(end); 
					m1.insert(sopair(o1.get_id(),o1));
				}
			}
	//	} // if data count is greater than 0
		  // there are no rows to fetch
		sqlite3_finalize( stmt );
		return m1;
	}





	template <typename a, typename b,typename o,  typename c,typename k,typename g,typename p,typename t,typename r, typename d, typename e> 
	c& stopTableRTD2(a& netdb ,b& sql,o& o1,c& m1,k& kstop,g& gc,p& ph,t& tp,r& srid, d& logFile, e& elist)
	{
		int ret=0 , i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		sqlite3_stmt *stmt2 = NULL;
		int rows=0, columns=0 , col=0;
		char *zErrMsg = NULL;
		string sstopi="",strCDim = "xy", strGType="Point";
		bool blnRet = false;
		int pkCount = 0 , fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0, iGDim=2, iMeas=0, iSrid=srid;
		long lFldVal=0;
		double dFldVal=0,mcumDist=0,mpctDist=0,dschlTime=0, dDwell=0,dpDwell=0, dcumDwell=0, drunTimeC=0,dpundRdTm=0;
		double mRteLen=0,dpprobStoph=0,dpDwellTm=0,dpArrDelay=0,dpDepDelay=0,dpArrTm=0,dpDepTm=0,dpsegRdTm=0,dpCRdTm=0,drDelay=0;
		// dwell - dwell delay from boarding and alighting cum (DepTime - ArrTime) minutes
		// dcumDwell - cumulative dwell delay from begining of route to the current stop minutes
		// dsegRunTm = dArrTm(current Stop) - dpDepTm(Prev Stop) - Dwell(Prev Stop)
		o* pStop;
		b tblStop="",tbLine="", tblBuf="", sqlQ="",tblBufx2="",stokens="", sstopid="", sFldVal="";
		typedef pair <int,o> sopair;
		typedef c::iterator InpIt; 
		InpIt stIt;

		vector<string> svecHMS, svecHAMPM;	
		vector<int> ivecHMS;	
		double dschlTm=0, dacArrTm =0,dwellTm =0,drunTm =0, dacDepTm=0,dundelTm=0, dschTmPt = 0, dbegArrTm=0, dbegDepTm=0,dpDepVol=0; 
		if (kstop.hOnSum() < kstop.hOffSum()) {
			dpDepVol = kstop.hOffSum() - kstop.hOnSum();
		}
		// prepare the table name 
		tblStop = "s" + kstop.route()+ kstop.schlName()+ kstop.dir()+ to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
		ReplaceAll2(tblStop,"__","_");
		//sqlQ="Drop Table " + tblStop + "; "; 
		//if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		//{
		//		//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
		//		logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		//}

		sqlQ="Drop Table " + tblStop +"_Buf" + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		
		//replace(tblStop.begin(),tblStop.end(),'__','_');
		// Create a stop table for this run if does exist already drop and recreate it
		sqlQ = " as " + sql;
		blnRet = createSpaTbl (sqlQ,tblStop,netdb,srid,logFile);
		if ( blnRet )
		{  // recover geometry
			// first get the geometry type and dimension details for the created table
			sqlQ = "SELECT Count(*) NumObj, GeometryType(\"Geometry\") GType, Srid(\"Geometry\") Srid, CoordDimension(\"Geometry\") CoorDim, ST_NDims(\"Geometry\") Dim , St_IsMeasured(\"Geometry\") Measured" 
				" from " + tblStop + " " ;
			if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
			{
					//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " Geometric information could not be read !"<<endl<<"Error Msg : "<<sqlite3_errmsg( netdb)<<endl;// assume generic
			} else { 

				// read geoemtric information if it is available information
				rows=0;

				while ( sqlite3_step( stmt ) == SQLITE_ROW )
				{
					rows++;
					// query the stop geometry detail data from the stop table
					col=0; // count of geometries
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
					}
					col = 1;  // Geometric type "Point", "LineString", "PolyGon" , etc
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strGType = sFldVal;
					}
					col=2; // SRID 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iSrid = iFldVal;
					}
					col=3; // Coordinate Dimension as a string XY, XYM, XYZ, XYZM
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strCDim = sFldVal;
					}
					col=4; // Coordinate Dimension as an integer 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iGDim = iFldVal;
					}
					col=5; // If this is measured data 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iMeas = iFldVal;
					}
				}
			}
			if (rows) {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
					"'" + strGType +"' , "  + to_string<int>(iGDim) + " );";
			} else  {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " ,'Point',2);";
			}
			blnRet = recoverSpatialGeometry(sqlQ,tblStop,netdb,srid,logFile);
			// create a line shape using the stop tbl just created to compute the cumulative distance to each stop 
			tbLine =  tblStop + "_RteGeom";
			sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
				" makeline(casttomultipoint( Geometry))  Geometry from " ;
			sqlQ.append(tblStop) ;
			sqlQ.append (" group by SchlName, RteName, DirName, TimePeriod, TripId;" );
			blnRet = createSpaTbl (sqlQ,tbLine,netdb,srid,logFile);

			if (blnRet) {// recover spatial geometry
				if (rows) {
					sqlQ = " SELECT recovergeometrycolumn('" +  tbLine  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
						" 'Linestring', "  + to_string<int>(iGDim) + " );";
				} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tbLine + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Linestring',2);";
				}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
			// create a buffer polygon using the stop tbl just created to select the demand distribution area 
				tblBuf =  tblStop + "_Buf";
				sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer(( Geometry) , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tbLine + "  "
					" group by SchlName, RteName, DirName, TimePeriod, TripId;" ;
				blnRet = createSpaTbl (sqlQ,tblBuf,netdb,srid,logFile);
				if ( blnRet )
				{ // recover spatial geometry
					if (rows) {
						sqlQ = " SELECT recovergeometrycolumn('" +  tblBuf  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
							" 'Polygon', "  + to_string<int>(iGDim) + " );";
					} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tblBuf + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
					}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
					// create a second buffer polygon using the polygon just created to select the edge distribution area 
					tblBufx2 = tblBuf + "x2" ;
					// Create a stop table for this run if doesn't exist already if it does drop an recreate it
					sqlQ  = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer( Geometry , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tblBuf + " " ;
					blnRet = createSpaTbl (sqlQ,tblBufx2,netdb,srid,logFile);
					if (blnRet) {
						if (rows) {
							sqlQ = " SELECT recovergeometrycolumn('" +  tblBufx2  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
								" 'Polygon', "  + to_string<int>(iGDim) + " );";
						} else  {
							sqlQ = " SELECT recovergeometrycolumn('" + tblBufx2 + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
						}
						blnRet = recoverSpatialGeometry(sqlQ,tblBufx2,netdb,srid,logFile);
					}
				} else
				{
				  // some error occurred
				  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				  logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
				}
			} // if make line succeeds then create the buffer polygons 
			else {
				  // some error occurred
				fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
			}
		} else
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
		  //fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"\nSQLite error: "<<sqlQ <<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		}

		// query the stop table for this run and store it in the stop table

//		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"SQLite prepare error: \tSQL error: "<< sql<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		  return m1;
		}
		m1.clear();

		// query the stop detail data from the stop table
		rows=0;
//		if(sqlite3_data_count(stmt) > 0 ) { 

			while ( sqlite3_step( stmt ) == SQLITE_ROW )
			{
				rows++;
				pStop = &o1;
				// query the stop detail data from the stop table
				col=0;
				iFldVal = sqlite3_column_int(stmt, col);
				if (iFldVal>0) {  // id
					pStop->set_id(rows);
					pStop->set_StOrdr(iFldVal);
				}
				// j==1 Schl Name
				// j==2 route name 
				// j==3 DirName
				// j==4 Time period 
				// j==5 set trip start time
				// j==6 Trip Key or Id
				col = 6;
				iFldVal = sqlite3_column_int(stmt, col);
				pStop->settripId(iFldVal);
				// j==7 Trip Number
				// j==8 schduled time  
				col = 8 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dschlTime = triphrs(svecHMS,iFldVal,dschlTime); 
					if (svecHAMPM.at(1) == "PM" && dschlTime < 12.0) {
						dschlTime += 12 ; 
					}
				}			//pStop->setschlTime(sFldVal);
				// j==9 Actual Arrival Time / Cumulative Ride time to stop
				col = 9 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dacArrTm = from_string<double>(sFldVal); 
					if (rows==1 )  {
						dbegArrTm = dacArrTm;
						dpArrTm = dacArrTm;
					}
						drunTm = (dacArrTm - dbegArrTm) ;  
				}			//pStop->setschlTime(sFldVal);
				// j==10 get Actual Departure time = actual arrival time + dwell time
				col = 10 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dundelTm = from_string<double>(sFldVal); 
					if (rows==1 )  {
						dbegDepTm = dundelTm;
						dpDepTm = dundelTm;
					}

				}			//pStop->setschlTime(sFldVal);
				// j==11 CumDist  
				col = 11 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dFldVal= from_string<double> (sFldVal);
					pStop->set_CumDist(dFldVal);
				} 
				// j==12 Stop Id
				col = 12;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopLbl(sFldVal);
					// calculate the cumulative distance from the begining of the route using the line geometry  
					//sqlQ = " select t1.Stop_Id , "
					//" st_length(st_splitLeft(makeline(casttomultipoint(t1.Geometry)),st_PointN(makeline(casttomultipoint(t1.Geometry)),t1.Id)))/5280 SplitLen, "
					//" ST_length(t2.Geometry) RteLength "
					sqlQ = " select t1.Stop_Id , ST_Line_Locate_Point(t2.Geometry, t1.Geometry) PctDist, st_length(t2.Geometry)/5280 RteLength "
						", st_length(t2.Geometry) * ST_Line_Locate_Point(t2.Geometry, t1.Geometry) CumDist "  
							" from " + tblStop + " t1 , " + tbLine + " t2  where t1.Stop_Id like '" + sFldVal + "' ; " ;
					// query the data 
					if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt2, NULL ) != SQLITE_OK )
					{
					  // some error occurred
					  fprintf(stderr, "\nStop Distance Query error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
					  logFile<<"Stop Distance Query: \tSQL error: "<< sqlQ<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
					  return m1;
					}
					// j==1 PctDist 
					while ( sqlite3_step( stmt2 ) == SQLITE_ROW )
					{
						col = 1;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
						mpctDist = sqlite3_column_double(stmt2, col);
						}
						// j==2 Rte Length  
						col = 2;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mRteLen = sqlite3_column_double(stmt2, col);
						}  
						// j==3 Cum Length  
						col = 3;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mcumDist = sqlite3_column_double(stmt2, col);
						}  
					}  
					mcumDist = mRteLen * mpctDist;
					
						//st_PointN(makeline(casttomultipoint(t1.Geometry)),3)) Point3,
						//st_length(makeline(casttomultipoint(t1.Geometry))) TotalLeng
						//pStop->set_CumDist(mcumDist);

				}
				// j==13 Main and Cross Street Intersection 
				col = 13;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopName(sFldVal);
				}
				// j==14 Trip Boardings  
				col = 14;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOns(dFldVal);
					pStop->set_Ons(dFldVal);
				}  
				// j==15 Trip alightings 
				col = 15;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOffs(dFldVal);
					pStop->set_Offs(dFldVal);
				}
				// j==16 location latitude , xc
				col = 16;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_xc(dFldVal);
				}   
				// j==17 location longitude , yc 
				col = 17;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_yc(dFldVal);
				} 
				// j==18 EdgeID 
				col = 18;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					lFldVal= sqlite3_column_int(stmt, col);
					pStop->set_Edgeid(lFldVal);
				} 
				// j==19 Pos Along  Edge  
				col = 19;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_posalong(dFldVal);
				} 
				// j==20 arrDelay  
				col = 20;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_arrDelay(dFldVal); // this is omitted as Arrival Ride time data includes the deceleration delay
					//pStop->set_arrDelay(0);
				} 
				col = 21; // departure delay
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_depDelay(dFldVal); // omitted due to double counting as it is included in the segment time for RTD data
					//pStop->set_depDelay(0);
				} 
				if (sqlite3_column_count(stmt) > 22) {
					col = 22; // schl run time
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dFldVal= sqlite3_column_double(stmt, col);
					} 
				}
				if (sqlite3_column_count(stmt) > 23) {
					col = 23; // Historic indicator
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnHist(iFldVal); // historic stop?
					} 
				}
				if (sqlite3_column_count(stmt) > 24) {
					col = 24; // Dwell Delay
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dwellTm= sqlite3_column_double(stmt, col);
						pStop->set_dwellDelay(dwellTm);
					} 
				}
				if (sqlite3_column_count(stmt) > 25) {
					col = 25; // ride Delay to thru passengers
					if (sqlite3_column_bytes(stmt, col) !=0) {
						drDelay= sqlite3_column_double(stmt, col);
						pStop->set_rideDelay(drDelay);
					} 
				}
				if (sqlite3_column_count(stmt) > 26) {
					col = 26; // eliminatable during segment ride time
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnElim(iFldVal); // Eliminate in 1 yes, 0 no DP?
						if (iFldVal == 0 ) {elist.insert(pStop->get_id());}
					} 
				}
				if (pStop->get_id())
				{
					if(rows>1) { // compute the ride time using the ons and offs
						pStop->cdwlDelayh(gc->get_unitontm(),gc->get_unitofftm(),ph.get_hdway());
						dacDepTm = dacArrTm + dwellTm;
						dpsegRdTm = (dacArrTm - dpDepTm)*60 ; // run time in minutes between prev stop and current stop
						if (dpsegRdTm<=0) { // wrong data
							fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dpDepTm );
							logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Prev. Stop Dep. Tm. " << dpDepTm <<endl;
							//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",);
						}
					}
					if (dacArrTm > dacDepTm) { // wrong data
						fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Dep. Tm. " << dacDepTm <<endl;
						dacDepTm = dacArrTm + pStop->get_dwellDelay()/3600 ;
					}
					if (tp) {
						// compute the dwell time and historic run time using ons and offs
						//pStop->set_dwellDelay(max(pStop->get_dwellDelay(),(dacDepTm-dacArrTm)*60));
						pStop->cprobstoph(ph.get_hdway());
						if(rows>1) { // compute the ride time using the ons and offs
							pStop->set_CRdTm( dpCRdTm + dpsegRdTm+ dpDwellTm +(dpDepDelay/60)*dpprobStoph + (pStop->get_arrDelay()/60)*pStop->get_probStoph()); // drunTm*60);
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
							pStop->set_dwellDelay(0);
						}
					} else {
						if(rows>1) { // compute the ride time using the ons and offs
						pStop->set_CRdTm(drunTm*60); // actual arrival time from beg of route in mins 
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
						}
					}
					drunTimeC = pStop->get_CRdTm();
					if (rows==1) {
						pStop->set_dwellDelay(0);
						dpDwellTm=0;
						pStop->set_undCRdTm(pStop->get_CRdTm());
						dpArrDelay = 0;
						dpDepTm = dacDepTm;
					} else {
						pStop->set_undCRdTm(dpundRdTm+dpsegRdTm);  // add the segment ride time including the unaccounted delay at prev stop
						dDwell = pStop->get_dwellDelay()/60; // (max<double>((dacDepTm - dacArrTm),pStop->get_dwellDelay()/60)*60);  // no trips / hr = (60/ph.get_hdway()
						dpArrDelay = pStop->get_arrDelay();
					}
					dpDepVol += (pStop->get_HistOns() - pStop->get_HistOffs());
					pStop->set_HistDepVol(dpDepVol);
					m1.insert(sopair(pStop->get_id(),o1));
					dpDwellTm = pStop->get_dwellDelay()*60 ; // minutes
					dpDwell = (dacDepTm-dacArrTm)*60;
					dDwell = max<double>(0.0,(dpDwell-dpDwellTm));  // dwell delay that is not accounted by demand
					dpsegRdTm = (dacArrTm - dpDepTm)*60 + dDwell ;  // segement ride time
					dcumDwell += dpDwellTm + dDwell;
					dpDepDelay = pStop->get_depDelay();
					dpprobStoph = pStop->get_probStoph();
					dpArrTm = dacArrTm;
					dpDepTm = dacDepTm;
					dpundRdTm = pStop->get_undCRdTm();
					dpCRdTm = pStop->get_CRdTm();
				}
			} // loop over all rows 
			// update the last stop depDelay (0) and probStopping (1)  
			if (m1.size()>1) {
				InpIt begin= m1.begin();
				InpIt end = m1.end();
				end--;
				end = m1.find(end->first);
				if(end!=m1.end()) { 
					o1 = end->second;
					o1.set_probStop(1.0);
					o1.set_depDelay(0.0);
					m1.erase(end); 
					m1.insert(sopair(o1.get_id(),o1));
				}
			}
	//	} // if data count is greater than 0
		  // there are no rows to fetch
		sqlite3_finalize( stmt );
		return m1;
	}







	template <typename a, typename b,typename o,  typename c,typename k,typename g,typename p,typename t,typename r, typename d> 
	c& stopTableRTD2(a& netdb ,b& sql,o& o1,c& m1,k& kstop,g& gc,p& ph,t& tp,r& srid, d& logFile)
	{
		int ret=0 , i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		sqlite3_stmt *stmt2 = NULL;
		int rows=0, columns=0 , col=0;
		char *zErrMsg = NULL;
		string sstopi="",strCDim = "xy", strGType="Point";
		bool blnRet = false;
		int pkCount = 0 , fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0, iGDim=2, iMeas=0, iSrid=srid;
		long lFldVal=0;
		double dFldVal=0,mcumDist=0,mpctDist=0,dschlTime=0, dDwell=0,dpDwell=0, dcumDwell=0, drunTimeC=0,dpundRdTm=0;
		double mRteLen=0,dpprobStoph=0,dpDwellTm=0,dpArrDelay=0,dpDepDelay=0,dpArrTm=0,dpDepTm=0,dpsegRdTm=0,dpCRdTm=0;
		// dwell - dwell delay from boarding and alighting cum (DepTime - ArrTime) minutes
		// dcumDwell - cumulative dwell delay from begining of route to the current stop minutes
		// dsegRunTm = dArrTm(current Stop) - dpDepTm(Prev Stop) - Dwell(Prev Stop)
		o* pStop;
		b tblStop="",tbLine="", tblBuf="", sqlQ="",tblBufx2="",stokens="", sstopid="", sFldVal="";
		typedef pair <int,o> sopair;
		typedef c::iterator InpIt; 
		InpIt stIt;

		vector<string> svecHMS, svecHAMPM;	
		vector<int> ivecHMS;	
		double dschlTm=0, dacArrTm =0,dwellTm =0,drunTm =0, dacDepTm=0,dundelTm=0, dschTmPt = 0, dbegArrTm=0, dbegDepTm=0,dpDepVol=0; 
		if (kstop.hOnSum() < kstop.hOffSum()) {
			dpDepVol = kstop.hOffSum() - kstop.hOnSum();
		}
		// prepare the table name 
		tblStop = "s" + kstop.route()+ kstop.schlName()+ kstop.dir()+ to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
		ReplaceAll2(tblStop,"__","_");
		//sqlQ="Drop Table " + tblStop + "; "; 
		//if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		//{
		//		//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
		//		logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		//}

		sqlQ="Drop Table " + tblStop +"_Buf" + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		
		//replace(tblStop.begin(),tblStop.end(),'__','_');
		// Create a stop table for this run if does exist already drop and recreate it
		sqlQ = " as " + sql;
		blnRet = createSpaTbl (sqlQ,tblStop,netdb,srid,logFile);
		if ( blnRet )
		{  // recover geometry
			// first get the geometry type and dimension details for the created table
			sqlQ = "SELECT Count(*) NumObj, GeometryType(\"Geometry\") GType, Srid(\"Geometry\") Srid, CoordDimension(\"Geometry\") CoorDim, ST_NDims(\"Geometry\") Dim , St_IsMeasured(\"Geometry\") Measured" 
				" from " + tblStop + " " ;
			if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
			{
					//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " Geometric information could not be read !"<<endl<<"Error Msg : "<<sqlite3_errmsg( netdb)<<endl;// assume generic
			} else { 

				// read geoemtric information if it is available information
				rows=0;

				while ( sqlite3_step( stmt ) == SQLITE_ROW )
				{
					rows++;
					// query the stop geometry detail data from the stop table
					col=0; // count of geometries
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
					}
					col = 1;  // Geometric type "Point", "LineString", "PolyGon" , etc
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strGType = sFldVal;
					}
					col=2; // SRID 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iSrid = iFldVal;
					}
					col=3; // Coordinate Dimension as a string XY, XYM, XYZ, XYZM
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strCDim = sFldVal;
					}
					col=4; // Coordinate Dimension as an integer 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iGDim = iFldVal;
					}
					col=5; // If this is measured data 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iMeas = iFldVal;
					}
				}
			}
			if (rows) {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
					"'" + strGType +"' , "  + to_string<int>(iGDim) + " );";
			} else  {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " ,'Point',2);";
			}
			blnRet = recoverSpatialGeometry(sqlQ,tblStop,netdb,srid,logFile);
			// create a line shape using the stop tbl just created to compute the cumulative distance to each stop 
			tbLine =  tblStop + "_RteGeom";
			sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
				" makeline(casttomultipoint( Geometry))  Geometry from " ;
			sqlQ.append(tblStop) ;
			sqlQ.append (" group by SchlName, RteName, DirName, TimePeriod, TripId;" );
			blnRet = createSpaTbl (sqlQ,tbLine,netdb,srid,logFile);

			if (blnRet) {// recover spatial geometry
				if (rows) {
					sqlQ = " SELECT recovergeometrycolumn('" +  tbLine  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
						" 'Linestring', "  + to_string<int>(iGDim) + " );";
				} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tbLine + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Linestring',2);";
				}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
			// create a buffer polygon using the stop tbl just created to select the demand distribution area 
				tblBuf =  tblStop + "_Buf";
				sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer(( Geometry) , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tbLine + "  "
					" group by SchlName, RteName, DirName, TimePeriod, TripId;" ;
				blnRet = createSpaTbl (sqlQ,tblBuf,netdb,srid,logFile);
				if ( blnRet )
				{ // recover spatial geometry
					if (rows) {
						sqlQ = " SELECT recovergeometrycolumn('" +  tblBuf  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
							" 'Polygon', "  + to_string<int>(iGDim) + " );";
					} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tblBuf + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
					}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
					// create a second buffer polygon using the polygon just created to select the edge distribution area 
					tblBufx2 = tblBuf + "x2" ;
					// Create a stop table for this run if doesn't exist already if it does drop an recreate it
					sqlQ  = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer( Geometry , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tblBuf + " " ;
					blnRet = createSpaTbl (sqlQ,tblBufx2,netdb,srid,logFile);
					if (blnRet) {
						if (rows) {
							sqlQ = " SELECT recovergeometrycolumn('" +  tblBufx2  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
								" 'Polygon', "  + to_string<int>(iGDim) + " );";
						} else  {
							sqlQ = " SELECT recovergeometrycolumn('" + tblBufx2 + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
						}
						blnRet = recoverSpatialGeometry(sqlQ,tblBufx2,netdb,srid,logFile);
					}
				} else
				{
				  // some error occurred
				  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				  logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
				}
			} // if make line succeeds then create the buffer polygons 
			else {
				  // some error occurred
				fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
			}
		} else
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
		  //fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"\nSQLite error: "<<sqlQ <<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		}

		// query the stop table for this run and store it in the stop table

//		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"SQLite prepare error: \tSQL error: "<< sql<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		  return m1;
		}
		m1.clear();

		// query the stop detail data from the stop table
		rows=0;
//		if(sqlite3_data_count(stmt) > 0 ) { 

			while ( sqlite3_step( stmt ) == SQLITE_ROW )
			{
				rows++;
				pStop = &o1;
				// query the stop detail data from the stop table
				col=0;
				iFldVal = sqlite3_column_int(stmt, col);
				if (iFldVal>0) {  // id
					pStop->set_id(rows);
					pStop->set_StOrdr(iFldVal);
				}
				// j==1 Schl Name
				// j==2 route name 
				// j==3 DirName
				// j==4 Time period 
				// j==5 set trip start time
				// j==6 Trip Key or Id
				col = 6;
				iFldVal = sqlite3_column_int(stmt, col);
				pStop->settripId(iFldVal);
				// j==7 Trip Number
				// j==8 schduled time  
				col = 8 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dschlTime = triphrs(svecHMS,iFldVal,dschlTime); 
					if (svecHAMPM.at(1) == "PM" && dschlTime < 12.0) {
						dschlTime += 12 ; 
					}
				}			//pStop->setschlTime(sFldVal);
				// j==9 Actual Arrival Time / Cumulative Ride time to stop
				col = 9 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dacArrTm = from_string<double>(sFldVal); 
					if (rows==1 )  {
						dbegArrTm = dacArrTm;
						dpArrTm = dacArrTm;
					}
						drunTm = (dacArrTm - dbegArrTm) ;  
				}			//pStop->setschlTime(sFldVal);
				// j==10 get Actual Departure time = actual arrival time + dwell time
				col = 10 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dundelTm = from_string<double>(sFldVal); 
					if (rows==1 )  {
						dbegDepTm = dundelTm;
						dpDepTm = dundelTm;
					}

				}			//pStop->setschlTime(sFldVal);
				// j==11 CumDist  
				col = 11 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					dFldVal= from_string<double> (sFldVal);
					pStop->set_CumDist(dFldVal);
				} 
				// j==12 Stop Id
				col = 12;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopLbl(sFldVal);
					// calculate the cumulative distance from the begining of the route using the line geometry  
					//sqlQ = " select t1.Stop_Id , "
					//" st_length(st_splitLeft(makeline(casttomultipoint(t1.Geometry)),st_PointN(makeline(casttomultipoint(t1.Geometry)),t1.Id)))/5280 SplitLen, "
					//" ST_length(t2.Geometry) RteLength "
					sqlQ = " select t1.Stop_Id , ST_Line_Locate_Point(t2.Geometry, t1.Geometry) PctDist, st_length(t2.Geometry)/5280 RteLength "
						", st_length(t2.Geometry) * ST_Line_Locate_Point(t2.Geometry, t1.Geometry) CumDist "  
							" from " + tblStop + " t1 , " + tbLine + " t2  where t1.Stop_Id like '" + sFldVal + "' ; " ;
					// query the data 
					if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt2, NULL ) != SQLITE_OK )
					{
					  // some error occurred
					  fprintf(stderr, "\nStop Distance Query error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
					  logFile<<"Stop Distance Query: \tSQL error: "<< sqlQ<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
					  return m1;
					}
					// j==1 PctDist 
					while ( sqlite3_step( stmt2 ) == SQLITE_ROW )
					{
						col = 1;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
						mpctDist = sqlite3_column_double(stmt2, col);
						}
						// j==2 Rte Length  
						col = 2;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mRteLen = sqlite3_column_double(stmt2, col);
						}  
						// j==3 Cum Length  
						col = 3;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mcumDist = sqlite3_column_double(stmt2, col);
						}  
					}  
					mcumDist = mRteLen * mpctDist;
					
						//st_PointN(makeline(casttomultipoint(t1.Geometry)),3)) Point3,
						//st_length(makeline(casttomultipoint(t1.Geometry))) TotalLeng
						//pStop->set_CumDist(mcumDist);

				}
				// j==13 Main and Cross Street Intersection 
				col = 13;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopName(sFldVal);
				}
				// j==14 Trip Boardings  
				col = 14;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOns(dFldVal);
					pStop->set_Ons(dFldVal);
				}  
				// j==15 Trip alightings 
				col = 15;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOffs(dFldVal);
					pStop->set_Offs(dFldVal);
				}
				// j==16 location latitude , xc
				col = 16;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_xc(dFldVal);
				}   
				// j==17 location longitude , yc 
				col = 17;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_yc(dFldVal);
				} 
				// j==18 EdgeID 
				col = 18;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					lFldVal= sqlite3_column_int(stmt, col);
					pStop->set_Edgeid(lFldVal);
				} 
				// j==19 Pos Along  Edge  
				col = 19;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_posalong(dFldVal);
				} 
				// j==20 arrDelay  
				col = 20;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_arrDelay(dFldVal); // this is omitted as Arrival Ride time data includes the deceleration delay
					//pStop->set_arrDelay(0);
				} 
				col = 21; // departure delay
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_depDelay(dFldVal); // omitted due to double counting as it is included in the segment time for RTD data
					//pStop->set_depDelay(0);
				} 
				if (sqlite3_column_count(stmt) > 22) {
					col = 22; // schl run time
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dFldVal= sqlite3_column_double(stmt, col);
					} 
				}
				if (sqlite3_column_count(stmt) > 23) {
					col = 23; // Historic indicator
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnHist(iFldVal); // historic stop?
					} 
				}
				if (sqlite3_column_count(stmt) > 24) {
					col = 24; // Dwell Delay
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dwellTm= sqlite3_column_double(stmt, col);
						pStop->set_dwellDelay(dFldVal);
					} 
				}
				if (sqlite3_column_count(stmt) > 25) {
					col = 25; // segment ride time
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dpsegRdTm= sqlite3_column_double(stmt, col);
						pStop->set_rideDelay(dpsegRdTm);
					} 
				}
				if (pStop->get_id())
				{
					if(rows>1) { // compute the ride time using the ons and offs
						pStop->cdwlDelayh(gc->get_unitontm(),gc->get_unitofftm(),ph.get_hdway());
						dacDepTm = dacArrTm + dwellTm;
						dpsegRdTm = (dpsegRdTm)*60 ; // run time in minutes between prev stop and current stop
						if (dpsegRdTm<=0) { // wrong data
							fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dpDepTm );
							logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Prev. Stop Dep. Tm. " << dpDepTm <<endl;
							//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",);
						}
					}
					if (dacArrTm > dacDepTm) { // wrong data
						fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Dep. Tm. " << dacDepTm <<endl;
						dacDepTm = dacArrTm + pStop->get_dwellDelay()/3600 ;
					}
					if (tp) {
						// compute the dwell time and historic run time using ons and offs
						//pStop->set_dwellDelay(max(pStop->get_dwellDelay(),(dacDepTm-dacArrTm)*60));
						pStop->cprobstoph(ph.get_hdway());
						if(rows>1) { // compute the ride time using the ons and offs
							pStop->set_CRdTm( dpCRdTm + dpsegRdTm+ dpDwellTm +(dpDepDelay/60)*dpprobStoph + (pStop->get_arrDelay()/60)*pStop->get_probStoph()); // drunTm*60);
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
							pStop->set_dwellDelay(0);
						}
					} else {
						if(rows>1) { // compute the ride time using the ons and offs
						pStop->set_CRdTm(drunTm*60); // actual arrival time from beg of route in mins 
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
						}
					}
					drunTimeC = pStop->get_CRdTm();
					if (rows==1) {
						pStop->set_dwellDelay(0);
						dpDwellTm=0;
						pStop->set_undCRdTm(pStop->get_CRdTm());
						dpArrDelay = 0;
						dpDepTm = dacDepTm;
					} else {
						pStop->set_undCRdTm(dpundRdTm+dpsegRdTm);  // add the segment ride time including the unaccounted delay at prev stop
						dDwell = pStop->get_dwellDelay()/60; // (max<double>((dacDepTm - dacArrTm),pStop->get_dwellDelay()/60)*60);  // no trips / hr = (60/ph.get_hdway()
						dpArrDelay = pStop->get_arrDelay();
					}
					dpDepVol += (pStop->get_HistOns() - pStop->get_HistOffs());
					pStop->set_HistDepVol(dpDepVol);
					m1.insert(sopair(pStop->get_id(),o1));
					dpDwellTm = pStop->get_dwellDelay()*60 ; // minutes
					dpDwell = (dacDepTm-dacArrTm)*60;
					dDwell = max<double>(0.0,(dpDwell-dpDwellTm));  // dwell delay that is not accounted by demand
					dpsegRdTm = (dacArrTm - dpDepTm)*60 + dDwell ;  // segement ride time
					dcumDwell += dpDwellTm + dDwell;
					dpDepDelay = pStop->get_depDelay();
					dpprobStoph = pStop->get_probStoph();
					dpArrTm = dacArrTm;
					dpDepTm = dacDepTm;
					dpundRdTm = pStop->get_undCRdTm();
					dpCRdTm = pStop->get_CRdTm();
				}
			} // loop over all rows 
			// update the last stop depDelay (0) and probStopping (1)  
			if (m1.size()>1) {
				InpIt begin= m1.begin();
				InpIt end = m1.end();
				end--;
				end = m1.find(end->first);
				if(end!=m1.end()) { 
					o1 = end->second;
					o1.set_probStop(1.0);
					o1.set_depDelay(0.0);
					m1.erase(end); 
					m1.insert(sopair(o1.get_id(),o1));
				}
			}
	//	} // if data count is greater than 0
		  // there are no rows to fetch
		sqlite3_finalize( stmt );
		return m1;
	}






	template <typename a, typename b,typename o,  typename c,typename k,typename g,typename p,typename t,typename r, typename d, typename e> 
	c& stopTableRTD(a& netdb ,b& sql,o& o1,c& m1,k& kstop,g& gc,p& ph,t& tp,r& srid, d& logFile,e& elist)
	{
		int ret=0 , i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		sqlite3_stmt *stmt2 = NULL;
		int rows=0, columns=0 , col=0;
		char *zErrMsg = NULL;
		string sstopi="",strCDim = "xy", strGType="Point";
		bool blnRet = false;
		int pkCount = 0 , fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0, iGDim=2, iMeas=0, iSrid=srid;
		long lFldVal=0;
		double dFldVal=0,mcumDist=0,mpctDist=0,dschlTime=0, dDwell=0,dpDwell=0, dcumDwell=0, drunTimeC=0,dpundRdTm=0;
		double mRteLen=0,dpprobStoph=0,dpDwellTm=0,dpArrDelay=0,dpDepDelay=0,dpArrTm=0,dpDepTm=0,dpsegRdTm=0,dpCRdTm=0;
		// dwell - dwell delay from boarding and alighting cum (DepTime - ArrTime) minutes
		// dcumDwell - cumulative dwell delay from begining of route to the current stop minutes
		// dsegRunTm = dArrTm(current Stop) - dpDepTm(Prev Stop) - Dwell(Prev Stop)
		o* pStop;
		b tblStop="",tbLine="", tblBuf="", sqlQ="",tblBufx2="",stokens="", sstopid="", sFldVal="";
		typedef pair <int,o> sopair;
		typedef c::iterator InpIt; 
		InpIt stIt;

		vector<string> svecHMS, svecHAMPM;	
		vector<int> ivecHMS;	
		double dschlTm=0, dacArrTm =0,drunTm =0, dacDepTm=0, dschTmPt = 0, dbegArrTm=0, dbegDepTm=0,dpDepVol=0; 
		if (kstop.hOnSum() < kstop.hOffSum()) {
			dpDepVol = kstop.hOffSum() - kstop.hOnSum();
		}
		// drop the old table named without the direction indicator
		tblStop = "s" + kstop.route()+ kstop.schlName()+  to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
		ReplaceAll2(tblStop,"__","_");

		sqlQ="Drop Table " + tblStop + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		sqlQ="Drop Table " + tblStop +"_Buf" + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		
		tblStop = "s" + kstop.route()+ kstop.schlName()+ kstop.dir()+ to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
		ReplaceAll2(tblStop,"__","_");

		//replace(tblStop.begin(),tblStop.end(),'__','_');
		// Create a stop table for this run if doesn't exist already if it does drop and recreate it
		sqlQ = " as " + sql;
		blnRet = createSpaTbl (sqlQ,tblStop,netdb,srid,logFile);
		if ( blnRet )
		{  // recover geometry
			// first get teh geometry type and dimension details for the created table
			sqlQ = "SELECT Count(*) NumObj, GeometryType(\"Geometry\") GType, Srid(\"Geometry\") Srid, CoordDimension(\"Geometry\") CoorDim, ST_NDims(\"Geometry\") Dim , St_IsMeasured(\"Geometry\") Measured" 
				" from " + tblStop + " " ;
			if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
			{
					//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
					logFile<<"Table "<<tblStop<< " Geometric information could not be read !"<<endl;// assume generic
			} else { 

				// read tabel geoemtric information if it is available information
				rows=0;

				while ( sqlite3_step( stmt ) == SQLITE_ROW )
				{
					rows++;
					// query the stop geometry detail data from the stop table
					col=0; // count of geometries
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
					}
					col = 1;  // Geometric type "Point", "LineString", "PolyGon" , etc
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strGType = sFldVal;
					}
					col=2; // SRID 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iSrid = iFldVal;
					}
					col=3; // Coordinate Dimension as a string XY, XYM, XYZ, XYZM
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strCDim = sFldVal;
					}
					col=4; // Coordinate Dimension as an integer 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iGDim = iFldVal;
					}
					col=5; // If this is measured data 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iMeas = iFldVal;
					}
				}
			}
			if (rows) {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
					"'" + strGType +"' , "  + to_string<int>(iGDim) + " );";
			} else  {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " ,'Point',2);";
			}
			blnRet = recoverSpatialGeometry(sqlQ,tblStop,netdb,srid,logFile);
			// create a line shape using the stop tbl just created to compute the cumulative distance to each stop 
			tbLine =  tblStop + "_RteGeom";
			sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
				" makeline(casttomultipoint( Geometry))  Geometry from " ;
			sqlQ.append(tblStop) ;
			sqlQ.append (" group by SchlName, RteName, DirName, TimePeriod, TripId;" );
			blnRet = createSpaTbl (sqlQ,tbLine,netdb,srid,logFile);

			if (blnRet) {// recover spatial geometry
				if (rows) {
					sqlQ = " SELECT recovergeometrycolumn('" +  tbLine  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
						" 'Linestring', "  + to_string<int>(iGDim) + " );";
				} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tbLine + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Linestring',2);";
				}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
			// create a buffer polygon using the stop tbl just created to select the demand distribution area 
				tblBuf =  tblStop + "_Buf";
				sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer(( Geometry) , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tbLine + "  "
					" group by SchlName, RteName, DirName, TimePeriod, TripId;" ;
				blnRet = createSpaTbl (sqlQ,tblBuf,netdb,srid,logFile);
				if ( blnRet )
				{ // recover spatial geometry
					if (rows) {
						sqlQ = " SELECT recovergeometrycolumn('" +  tblBuf  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
							" 'Polygon', "  + to_string<int>(iGDim) + " );";
					} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tblBuf + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
					}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
					// create a second buffer polygon using the polygon just created to select the edge distribution area 
					tblBufx2 = tblBuf + "x2" ;
					// Create a stop table for this run if doesn't exist already if it does drop an recreate it
					sqlQ  = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer( Geometry , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tblBuf + " " ;
					blnRet = createSpaTbl (sqlQ,tblBufx2,netdb,srid,logFile);
					if (blnRet) {
						if (rows) {
							sqlQ = " SELECT recovergeometrycolumn('" +  tblBufx2  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
								" 'Polygon', "  + to_string<int>(iGDim) + " );";
						} else  {
							sqlQ = " SELECT recovergeometrycolumn('" + tblBufx2 + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
						}
						blnRet = recoverSpatialGeometry(sqlQ,tblBufx2,netdb,srid,logFile);
					}
				} else
				{
				  // some error occurred
				  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				  logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
				}
			} // if make line succeeds then create the buffer polygons 
			else {
				  // some error occurred
				fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
			}
		} else
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
		  //fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"\nSQLite error: "<<sqlQ <<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		}

		// query the stop table for this run and store it in the stop table

//		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"SQLite prepare error: \tSQL error: "<< sql<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		  return m1;
		}
		m1.clear();

		// query the stop detail data from the stop table
		rows=0;
//		if(sqlite3_data_count(stmt) > 0 ) { 

			while ( sqlite3_step( stmt ) == SQLITE_ROW )
			{
				rows++;
				pStop = &o1;
				// query the stop detail data from the stop table
				col=0;
				iFldVal = sqlite3_column_int(stmt, col);
				if (iFldVal>0) {  // id
					pStop->set_id(rows);
					pStop->set_StOrdr(iFldVal);
				}
				// j==1 Schl Name
				// j==2 route name 
				// j==3 DirName
				// j==4 Time period 
				// j==5 set trip start time
				// j==6 Trip Key or Id
				col = 6;
				iFldVal = sqlite3_column_int(stmt, col);
				pStop->settripId(iFldVal);
				// j==7 Trip Number
				// j==8 scheduled start time time
				col = 8 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dschlTime = triphrs(svecHMS,iFldVal,dschlTime); 
					if (svecHAMPM.at(1) == "PM" && dschlTime < 12.0) {
						dschlTime += 12 ; 
					}
				}			//pStop->setschlTime(sFldVal);
				// j==9 Actual Arrival Time
				col = 9 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dacArrTm = triphrs(svecHMS,iFldVal,dacArrTm); 
					if(svecHAMPM.size()>1) { 
						if (svecHAMPM.at(1) == "PM" && dacArrTm < 12.0) {
							dacArrTm  += 12 ; 
						}
					}
					if (rows==1 )  {
						dbegArrTm = dacArrTm;
						dpArrTm = dacArrTm;
					}
						drunTm = (dacArrTm - dbegArrTm) ;  
				}			//pStop->setschlTime(sFldVal);
				// j==10 Actual Departure Time
				col = 10 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dacDepTm = triphrs(svecHMS,iFldVal,dacDepTm); 
					if(svecHAMPM.size()>1) { 
						if (svecHAMPM.at(1) == "PM" && dacDepTm < 12.0) {
							dacDepTm += 12 ; 
						}
					}
					if (rows==1 )  {
						dbegDepTm = dacDepTm;
						dpDepTm = dacDepTm;
					}

				}			//pStop->setschlTime(sFldVal);
				// j==11 Schedule Run Time  
				col = 11 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					iFldVal= from_string<int> (sFldVal);
				} 
				// j==12 Stop Id
				col = 12;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopLbl(sFldVal);
					// calculate the cumulative distance from the begining of the route using the line geometry  
					//sqlQ = " select t1.Stop_Id , "
					//" st_length(st_splitLeft(makeline(casttomultipoint(t1.Geometry)),st_PointN(makeline(casttomultipoint(t1.Geometry)),t1.Id)))/5280 SplitLen, "
					//" ST_length(t2.Geometry) RteLength "
					sqlQ = " select t1.Stop_Id , ST_Line_Locate_Point(t2.Geometry, t1.Geometry) PctDist, st_length(t2.Geometry)/5280 RteLength "
						", st_length(t2.Geometry) * ST_Line_Locate_Point(t2.Geometry, t1.Geometry) CumDist "  
							" from " + tblStop + " t1 , " + tbLine + " t2  where t1.Stop_Id like '" + sFldVal + "' ; " ;
					// query the data 
					if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt2, NULL ) != SQLITE_OK )
					{
					  // some error occurred
					  fprintf(stderr, "\nStop Distance Query error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
					  logFile<<"Stop Distance Query: \tSQL error: "<< sqlQ<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
					  return m1;
					}
					// j==1 PctDist 
					while ( sqlite3_step( stmt2 ) == SQLITE_ROW )
					{
						col = 1;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
						mpctDist = sqlite3_column_double(stmt2, col);
						}
						// j==2 Rte Length  
						col = 2;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mRteLen = sqlite3_column_double(stmt2, col);
						}  
						// j==3 Cum Length  
						col = 3;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mcumDist = sqlite3_column_double(stmt2, col);
						}  
					}  
					mcumDist = mRteLen * mpctDist;
					
						//st_PointN(makeline(casttomultipoint(t1.Geometry)),3)) Point3,
						//st_length(makeline(casttomultipoint(t1.Geometry))) TotalLeng
					//	pStop->set_CumDist(dacDepTm);
						pStop->set_CumDist(mcumDist);

				}
				// j==13 Main and Cross Street Intersection 
				col = 13;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopName(sFldVal);
				}
				// j==14 Trip Boardings  
				col = 14;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOns(dFldVal);
					pStop->set_Ons(dFldVal);
				}  
				// j==15 Trip alightings 
				col = 15;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOffs(dFldVal);
					pStop->set_Offs(dFldVal);
				}
				// j==16 location latitude , xc
				col = 16;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_xc(dFldVal);
				}   
				// j==17 location longitude , yc 
				col = 17;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_yc(dFldVal);
				} 
				// j==18 EdgeID 
				col = 18;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					lFldVal= sqlite3_column_int(stmt, col);
					pStop->set_Edgeid(lFldVal);
				} 
				// j==19 Pos Along  Edge  
				col = 19;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_posalong(dFldVal);
				} 
				// j==20 arrDelay  
				col = 20;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_arrDelay(dFldVal); // this is omitted as Arrival Ride time data includes the deceleration delay
					//pStop->set_arrDelay(0);
				} 
				col = 21; // departure delay
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_depDelay(dFldVal); // omitted due to double counting as it is included in the segment time for RTD data
					//pStop->set_depDelay(0);
				} 
				if (sqlite3_column_count(stmt) > 22) {
				col = 22; // departure delay
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dFldVal= sqlite3_column_double(stmt, col);
						//pStop->set_depDelay(0);
					} 
				}
				if (sqlite3_column_count(stmt) > 23) {
					col = 23; // Historic indicator
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnHist(iFldVal); // historic stop?
					} 
				}
				if (sqlite3_column_count(stmt) > 24) {
					col = 24; // Eliminatable indicator
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnElim(iFldVal); // Eliminate in 1 yes, 0 no DP?
						if (iFldVal == 0 ) {elist.insert(pStop->get_id());}
					} 
				}
				if (pStop->get_id())
				{
					if(rows>1) { // compute the ride time using the ons and offs
						pStop->cdwlDelayh(gc->get_unitontm(),gc->get_unitofftm(),ph.get_hdway());
						dpsegRdTm = (dacArrTm - dpDepTm)*60 ; // run time in minutes between prev stop and current stop
						if (dpsegRdTm<=0) { // wrong data
							fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dpDepTm );
							logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Prev. Stop Dep. Tm. " << dpDepTm <<endl;
							//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",);
						}
					}
					if (dacArrTm > dacDepTm) { // wrong data
						fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Dep. Tm. " << dacDepTm <<endl;
						dacDepTm = dacArrTm + pStop->get_dwellDelay()/3600 ;
					}
					if (tp) {
						// compute the dwell time and historic run time using ons and offs
						//pStop->set_dwellDelay(max(pStop->get_dwellDelay(),(dacDepTm-dacArrTm)*60));
						pStop->cprobstoph(ph.get_hdway());
						if(rows>1) { // compute the ride time using the ons and offs
							pStop->set_CRdTm( dpCRdTm + dpsegRdTm+ dpDwellTm +(dpDepDelay/60)*dpprobStoph + (pStop->get_arrDelay()/60)*pStop->get_probStoph()); // drunTm*60);
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
							pStop->set_dwellDelay(0);
						}
					} else {
						if(rows>1) { // compute the ride time using the ons and offs
						pStop->set_CRdTm(drunTm*60); // actual arrival time from beg of route in mins 
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
						}
					}
					drunTimeC = pStop->get_CRdTm();
					if (rows==1) {
						pStop->set_dwellDelay(0);
						dpDwellTm=0;
						pStop->set_undCRdTm(pStop->get_CRdTm());
						dpArrDelay = 0;
						dpDepTm = dacDepTm;
					} else {
						pStop->set_undCRdTm(dpundRdTm+dpsegRdTm);  // add the segment ride time including the unaccounted delay at prev stop
						dDwell = pStop->get_dwellDelay()/60; // (max<double>((dacDepTm - dacArrTm),pStop->get_dwellDelay()/60)*60);  // no trips / hr = (60/ph.get_hdway()
						dpArrDelay = pStop->get_arrDelay();
					}
					dpDepVol += (pStop->get_HistOns() - pStop->get_HistOffs());
					pStop->set_HistDepVol(dpDepVol);
					m1.insert(sopair(pStop->get_id(),o1));
					dpDwellTm = pStop->get_dwellDelay()/60 ; // minutes
					dpDwell = (dacDepTm-dacArrTm)*60;
					dDwell = max<double>(0.0,(dpDwell-dpDwellTm));  // dwell delay that is not accounted by demand
					dpsegRdTm = (dacArrTm - dpDepTm)*60 + dDwell ;  // segement ride time
					dcumDwell += dpDwellTm + dDwell;
					dpDepDelay = pStop->get_depDelay();
					dpprobStoph = pStop->get_probStoph();
					dpArrTm = dacArrTm;
					dpDepTm = dacDepTm;
					dpundRdTm = pStop->get_undCRdTm();
					dpCRdTm = pStop->get_CRdTm();
				}
			} // loop over all rows 
			// update the last stop depDelay (0) and probStopping (1)  
			InpIt begin= m1.begin();
			InpIt end = m1.end();
			end--;
			end = m1.find(end->first);
			if(end!=m1.end()) { 
				o1 = end->second;
				o1.set_probStop(1.0);
				o1.set_depDelay(0.0);
				m1.erase(end); 
				m1.insert(sopair(o1.get_id(),o1));
			}
	//	} // if data count is greater than 0
		  // there are no rows to fetch
		sqlite3_finalize( stmt );
		return m1;
	}

	template <typename a, typename b,typename o,  typename c,typename k,typename g,typename p,typename t,typename r, typename d> 
	c& stopTableRTD(a& netdb ,b& sql,o& o1,c& m1,k& kstop,g& gc,p& ph,t& tp,r& srid, d& logFile)
	{
		int ret=0 , i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		sqlite3_stmt *stmt2 = NULL;
		int rows=0, columns=0 , col=0;
		char *zErrMsg = NULL;
		string sstopi="",strCDim = "xy", strGType="Point";
		bool blnRet = false;
		int pkCount = 0 , fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0, iGDim=2, iMeas=0, iSrid=srid;
		long lFldVal=0;
		double dFldVal=0,mcumDist=0,mpctDist=0,dschlTime=0, dDwell=0,dpDwell=0, dcumDwell=0, drunTimeC=0,dpundRdTm=0;
		double mRteLen=0,dpprobStoph=0,dpDwellTm=0,dpArrDelay=0,dpDepDelay=0,dpArrTm=0,dpDepTm=0,dpsegRdTm=0,dpCRdTm=0;
		// dwell - dwell delay from boarding and alighting cum (DepTime - ArrTime) minutes
		// dcumDwell - cumulative dwell delay from begining of route to the current stop minutes
		// dsegRunTm = dArrTm(current Stop) - dpDepTm(Prev Stop) - Dwell(Prev Stop)
		o* pStop;
		b tblStop="",tbLine="", tblBuf="", sqlQ="",tblBufx2="",stokens="", sstopid="", sFldVal="";
		typedef pair <int,o> sopair;
		typedef c::iterator InpIt; 
		InpIt stIt;

		vector<string> svecHMS, svecHAMPM;	
		vector<int> ivecHMS;	
		double dschlTm=0, dacArrTm =0,drunTm =0, dacDepTm=0, dschTmPt = 0, dbegArrTm=0, dbegDepTm=0,dpDepVol=0; 
		if (kstop.hOnSum() < kstop.hOffSum()) {
			dpDepVol = kstop.hOffSum() - kstop.hOnSum();
		}
		// drop the old table named without the direction indicator
		tblStop = "s" + kstop.route()+ kstop.schlName()+  to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
		ReplaceAll2(tblStop,"__","_");

		sqlQ="Drop Table " + tblStop + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		sqlQ="Drop Table " + tblStop +"_Buf" + "; "; 
		if ( sqlite3_exec( netdb, sqlQ.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblStop<< " dropped!"<<endl;
		}
		
		tblStop = "s" + kstop.route()+ kstop.schlName()+ kstop.dir()+ to_string<long>(kstop.tripId())+kstop.tripPeriod();
		replace(tblStop.begin(),tblStop.end(),' ','_');
		replace(tblStop.begin(),tblStop.end(),'(','_');
		replace(tblStop.begin(),tblStop.end(),')','_');
		replace(tblStop.begin(),tblStop.end(),':','_');
		replace(tblStop.begin(),tblStop.end(),'-','_');
		replace(tblStop.begin(),tblStop.end(),'/','x');
		ReplaceAll2(tblStop,"__","_");

		//replace(tblStop.begin(),tblStop.end(),'__','_');
		// Create a stop table for this run if doesn't exist already if it does drop and recreate it
		sqlQ = " as " + sql;
		blnRet = createSpaTbl (sqlQ,tblStop,netdb,srid,logFile);
		if ( blnRet )
		{  // recover geometry
			// first get teh geometry type and dimension details for the created table
			sqlQ = "SELECT Count(*) NumObj, GeometryType(\"Geometry\") GType, Srid(\"Geometry\") Srid, CoordDimension(\"Geometry\") CoorDim, ST_NDims(\"Geometry\") Dim , St_IsMeasured(\"Geometry\") Measured" 
				" from " + tblStop + " " ;
			if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
			{
					//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
					logFile<<"Table "<<tblStop<< " Geometric information could not be read !"<<endl;// assume generic
			} else { 

				// read tabel geoemtric information if it is available information
				rows=0;

				while ( sqlite3_step( stmt ) == SQLITE_ROW )
				{
					rows++;
					// query the stop geometry detail data from the stop table
					col=0; // count of geometries
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
					}
					col = 1;  // Geometric type "Point", "LineString", "PolyGon" , etc
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strGType = sFldVal;
					}
					col=2; // SRID 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iSrid = iFldVal;
					}
					col=3; // Coordinate Dimension as a string XY, XYM, XYZ, XYZM
					if (sqlite3_column_bytes(stmt, col) !=0) {
						sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						strCDim = sFldVal;
					}
					col=4; // Coordinate Dimension as an integer 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iGDim = iFldVal;
					}
					col=5; // If this is measured data 
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal = sqlite3_column_int(stmt, col);
						iMeas = iFldVal;
					}
				}
			}
			if (rows) {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
					"'" + strGType +"' , "  + to_string<int>(iGDim) + " );";
			} else  {
				sqlQ = " SELECT recovergeometrycolumn('" + tblStop + "', 'Geometry'," + (to_string<int>(iSrid)) + " ,'Point',2);";
			}
			blnRet = recoverSpatialGeometry(sqlQ,tblStop,netdb,srid,logFile);
			// create a line shape using the stop tbl just created to compute the cumulative distance to each stop 
			tbLine =  tblStop + "_RteGeom";
			sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
				" makeline(casttomultipoint( Geometry))  Geometry from " ;
			sqlQ.append(tblStop) ;
			sqlQ.append (" group by SchlName, RteName, DirName, TimePeriod, TripId;" );
			blnRet = createSpaTbl (sqlQ,tbLine,netdb,srid,logFile);

			if (blnRet) {// recover spatial geometry
				if (rows) {
					sqlQ = " SELECT recovergeometrycolumn('" +  tbLine  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
						" 'Linestring', "  + to_string<int>(iGDim) + " );";
				} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tbLine + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Linestring',2);";
				}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
			// create a buffer polygon using the stop tbl just created to select the demand distribution area 
				tblBuf =  tblStop + "_Buf";
				sqlQ = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer(( Geometry) , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tbLine + "  "
					" group by SchlName, RteName, DirName, TimePeriod, TripId;" ;
				blnRet = createSpaTbl (sqlQ,tblBuf,netdb,srid,logFile);
				if ( blnRet )
				{ // recover spatial geometry
					if (rows) {
						sqlQ = " SELECT recovergeometrycolumn('" +  tblBuf  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
							" 'Polygon', "  + to_string<int>(iGDim) + " );";
					} else  {
						sqlQ = " SELECT recovergeometrycolumn('" + tblBuf + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
					}
					blnRet = recoverSpatialGeometry(sqlQ,tblBuf,netdb,srid,logFile);
					// create a second buffer polygon using the polygon just created to select the edge distribution area 
					tblBufx2 = tblBuf + "x2" ;
					// Create a stop table for this run if doesn't exist already if it does drop an recreate it
					sqlQ  = " as select SchlName, RteName, DirName, TimePeriod, TripId, " 
					" st_buffer( Geometry , " + to_string<long>(gc->get_maxwalkdist()) + " ) " 
					" Geometry from "  + tblBuf + " " ;
					blnRet = createSpaTbl (sqlQ,tblBufx2,netdb,srid,logFile);
					if (blnRet) {
						if (rows) {
							sqlQ = " SELECT recovergeometrycolumn('" +  tblBufx2  + "', 'Geometry'," + (to_string<int>(iSrid)) + " , " 
								" 'Polygon', "  + to_string<int>(iGDim) + " );";
						} else  {
							sqlQ = " SELECT recovergeometrycolumn('" + tblBufx2 + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Polygon',2);";
						}
						blnRet = recoverSpatialGeometry(sqlQ,tblBufx2,netdb,srid,logFile);
					}
				} else
				{
				  // some error occurred
				  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				  logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
				}
			} // if make line succeeds then create the buffer polygons 
			else {
				  // some error occurred
				fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
				logFile<< "\nSQLite error: "<< sqlQ <<endl<<" DB Err. Msg "<<to_string<const char *>(sqlite3_errmsg( netdb)) ;
			}
		} else
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
		  //fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"\nSQLite error: "<<sqlQ <<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		}

		// query the stop table for this run and store it in the stop table

//		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  logFile<<"SQLite prepare error: \tSQL error: "<< sql<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
		  return m1;
		}
		m1.clear();

		// query the stop detail data from the stop table
		rows=0;
//		if(sqlite3_data_count(stmt) > 0 ) { 

			while ( sqlite3_step( stmt ) == SQLITE_ROW )
			{
				rows++;
				pStop = &o1;
				// query the stop detail data from the stop table
				col=0;
				iFldVal = sqlite3_column_int(stmt, col);
				if (iFldVal>0) {  // id
					pStop->set_id(rows);
					pStop->set_StOrdr(iFldVal);
				}
				// j==1 Schl Name
				// j==2 route name 
				// j==3 DirName
				// j==4 Time period 
				// j==5 set trip start time
				// j==6 Trip Key or Id
				col = 6;
				iFldVal = sqlite3_column_int(stmt, col);
				pStop->settripId(iFldVal);
				// j==7 Trip Number
				// j==8 scheduled start time time
				col = 8 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dschlTime = triphrs(svecHMS,iFldVal,dschlTime); 
					if (svecHAMPM.at(1) == "PM" && dschlTime < 12.0) {
						dschlTime += 12 ; 
					}
				}			//pStop->setschlTime(sFldVal);
				// j==9 Actual Arrival Time
				col = 9 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dacArrTm = triphrs(svecHMS,iFldVal,dacArrTm); 
					if(svecHAMPM.size()>1) { 
						if (svecHAMPM.at(1) == "PM" && dacArrTm < 12.0) {
							dacArrTm  += 12 ; 
						}
					}
					if (rows==1 )  {
						dbegArrTm = dacArrTm;
						dpArrTm = dacArrTm;
					}
						drunTm = (dacArrTm - dbegArrTm) ;  
				}			//pStop->setschlTime(sFldVal);
				// j==10 Actual Departure Time
				col = 10 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					svecHMS.clear();
					svecHAMPM.clear();
					Tokenize(sFldVal, svecHAMPM, " ");
					Tokenize(svecHAMPM.at(0), svecHMS, ":");
					dacDepTm = triphrs(svecHMS,iFldVal,dacDepTm); 
					if(svecHAMPM.size()>1) { 
						if (svecHAMPM.at(1) == "PM" && dacDepTm < 12.0) {
							dacDepTm += 12 ; 
						}
					}
					if (rows==1 )  {
						dbegDepTm = dacDepTm;
						dpDepTm = dacDepTm;
					}

				}			//pStop->setschlTime(sFldVal);
				// j==11 Schedule Run Time  
				col = 11 ;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					iFldVal= from_string<int> (sFldVal);
				} 
				// j==12 Stop Id
				col = 12;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopLbl(sFldVal);
					// calculate the cumulative distance from the begining of the route using the line geometry  
					//sqlQ = " select t1.Stop_Id , "
					//" st_length(st_splitLeft(makeline(casttomultipoint(t1.Geometry)),st_PointN(makeline(casttomultipoint(t1.Geometry)),t1.Id)))/5280 SplitLen, "
					//" ST_length(t2.Geometry) RteLength "
					sqlQ = " select t1.Stop_Id , ST_Line_Locate_Point(t2.Geometry, t1.Geometry) PctDist, st_length(t2.Geometry)/5280 RteLength "
						", st_length(t2.Geometry) * ST_Line_Locate_Point(t2.Geometry, t1.Geometry) CumDist "  
							" from " + tblStop + " t1 , " + tbLine + " t2  where t1.Stop_Id like '" + sFldVal + "' ; " ;
					// query the data 
					if ( sqlite3_prepare_v2( netdb, sqlQ.c_str(), -1, &stmt2, NULL ) != SQLITE_OK )
					{
					  // some error occurred
					  fprintf(stderr, "\nStop Distance Query error: %s\n\nSQL: %s \n", sqlQ.c_str() ,  sqlite3_errmsg( netdb) );
					  logFile<<"Stop Distance Query: \tSQL error: "<< sqlQ<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( netdb))<<endl;
					  return m1;
					}
					// j==1 PctDist 
					while ( sqlite3_step( stmt2 ) == SQLITE_ROW )
					{
						col = 1;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
						mpctDist = sqlite3_column_double(stmt2, col);
						}
						// j==2 Rte Length  
						col = 2;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mRteLen = sqlite3_column_double(stmt2, col);
						}  
						// j==3 Cum Length  
						col = 3;
						if (sqlite3_column_bytes(stmt2, col) !=0) {
							mcumDist = sqlite3_column_double(stmt2, col);
						}  
					}  
					mcumDist = mRteLen * mpctDist;
					
						//st_PointN(makeline(casttomultipoint(t1.Geometry)),3)) Point3,
						//st_length(makeline(casttomultipoint(t1.Geometry))) TotalLeng
					//	pStop->set_CumDist(dacDepTm);
						pStop->set_CumDist(mcumDist);

				}
				// j==13 Main and Cross Street Intersection 
				col = 13;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
					pStop->set_StopName(sFldVal);
				}
				// j==14 Trip Boardings  
				col = 14;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOns(dFldVal);
					pStop->set_Ons(dFldVal);
				}  
				// j==15 Trip alightings 
				col = 15;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_HistOffs(dFldVal);
					pStop->set_Offs(dFldVal);
				}
				// j==16 location latitude , xc
				col = 16;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_xc(dFldVal);
				}   
				// j==17 location longitude , yc 
				col = 17;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_yc(dFldVal);
				} 
				// j==18 EdgeID 
				col = 18;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					lFldVal= sqlite3_column_int(stmt, col);
					pStop->set_Edgeid(lFldVal);
				} 
				// j==19 Pos Along  Edge  
				col = 19;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_posalong(dFldVal);
				} 
				// j==20 arrDelay  
				col = 20;
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_arrDelay(dFldVal); // this is omitted as Arrival Ride time data includes the deceleration delay
					//pStop->set_arrDelay(0);
				} 
				col = 21; // departure delay
				if (sqlite3_column_bytes(stmt, col) !=0) {
					dFldVal= sqlite3_column_double(stmt, col);
					pStop->set_depDelay(dFldVal); // omitted due to double counting as it is included in the segment time for RTD data
					//pStop->set_depDelay(0);
				} 
				if (sqlite3_column_count(stmt) > 22) {
				col = 22; // departure delay
					if (sqlite3_column_bytes(stmt, col) !=0) {
						dFldVal= sqlite3_column_double(stmt, col);
						//pStop->set_depDelay(0);
					} 
				}
				if (sqlite3_column_count(stmt) > 23) {
					col = 23; // Historic indicator
					if (sqlite3_column_bytes(stmt, col) !=0) {
						iFldVal= sqlite3_column_int(stmt, col);
						pStop->set_blnHist(iFldVal); // historic stop?
					} 
				}
				if (pStop->get_id())
				{
					if(rows>1) { // compute the ride time using the ons and offs
						pStop->cdwlDelayh(gc->get_unitontm(),gc->get_unitofftm(),ph.get_hdway());
						dpsegRdTm = (dacArrTm - dpDepTm)*60 ; // run time in minutes between prev stop and current stop
						if (dpsegRdTm<=0) { // wrong data
							fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dpDepTm );
							logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Prev. Stop Dep. Tm. " << dpDepTm <<endl;
							//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f < Previous stop Dep. Time %8.5f\n",);
						}
					}
					if (dacArrTm > dacDepTm) { // wrong data
						fprintf(stderr, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						//fprintf(logFile, "\nStop %d name %s Data error Arr. Time: %8.5f > Dep. Time %8.5f\n",pStop->get_id(),pStop->get_StopName().c_str(),dacArrTm ,dacDepTm );
						logFile<<"Stop "<<pStop->get_id()<<", \t  " <<pStop->get_StopName()<< ", \t Data Error Arrival Tm " << dacArrTm<< ", \t Dep. Tm. " << dacDepTm <<endl;
						dacDepTm = dacArrTm + pStop->get_dwellDelay()/3600 ;
					}
					if (tp) {
						// compute the dwell time and historic run time using ons and offs
						//pStop->set_dwellDelay(max(pStop->get_dwellDelay(),(dacDepTm-dacArrTm)*60));
						pStop->cprobstoph(ph.get_hdway());
						if(rows>1) { // compute the ride time using the ons and offs
							pStop->set_CRdTm( dpCRdTm + dpsegRdTm+ dpDwellTm +(dpDepDelay/60)*dpprobStoph + (pStop->get_arrDelay()/60)*pStop->get_probStoph()); // drunTm*60);
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
							pStop->set_dwellDelay(0);
						}
					} else {
						if(rows>1) { // compute the ride time using the ons and offs
						pStop->set_CRdTm(drunTm*60); // actual arrival time from beg of route in mins 
						} else {
							pStop->set_CRdTm(0); // actual arrival time from beg of route in mins 
							pStop->set_probStoph(1.0);
						}
					}
					drunTimeC = pStop->get_CRdTm();
					if (rows==1) {
						pStop->set_dwellDelay(0);
						dpDwellTm=0;
						pStop->set_undCRdTm(pStop->get_CRdTm());
						dpArrDelay = 0;
						dpDepTm = dacDepTm;
					} else {
						pStop->set_undCRdTm(dpundRdTm+dpsegRdTm);  // add the segment ride time including the unaccounted delay at prev stop
						dDwell = pStop->get_dwellDelay()/60; // (max<double>((dacDepTm - dacArrTm),pStop->get_dwellDelay()/60)*60);  // no trips / hr = (60/ph.get_hdway()
						dpArrDelay = pStop->get_arrDelay();
					}
					dpDepVol += (pStop->get_HistOns() - pStop->get_HistOffs());
					pStop->set_HistDepVol(dpDepVol);
					m1.insert(sopair(pStop->get_id(),o1));
					dpDwellTm = pStop->get_dwellDelay()/60 ; // minutes
					dpDwell = (dacDepTm-dacArrTm)*60;
					dDwell = max<double>(0.0,(dpDwell-dpDwellTm));  // dwell delay that is not accounted by demand
					dpsegRdTm = (dacArrTm - dpDepTm)*60 + dDwell ;  // segement ride time
					dcumDwell += dpDwellTm + dDwell;
					dpDepDelay = pStop->get_depDelay();
					dpprobStoph = pStop->get_probStoph();
					dpArrTm = dacArrTm;
					dpDepTm = dacDepTm;
					dpundRdTm = pStop->get_undCRdTm();
					dpCRdTm = pStop->get_CRdTm();
				}
			} // loop over all rows 
			// update the last stop depDelay (0) and probStopping (1)  
			InpIt begin= m1.begin();
			InpIt end = m1.end();
			end--;
			end = m1.find(end->first);
			if(end!=m1.end()) { 
				o1 = end->second;
				o1.set_probStop(1.0);
				o1.set_depDelay(0.0);
				m1.erase(end); 
				m1.insert(sopair(o1.get_id(),o1));
			}
	//	} // if data count is greater than 0
		  // there are no rows to fetch
		sqlite3_finalize( stmt );
		return m1;
	}






	template <typename a, typename b,typename o,  typename c, typename d> 
	c& stopTripTableData(a& netdb ,b& sql,o& o1, c& mapData,d& logFile)
	{
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows;
		int columns;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		long lFldVal=0;
		o* pStop;
		typedef pair <int,o> sopair;
		string svertid="", sstopid="", sFldVal="";
		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( ret != SQLITE_OK )
			return mapData;
		if ( rows < 1 )
			;
		else
		{
			// query the stop detail data from the stop table
			mapData.clear();
			fprintf(logFile," Sql Query %s rows %d cols %d !",sql.c_str(),rows,columns);
			for ( i = 1; i <= rows; i++ )
			{
				pStop = &o1;
				for (j = 0; j <= columns-1; j++ )
				{
					if (j==0) // id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						pStop->set_id(iFldVal);

					} else if (j==1) // Schl Name
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setAgencyId(sFldVal);
					} else if (j==2) // route name 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setRouteName(sFldVal);
					} else if (j==3) // DirName
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setdirName(sFldVal);
					} else if (j==4) // Time period 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->settimePeriod(sFldVal);
					} else if (j==5) // set trip start time
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//dFldVal= from_string<double> (sFldVal);
						pStop->settripSTime(sFldVal);
					} else if (j==6) // Trip Id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->settripId(sFldVal);
					} else if (j==7) // Trip Number
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						pStop->settripNumber(iFldVal);
					} else if (j==8) // schedule time
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setschlTime(sFldVal);
					} else if (j==9) // Actual Arrival Time
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setacTimeArr(sFldVal);
					} else if (j==10) // Actual Departure Time
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setacTimeDep(sFldVal);
					} else if (j==11) // Schedule Run Time  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						pStop->setschlRunTm(iFldVal);
					} else if (j==12) // Stop Id  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						lFldVal= from_string<long> (sFldVal);
						pStop->setstopId(lFldVal);
					}  else if (j==13) // Parent Station 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						pStop->setstrMainCrs(sFldVal);
					}  else if (j==14) // Trip Boardings  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setons(dFldVal);
					}  else if (j==15) // Trip alightings 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setoffs(dFldVal);
					}  else if (j==16) // location latitude , yc
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setyc(dFldVal);
					}   else if (j==17) // location longitude , xc 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						pStop->setxc(dFldVal);
					} 

				} // loop over cols
				if (pStop->stOrdr())
				{
					mapData.insert(sopair(pStop->stOrdr(),o1));
				}
			} // loop over rows
			sqlite3_free_table( results );
		} // if rows are found
		return mapData;
	}


	template <typename a, typename b, typename c, typename d> 
	c& stopTripData(a& netdb ,b& sql, c& mapData,d& logFile)
	{
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows;
		int columns;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		string svertid="", sstopid="", sFldVal="";
		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( ret != SQLITE_OK )
			return mapData;
		if ( rows < 1 )
			;
		else
		{
			// query the stop detail data from the stop table
			mapData.clear();
			fprintf(logFile," Sql Query %s rows %d cols %d !",sql.c_str(),rows,columns);
			for ( i = 1; i <= rows; i++ )
			{
				sDmnd mstop;
				for (j = 0; j <= columns-1; j++ )
				{
					if (j==0) // ADay
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						mstop.setaDay(iFldVal);

					} else if (j==1) // Route
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						mstop.setRoute(sFldVal);
					} else if (j==2) // Direction  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						mstop.setDirn(iFldVal);
					} else if (j==3) // stop sequence
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						mstop.set_StOrdr(iFldVal);
					} else if (j==4) // Agency Name
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						mstop.setAgencyId(sFldVal);
					} else if (j==5) // Trip
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						mstop.setTrip(sFldVal);
					} else if (j==6) // Trip Block Id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						mstop.setTripBlkId(sFldVal);
					} else if (j==7) // Stop ID/QStop(demand Table) from Stop table
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						mstop.setid(iFldVal);
					} else if (j==8) // Num of Samples
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						mstop.setnSamples(iFldVal);
					} else if (j==9) // Number of Ons
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						mstop.setOns(dFldVal);
					} else if (j==10) // Number of Ons
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						mstop.setOffs(dFldVal);
					} else if (j==11) // Stop Latitude 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						mstop.setStopLat(sFldVal);
					} else if (j==12) // Stop Longitude 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						mstop.setStopLong(sFldVal);
					}
				} // loop over cols
				if (mstop.stOrdr())
				{
					string key = to_string <int> (mstop.stOrdr()) + 'x' +  (mstop.route()) + 'x'
						+ (mstop.trip()) + 'x' + to_string <short> (mstop.dirn()) ;
					mapData.insert(strDmndpair(key,mstop));
				}
			} // loop over rows
			sqlite3_free_table( results );
		} // if rows are found
		return mapData;
	}


	template <typename a, typename b,typename c,  typename d, typename e> 
	d& stopTableCount(a& netdb ,b& sql, c& c1, d& mkTrip, e& logFile)
	{

		int ret=0 , i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows=0, columns=0 , col=0;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0 , fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		long lFldVal=0;
		double dFldVal=0,mcost=0,dschlTime=0, dDwell=0, dcumDwell=0, dblHdway=0, dblnextschlTime=0;
		// dwell - dwell delay from boarding and alighting cum (DepTime - ArrTime) minutes
		// dcumDwell - cumulative dwell delay from begining of route to the current stop minutes
		c* pkTrip;
		typedef pair <int,c> sopair;
		string stokens="", sstopid="", sFldVal="";
		vector<string> svecHMS, svecHAMPM;	
		vector<int> ivecHMS;	
		double dschlTm=0; 
//		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  return mkTrip;
		}
		mkTrip.clear();
		// query the stop count summary data into stop key set 
		//	c  kTrip;

		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			rows++;
			pkTrip = &c1;
			// query the trip table summary data using the Key Data
			// (j==0) // Route Name
			col=0;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pkTrip->setRoute(sFldVal);
			}
			// j==1 Schl Name
			col=1;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pkTrip->setschlName(sFldVal);
			}
			// (j==2) // Trip Period
			col=2;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						pkTrip->settripPeriod(sFldVal);
			} 
			// (j==3) // Trip Start Time
			col=3;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pkTrip->settripSTime(sFldVal);
				svecHMS.clear();
				svecHAMPM.clear();
				Tokenize(sFldVal, svecHAMPM, " ");
				Tokenize(svecHAMPM.at(0), svecHMS, ":");
				dschlTime = triphrs(svecHMS,iFldVal,dschlTime); 
				if (svecHAMPM.size()>1) {
					if (svecHAMPM.at(1) == "PM") {
						dschlTime += 12 ; 
					}
				}
				pkTrip->settripTime(dschlTime);
				if (rows==1) { // greatest trip number descending for headway calc using schedule time (time till next trip) 
					dblnextschlTime = dschlTime;
				} else {
					dblHdway = (dblnextschlTime - dschlTime)*60;
					pkTrip->settimeNexTrip(dblHdway);
					dblnextschlTime = dschlTime;
				}
			}		

				// else if (j==4) // Trip Id / Key  
			col=4;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal = (sqlite3_column_int(stmt, col));
				pkTrip->settripId(lFldVal);
			} 
			// (j==5) // TripNumber
			col=5;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal = (sqlite3_column_int(stmt, col));
				pkTrip->settripNumber(lFldVal);
					}
			// (j==6) // Trip stop count 
			col=6;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal = (sqlite3_column_int(stmt, col));
				pkTrip->setstopCount(lFldVal);
			}
			// (j==7) // Trip sum of Boardings 
			col=7;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pkTrip->sethOnSum(dFldVal);
			}
			// (j==8) // Trip sum of Alightings 
			col=8;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pkTrip->sethOffSum(dFldVal);
			}
			// (j==9) // Route Direction 
			col=9;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pkTrip->setdir(sFldVal);
			}
			if (pkTrip->stopCount())
			{
				mkTrip.insert(tskidpair(pkTrip->tripNumber(),*pkTrip));
			}
		} // loop over rows
			sqlite3_finalize( stmt );
		return mkTrip;
	}

template <typename a, typename q, typename k,typename o,typename m,  typename d, typename e >
m& readLUCodesTable(a& netdb, q& q1, k& k1,o& o1, m& m1, d& logFile,e& tmPd)
{
	o* op1;
	op1=&o1;
    typedef pair <k,o> qopair;
    string rec1;
//	char *seps = "\t";
	long ip=0;
	double ecost=0;
	int ret;
	int i=0, j=0;
	sqlite3_stmt *stmt = NULL;
	char **results;
	int rows=0;
	int columns=0 , col=0;
	char *zErrMsg = NULL;
	string pkName="", sstopi="", sFldVal="";
	int pkCount = 0;
	int fldNo = 0, isOveRide=0;
	int iFldVal=0, stopid=0,stopi=0;
	long lFldVal=0;
	double dFldVal=0,mcost=0,dschlTime=0;
	size_t isAM = 0, isPM = 0;

	isAM = tmPd.find("AM");
	isPM = tmPd.find("PM");
	
	//ret = sqlite3_get_table( netdb, q1.c_str(), &results, &rows, &columns, &zErrMsg );
	if ( sqlite3_prepare_v2( netdb, q1.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
	{
	  // some error occurred
	  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
	  fprintf(logFile,"\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );  //"SQLite error:  "<<q1<< " !"<<endl<<to_string<const char *>(sqlite3_errmsg( netdb)) <<endl;
	  return m1;
	}
	m1.clear();
	// query the stop detail data from the stop table
	
	while ( sqlite3_step( stmt ) == SQLITE_ROW )
	{
		rows++;
		op1 = &o1;
		// query the stop detail data from the stop table
		// (j==0) // Ptype
		col=0;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						op1->set_pType(sFldVal);
		}
		// (j==1) // vx from 
		col=1; // description
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						//o1.set_desc(sFldVal);
		}
		// (j==2) // Vx to
		col=2; // Land Use Code
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						op1->set_LUC(sFldVal);
		}
		// (j==3) // Property Type Key text Code
		col=3;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
						op1->set_KeyProp(sFldVal);
		}

		col=6;
		if (sqlite3_column_bytes(stmt, col) !=0) {
			isOveRide = sqlite3_column_double(stmt, col);
			op1->set_PdOveRide(isOveRide);
		}
		// (j==4) // Unit On rate
		col=4;

		dFldVal = sqlite3_column_double(stmt, col);
		if (dFldVal>0  ) {  // id
			if (isAM==0 && op1->get_LUC()=="COM" && isOveRide == 0 ) {
				op1->set_OnCoeff(0);
			} else if (isPM==0 && op1->get_LUC()=="RES" && isOveRide == 0) {
				op1->set_OnCoeff(0);
			} else {
				op1->set_OnCoeff(dFldVal);
			}

		} 
		// (j==5) // Unit Off rate
		col=5;

		dFldVal = sqlite3_column_double(stmt, col);
		if (dFldVal>0) {  // Off Rate
			if (isPM==0 && op1->get_LUC()=="COM" && isOveRide == 0) {
				op1->set_OffCoeff(0);
			} else if (isAM==0 && op1->get_LUC()=="RES" && isOveRide == 0) {
				op1->set_OffCoeff(0);
			} else {
				op1->set_OffCoeff(dFldVal);
			}
		
		}


		ip++;
		if (op1->get_pType().length())
		{
			m1.insert(lucpair(o1.get_pType(),o1));
			o1.show_lucodes(cout);
		}
	} // loop over rows
	sqlite3_finalize(stmt );

	return  m1;

}



	template <typename a, typename b,typename o,  typename c, typename d> 
	c& luCodes(a& netdb ,b& sql,o& o0, c& mapData,d& logFile)
	{
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows;
		int columns;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		string svertid="", sstopid="", sFldVal="";
		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( ret != SQLITE_OK )
			return mapData;
		if ( rows < 1 )
			;
		else
		{
			// query the stop detail data from the stop table
			mapData.clear();
			fprintf(logFile," Sql Query %s rows %d cols %d !",sql.c_str(),rows,columns);
			for ( i = 1; i <= rows; i++ )
			{
				o o1;
				for (j = 0; j <= columns-1; j++ )
				{
					if (j==0) // PK_UID
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						//o1.set_id(iFldVal);

					} else if (j==1) // PType
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						o1.set_pType(sFldVal);
					} else if (j==2) // Description  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//o1.set_Desc(sFldVal);
					} else if (j==3) // Land Use Code
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//iFldVal= from_string<int> (sFldVal);
						o1.set_LUC(sFldVal);
					} else if (j==4) // Property Type Text Code
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						o1.set_KeyProp(sFldVal);
					} else if (j==5) // Unit On rate per key field
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						o1.set_OnCoeff(dFldVal);
					} else if (j==6) // Unit Off rate per key field
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						o1.set_OffCoeff(dFldVal);
					} 
				} // loop over cols
					mapData.insert(lucpair(o1.get_pType(),o1));
					o1.show_lucodes(cout);
			} // loop over rows
			sqlite3_free_table( results );
		} // if rows are found
		return mapData;
	}


	template <typename a, typename q,typename o,  typename c, typename d> 
	c& pdHdwy(a& netdb , q& q1,o& o0, c& m1,d& logFile)
	{
//FLDPERD,FLDBEGT,FLDHDWAY,FLDPDLEN,COSTOPER
//SELECT FLDPERD, FLDBEGT, FLDHDWAY, FLDPDLEN, COSTOPER FROM pdhdwayIB ORDER BY FLDPERD, FLDBEGT
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows=0;
		int columns;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0, col=0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		string svertid="", sstopid="", sFldVal="";
		if ( sqlite3_prepare_v2( netdb, q1.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(logFile,"\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );  //"SQLite error:  "<<q1<< " !"<<endl<<to_string<const char *>(sqlite3_errmsg( netdb)) <<endl;
		  return m1;
		}
		m1.clear();
		// query the stop detail data from the stop table

		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			rows++;
			// query the period detail data from the period table
			o o1;
			// (j==0) // FLDPERD
			col=0;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal = sqlite3_column_int(stmt, col);
				o1.set_pdId(iFldVal);
			}
			// if (j==1) // FLDBEGT
			col=1;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				o1.set_begtm(dFldVal);
			}
			// (j==2) // FLDHDWAY  
			col=2;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				o1.set_hdway(dFldVal);
			} 
			// (j==3) // FLDPDLEN
			col=3;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				o1.set_pdlen(dFldVal);
			} 
			// if (j==4) // COSTOPER
			col=4;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				o1.set_opercost(dFldVal);
			} 
			// if (j==5) // period key
			col=5;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				o1.set_pdKey(sFldVal);
			} 
			// if (j==6) // number of trips
			col=6;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal = sqlite3_column_double(stmt, col);
				o1.set_numTrips(iFldVal);
			} 
			// if (j==7) // Include
			col=7;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal = sqlite3_column_double(stmt, col);
				o1.set_include(iFldVal);
			} 
				m1.insert(phwpair(o1.get_pdId(),o1));
		} // loop over rows
		sqlite3_finalize(stmt );
		return m1;
	}

	template <typename a, typename b,typename o,  typename c, typename d> 
	c& globalCost(a& netdb ,b& sql,o& o1, c& m1,d& logFile)
	{
//SELECT  COSTWALK, COSTRIDE, UNITONTM, UNITOFFTM, MAXWLKDIST, PROPENSITY, FILESTEM, NOPERIODS, 
//WALKSPD, FILEPATH, MAXSKIP, DPDIMENSION , unitConv FROM GlobCostPXIBD5
		o* op1;
		op1=&o1;
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		int dp=5, col=0;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		long lFldVal=0;
		typedef pair <long,o> lopair;
		string svertid="", sstopid="", sFldVal="";
		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  return m1;
		}

		// query the stop detail data from the stop table
		m1.clear();
		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			//if (j==0) // COSTWALK
			col=0;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_walkcost(dFldVal);
			} 
			col=1; // COSTRIDE
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_ridecost(dFldVal);
			} 
			col=2; // UNITONTM
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_unitontm(dFldVal);
			} 
			col=3; // UNITOFFTM  
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_unitofftm(dFldVal);
			} 
			col=4; // MAXWLKDIST 
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_maxwalkdist(dFldVal);
			} 
			col=5; // PROPENSITY
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_propensity(dFldVal);
			} 
			col=6; //FILESTEM 
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				op1->set_filestem(sFldVal);
			} 
			col=7; // NOPERIODS
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal = sqlite3_column_int(stmt, col);
				op1->set_nopds(iFldVal);
			} 
			col=8; // WALKSPD
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_walkspd(dFldVal);
			} 
			col=9; // FILEPATH
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				op1->set_filepath(sFldVal);
			} 
			col=10; // MAXSKIP
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal = sqlite3_column_int(stmt, col);
				op1->set_maxskip(iFldVal);
			} 
			col=11; // DPDIMENSION
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal = sqlite3_column_int(stmt, col);
				op1->set_dpdimension(iFldVal);
			} 
			col=12; // unitConv
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = sqlite3_column_double(stmt, col);
				op1->set_unitconv(dFldVal);
			} 
			if (op1->get_dpdimension()>0) {
				m1.insert(lopair(op1->get_dpdimension(),*op1));
				} else {
					// assume default value from dp variable
					m1.insert(lopair(dp,*op1));
				}
			} // loop over rows
		sqlite3_finalize( stmt );
		return m1;
	}


// Split Strings using a delimiter and put them in an argument vector string

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// Split Strings using a delimiter and return a vector string

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


template <typename a, typename b, typename c> 
c&  triphrs(a& sVectHMS, b& b1, c& tm)
{
	c hr = from_string<b>(sVectHMS.at(0));
	c min = from_string<b>(sVectHMS.at(1));
	c sec = from_string<b>(sVectHMS.at(2));
	tm = hr + min/60 + sec/3600 ;
    return tm;
}


	template <typename a, typename b, typename c, typename d, typename e> 
	d& parcelMapTable(a& netdb ,b& sql, c& c1, d& mapData,e& logFile)
	{
		c* op1;
		op1=&c1;
		typedef pair <b,c> qopair;
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows;
		int columns;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		long lFldVal=0;
		string sTract="", sBlock="", sFldVal="";
//		StrSQLstoptable = " SELECT  PK_UID,GEOID10, NAME10, STATE, COUNTY, PLACE, TRACT, BLOCK, "
//							" TRACTBLOCK, POP, HU, PLACENAME, EdgeID, Geometry "
//							" FROM RTD_Blocks_CD "
//							" ORDER BY PK_UID; ";

		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &zErrMsg );
		if ( ret != SQLITE_OK )
			return mapData;
		if ( rows < 1 )
			;
		else
		{
			// query the stop detail data from the stop table
			mapData.clear();
			fprintf(logFile," Sql Query %s rows %d cols %d !",sql.c_str(),rows,columns);
			for ( i = 1; i <= rows; i++ )
			{
				c blk;
				for (j = 0; j <= columns-1; j++ )
				{
					if (j==0) // pk uid
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						blk.set_id(iFldVal);

					} else if (j==1) // Block Id
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						blk.set_pacid(sFldVal);
					} else if (j==2) // Block Name  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//iFldVal= from_string<int> (sFldVal);
						//blk.setDirn(iFldVal);
					} else if (j==3) // state
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//iFldVal= from_string<int> (sFldVal);
						//blk.setStOrdr(iFldVal);
					} else if (j==4) // county 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//blk.setAgencyId(sFldVal);
					} else if (j==5) // Place
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						//blk.setTrip(sFldVal);
					} else if (j==6) // Tract 
					{
						sTract = to_string<char * >(results[( i * columns ) + j]);
						//blk.setTripBlkId(sFldVal);
					} else if (j==7) // block number
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						lFldVal= from_string<long> ( sTract + sFldVal);
						blk.set_id(lFldVal);
					} else if (j==8) // Pop 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						iFldVal= from_string<int> (sFldVal);
						blk.set_pval(iFldVal);
					} else if (j==9) // Number of Household Units
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						blk.set_grarea(dFldVal);
					} else if (j==10) // Edgeid
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						lFldVal= from_string<long> (sFldVal);
						blk.set_eoid(lFldVal);
						dFldVal= 0.5; 
						blk.set_palong(dFldVal);
					} else if (j==11) // property code 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						blk.set_ptype(sFldVal);
					} else if (j==12) //Geometry get x coordinates 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						blk.set_xc(dFldVal);
					} else if (j==13) //Geometry get y Coordinates  
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						blk.set_yc(dFldVal);
					} else if (j==14) //Geometric area 
					{
						sFldVal = to_string<char * >(results[( i * columns ) + j]);
						dFldVal= from_string<double> (sFldVal);
						blk.set_yc(dFldVal);
					}
				} // loop over cols
				if (blk.get_id())
				{
					mapData.insert(PE_Pair(blk.get_id(),blk));
				}
			} // loop over rows
			sqlite3_free_table( results );
		} // if rows are found
		return mapData;
	}

	template <typename a, typename b, typename c, typename d, typename e> 
	d& parcelMapTableState(a& netdb ,b& sql, c& c1, d& m1,e& logFile)
	{
		c* pBlk;
		typedef pair <long,c> parpair;
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows=0;
		int columns=0, col=0;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		long lFldVal=0;
		string sTract="", sBlock="", sFldVal="";
				//strSQLTable = " SELECT  PK_UID,GEOID10, NAME10, STATE, COUNTY, PLACE, TRACT, BLOCK, "
				//	" POP, HU, EdgeID,Ptype, trim(Round(st_x(st_centroid(Geometry)),3)) xC, "
				//	" trim(Round(st_y(st_centroid(Geometry)),3)) yC, " 
				//	" trim(round(st_area(Geometry),3)) Area, PLACENAME "
				//				" FROM RTD_Blocks_CD  ORDER BY PK_UID; ";

		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stdout, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  return m1;
		}
		m1.clear();
		// query the Parcel detail data from the stop table

		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			pBlk = &c1;
			col=0; // Parcel/Block Id
			iFldVal = sqlite3_column_int(stmt, col);
			if (iFldVal>0) {  // 
				pBlk->set_id(iFldVal);
			}
			// j==1 // Parcel/Block Id
			col = 1 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pBlk->set_pacid(sFldVal);
			} 
			//  j==2 // Parcel/block number
			col = 2 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal= sqlite3_column_int(stmt, col);
				//pBlk->set_blockid(lFldVal);
			} 
			// 	j==3 // Pop 
			col = 3 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal= sqlite3_column_double(stmt, col);
				pBlk->set_pval(dFldVal);
			} 
			// j==4 // Edge Id 
			col = 4 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal= sqlite3_column_int(stmt, col);
				pBlk->set_eoid(lFldVal);
			} 
			// j==5 // property code 
			col = 5 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pBlk->set_ptype(sFldVal);
			} 
			// j==6 // Land use code 
			col = 6 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pBlk->set_lucd(sFldVal);
			} 
			// 	j==7  //Geometry get x coordinates 
			col = 7 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_xc(dFldVal);
			} 
			// 	j==8 //Geometry get y Coordinates  
			col = 8 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_yc(dFldVal);
			} 
			//	j==9 //Geometric area 
			col = 9 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_grarea(dFldVal);
			} 
			col = 10 ; // position along
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_palong(dFldVal);
			} 
			col = 11 ; // nearest dist to network 
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				//pBlk->set_palong(dFldVal);
			} 
			// 	j==12 // PopOff 
			col = 12 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal= sqlite3_column_double(stmt, col);
				pBlk->set_aval(dFldVal);
			} 
			if (pBlk->get_id())
			{
				m1.insert(PE_Pair(pBlk->get_id(),*pBlk));
			}
		} // loop over all rows
		  // there are no rows to fetch
			sqlite3_finalize( stmt );
		return m1;
	}


	template <typename a, typename b, typename c, typename d, typename e> 
	d& blockMapTableState(a& netdb ,b& sql, c& c1, d& m1,e& logFile)
	{
		c* pBlk;
		typedef pair <long,c> parpair;
		int ret;
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		char **results;
		int rows=0;
		int columns=0, col=0;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0;
		int iFldVal=0, stopid=0,stopi=0;
		double dFldVal=0,mcost=0;
		long lFldVal=0;
		string sTract="", sBlock="", sFldVal="";
				//strSQLTable = " SELECT  PK_UID,GEOID10, NAME10, STATE, COUNTY, PLACE, TRACT, BLOCK, "
				//	" POP, HU, EdgeID,Ptype, trim(Round(st_x(st_centroid(Geometry)),3)) xC, "
				//	" trim(Round(st_y(st_centroid(Geometry)),3)) yC, " 
				//	" trim(round(st_area(Geometry),3)) Area, PLACENAME "
				//				" FROM RTD_Blocks_CD  ORDER BY PK_UID; ";

		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stdout, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(logFile, "\nSQLite error: %s\n\nSQL: %s \n", sql.c_str() ,  sqlite3_errmsg( netdb) );
		  return m1;
		}
		m1.clear();
		// query the Parcel detail data from the stop table

		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			pBlk = &c1;
			col=0;
			iFldVal = sqlite3_column_int(stmt, col);
			if (iFldVal>0) {  // pk uid
				pBlk->set_id(iFldVal);
			}
			// j==1 // Parcel/Block Id
			col = 1 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pBlk->set_pacid(sFldVal);
			} 
			//  j==2  // state
			col = 2 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
			} 
			// j==3  // county 
			col = 3 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
			} 
			// j==4 // Place
			col = 4 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				//pBlk->set_placeid(sFldVal);
			} 
			//  j==5  // Tract 
			col = 5 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				//pBlk->set_tractid(sFldVal);
			} 
			//  j==6 // block number
			col = 6 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal= sqlite3_column_int(stmt, col);
				//pBlk->set_blockid(lFldVal);
			} 
			// 	j==7 // Pop 
			col = 7 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				iFldVal= sqlite3_column_int(stmt, col);
				pBlk->set_pval(iFldVal);
			} 
			// j==8 // Number of Household Units
			col = 8 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				//pBlk->set_grarea(dFldVal);
			} 
			// (j==9) // Edgeid
			col = 9;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal= sqlite3_column_int(stmt, col);
				pBlk->set_eoid(lFldVal);
			} 
			// j==10 // property code 
			col = 10 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pBlk->set_ptype(sFldVal);
			} 
			// j==11 // Land use code 
			col = 11 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				sFldVal = to_string<const unsigned char * > (sqlite3_column_text(stmt, col));
				pBlk->set_lucd(sFldVal);
			} 
			// 	j==12  //Geometry get x coordinates 
			col = 12 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_xc(dFldVal);
			} 
			// 	j==13 //Geometry get y Coordinates  
			col = 13 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_yc(dFldVal);
			} 
			//	j==14 //Geometric area 
			col = 14 ;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_grarea(dFldVal);
			} 
			col = 15 ; // position along
			if (sqlite3_column_bytes(stmt, col) !=0) {
				dFldVal = (sqlite3_column_double(stmt, col));
				pBlk->set_palong(dFldVal);
			} 
			if (pBlk->get_id())
			{
				m1.insert(PE_Pair(pBlk->get_id(),*pBlk));
			}
		} // loop over all rows
		  // there are no rows to fetch
			sqlite3_finalize( stmt );
		return m1;
	}



void initializeSpatialMetadata( sqlite3 *sqlite_handle , FILE* logFile)
{
// attempting to perform self-initialization for a newly created DB
  int ret;
  char sql[1024];
  char *zErrMsg = NULL;
  int count = 0;
  int i;
  char **results;
  int rows;
  int columns;

  if ( sqlite_handle == NULL )
    return;
  // checking if this DB is really empty
  strcpy( sql, "SELECT Count(*) from sqlite_master" );
  ret = sqlite3_get_table( sqlite_handle, sql, &results, &rows, &columns, NULL );
  if ( ret != SQLITE_OK )
    return;
  if ( rows < 1 )
    ;
  else
  {
    for ( i = 1; i <= rows; i++ )
      count = atoi( results[( i * columns ) + 0] );
  }
  sqlite3_free_table( results );

  if ( count > 0 )
    return;

  // all right, it's empty: proceding to initialize
  strcpy( sql, "SELECT InitSpatialMetadata()" );
  ret = sqlite3_exec( sqlite_handle, sql, NULL, NULL, &zErrMsg );
  if ( ret != SQLITE_OK )
  {
    string errCause = "Unable to initialize SpatialMetadata:\n" ;
    errCause += ( zErrMsg );
    fprintf( logFile, "SpatiaLite Database %s \n" , errCause );
    sqlite3_free( zErrMsg );
    return;
  }
  spatial_ref_sys_init( sqlite_handle, 0 );
}


bool createDb(string dbfName, sqlite3* &sqlite_handle, FILE* logFile)
{
  int ret;
  char *zErrMsg = NULL;

  if ( dbfName.length()==0 )
    return false;

    fprintf(logFile, "creating a spatialite db %s \n" ,dbfName.c_str() );


    // creating/opening the new database
    spatialite_init( 0 );

    ret = sqlite3_open_v2( dbfName.c_str(), &sqlite_handle, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );
    if ( ret!=SQLITE_OK )
    {
      // an error occurred
      string errCause =  "Could not create a new database \n" ;
      errCause += ( sqlite3_errmsg( sqlite_handle ) );
      sqlite3_close( sqlite_handle );
      fprintf(logFile,  "SpatiaLite Database %s \n" , errCause.c_str() );
      return false;
    }
    // activating Foreign Key constraints
    ret = sqlite3_exec( sqlite_handle, "PRAGMA foreign_keys = 1", NULL, 0, &zErrMsg );
    if ( ret != SQLITE_OK )
    {
      fprintf( logFile,  "SpatiaLite Database Unable to activate FOREIGN_KEY constraints %s \n", zErrMsg  );
      sqlite3_free( zErrMsg );
      sqlite3_close( sqlite_handle );
      return false;
    }
    initializeSpatialMetadata( sqlite_handle, logFile );

    // all done: closing the DB connection
    sqlite3_close( sqlite_handle );

  return true;
}


// create a spatial table using a SQL Statement 

template <typename q, typename t, typename d,typename s,typename o >
bool createSpaTbl(q& strSQL, t& tblName, d& db,s& srid, o& logFile)
{
	int ret;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmt = NULL;

	if ( strSQL.length()==0 )
    return false;

    logFile<< "creating a spatialite table %s \n" <<tblName<<endl;
	//fprintf(logFile, "\nCreating a spatialite table : %s , in Database  %s ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );

// CREATE TABLE "RTD_Streets_DC_Vx" ( "PK_UID" INTEGER PRIMARY KEY AUTOINCREMENT,
// "VERTID" INTEGER  PRIMARY KEY, "EDGEID" INTEGER, "IDP" INTEGER, "COST" DOUBLE,
// "SIGNAL" INTEGER, "Stop" INTEGER, "Geometry" POINT)

//CREATE TABLE "RTD_Streets_DC_Edge" ( "EDGEID" INTEGER PRIMARY KEY, "EDGESID" INTEGER,
//"VX1" INTEGER, "VX2" INTEGER, "LENGTH" DOUBLE, "STREETNAME" TEXT, "COST" DOUBLE, 
//"SCOST" DOUBLE, "ECOST" DOUBLE, "COMMENT" TEXT, "Stop" INTEGER, "Geometry" LINESTRING)
		// create collection table
		string strtblGeo="Drop Table " + tblName + "; "; 
		if ( sqlite3_exec( db, strtblGeo.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				//fprintf(logFile, "\nTable : %s\n is dropped ! \n", tblName.c_str() ,  sqlite3_errmsg( db) );
				logFile<<"Table "<<tblName<< " dropped!"<<endl;
		} else {
				logFile<<"\nTable : "<<tblName<<" is not dropped ! \n"<< sqlite3_errmsg( db) ;
		}
		strtblGeo="Create Table " + tblName + "  " +  strSQL ;
		if ( sqlite3_prepare_v2( db, strtblGeo.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL error: "<< strtblGeo<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			//fprintf(logFile,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", strtblGeo.c_str(), sqlite3_errmsg( db));
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", strtblGeo.c_str(), sqlite3_errmsg( db));
			return false;
		}
		if ( sqlite3_exec( db, strtblGeo.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
			logFile<<"\tSQL executed successfully: \t "<< strtblGeo<< " \n "<<endl;
			//fprintf(logFile,"\tSQL executed successfully: \t %s \n ", strtblGeo.c_str());
			fprintf(stdout,"\tSQL executed successfully: \t %s \n ", strtblGeo.c_str());

		} else {
				logFile<<"\tSQL %s failed to execute successfully: \n\t "<<strSQL<<endl<<to_string<const char *>( sqlite3_errmsg( db));
				//fprintf(logFile,"\tSQL %s failed to execute successfully: \n\t %s \n ", strSQL.c_str(), sqlite3_errmsg( db));
				fprintf(stdout,"\tSQL %s failed to execute successfully: \n\t %s \n ", strSQL.c_str(), sqlite3_errmsg( db));

		}

	sqlite3_finalize(stmt );			
	return true;
}


// create a View given a view name and SQL Statement 

template <typename q, typename t, typename d,typename o >
bool createView(q& strSQL, t& vwName, d& db, o& logFile)
{
	int ret;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmt = NULL;

	if ( strSQL.length()==0 )
    return false;

    logFile<< "creating view  %s \n" <<vwName<<endl;

// CREATE View "RTD_Streets_DC_Vx" as Select "PK_UID" , "VERTID" , "EDGEID" , "COST" ,
// "SIGNAL" , "Stop" , "Geometry" )

		// create collection table
		string strtblGeo="Drop View " + vwName + "; "; 
		if ( sqlite3_exec( db, strtblGeo.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
				logFile<<"View "<<vwName<< " dropped!"<<endl;
		}
		strtblGeo="Create View " + vwName + " as " +  strSQL ;
		if ( sqlite3_prepare_v2( db, strtblGeo.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL error: "<< strtblGeo<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", strtblGeo.c_str(), sqlite3_errmsg( db));
			return false;
		}
		if ( sqlite3_exec( db, strtblGeo.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
		{
			logFile<<"\tSQL executed successfully: \t "<< strtblGeo<< " \n "<<endl;
			fprintf(stdout,"\tSQL executed successfully: \t %s \n ", strtblGeo.c_str());

		} else {
				logFile<<"\tSQL %s failed to execute successfully: \n\t "<<strSQL<<endl<<to_string<const char *>( sqlite3_errmsg( db));
				fprintf(stdout,"\tSQL %s failed to execute successfully: \n\t %s \n ", strSQL.c_str(), sqlite3_errmsg( db));

		}

	sqlite3_finalize(stmt );			
	return true;
}



// Execute spatial query using a SQL Statement 

template <typename q, typename t, typename d,typename s,typename o >
bool execSpatialQuery(q& strSQL, t& tblName, d& db,s& srid, o& logFile)
{
	int ret;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmt = NULL;

	if ( strSQL.length()==0 )
    return false;

    logFile<< "Executing query %s for table  %s \n" <<strSQL<<tblName<<endl;

	//strSQL = " SELECT recovergeometrycolumn('" + tblName + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Point',2);";
		if ( sqlite3_prepare_v2( db, strSQL.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite %s prepare error. \tSQL %s error:" << strSQL<<endl<<to_string<const char *>( sqlite3_errmsg( db));
			fprintf(stderr,"SQLite %s prepare error: \tSQL  error : %s .\n", strSQL.c_str(), sqlite3_errmsg( db));
			return false;
		} else {
			logFile<<"\tSQL %s prepared successfully. \t "<<strSQL<<endl<< to_string<const char *>(sqlite3_errmsg( db));
			fprintf(stdout,"\tSQL %s prepared successfully. \t %s \n ", strSQL.c_str(), sqlite3_errmsg( db));

			if ( sqlite3_exec( db, strSQL.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
			{
				logFile<<"\tSQL %s executed successfully. "<< strSQL<<endl<<to_string<const char *>( sqlite3_errmsg( db));
				fprintf(stdout,"\tSQL %s executed successfully. %s . \n ", strSQL.c_str(), sqlite3_errmsg( db));
			}

		}

	sqlite3_finalize(stmt );			
	return true;
}


// Recover spatial geometry data table using a SQL Statement 

template <typename q, typename t, typename d,typename s,typename o >
bool recoverSpatialGeometry(q& strSQL, t& tblName, d& db,s& srid, o& logFile)
{
	int ret;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmt = NULL;

	if ( strSQL.length()==0 )
    return false;

    logFile<< "recovering spatialite geometry for table  %s \n" <<tblName<<endl;

	//strSQL = " SELECT recovergeometrycolumn('" + tblName + "', 'Geometry'," + (to_string<int>(srid)) + " ,'Point',2);";
		if ( sqlite3_prepare_v2( db, strSQL.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite %s prepare error. \tSQL %s error:" << strSQL<<endl<<to_string<const char *>( sqlite3_errmsg( db));
			fprintf(stderr,"SQLite %s prepare error: \tSQL  error : %s .\n", strSQL.c_str(), sqlite3_errmsg( db));
			return false;
		} else {
			logFile<<"\tSQL %s prepared successfully. \t "<<strSQL<<endl<< to_string<const char *>(sqlite3_errmsg( db));
			fprintf(stdout,"\tSQL %s prepared successfully. \t %s \n ", strSQL.c_str(), sqlite3_errmsg( db));

			if ( sqlite3_exec( db, strSQL.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
			{
				logFile<<"\tSQL %s executed successfully. "<< strSQL<<endl<<to_string<const char *>( sqlite3_errmsg( db));
				fprintf(stdout,"\tSQL %s executed successfully. %s . \n ", strSQL.c_str(), sqlite3_errmsg( db));
			}

		}

	sqlite3_finalize(stmt );			
	return true;
}




// populate a spatial table with result Vertex data  

template <typename q,  typename d, typename o,typename m,typename t, typename k, typename s, typename f >
bool inSpaTblVxRange(q& qrySQL, d& db,o& o1, m& m0,m& m1,t& tblNm,k& keyFld,s& srid, f& logFile)
{
	int ret=0, col=0;
	bool blnIns=false;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	op1=&o1;
	o o2;
	o o3;
	m::iterator it0;
	m::iterator it1;
	typedef m::iterator mit;
	mit mit1;
	mit mit2;

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		int id = 0;
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			o1 = it0->second;
			id = it0->second.get_id();
			//it1 = m1.find(id);
			//if(it1 !=m1.end()) {
			//	o2 = it1->second;
			//}
		    pair<mit, mit> vx1Range = m1.equal_range(id);
			size_t j = distance(vx1Range.first,vx1Range.second);
			recount++;
			checount++;
			if (j>=1) { 
				for (mit1=vx1Range.first; mit1!=vx1Range.second;mit1++)
				{
					o2 = mit1->second;
					// call the insert routine
					blnIns = exInSpaTblVx(db,stmti,o1,o2,logFile);
					if (blnIns) {
						qrySQL = "Update " + tblNm + " set geometry = MakePoint(" 
							" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
							" " + to_string<s> (srid ) + " ); " ; // where " + keyFld + " = " + to_string<int> (id ) + " ;";
						if ( sqlite3_exec( db, qrySQL.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
						{
							fprintf(stdout,"\n Geometry field MakePoint failed for SQL %s  \n Error :  %s \n", qrySQL.c_str(), sqlite3_errmsg(db)) ;
							logFile<<"Geometry field MakePoint update failed!" <<endl<<" Sql command : "<<qrySQL<<endl<<"Error : "<<to_string<const char *>(sqlite3_errmsg(db))<<endl;
						}

				// if insert executed without error	
					} else { // if insert casued error
						logFile<<"\n\tInsert SQL "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
						sqlite3_free(zErrMsg);
						return false;
					}
				} // possible mulitple on stops per off stop

			} // Check if there are any off stops for the current stop for each on vertex
			else
			{
				o2 = o3; // reinitialize o2 to a new object
				blnIns = exInSpaTblVx(db,stmti,o1,o2,logFile);
				if (blnIns) {
					qrySQL = "Update " + tblNm + " set geometry = MakePoint(" 
						" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
						" " + to_string<s> (srid ) + " ); " ; // where " + keyFld + " = " + to_string<int> (id ) + " ;";
					if ( sqlite3_exec( db, qrySQL.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
					{
						fprintf(stdout,"\n Geometry field MakePoint failed for SQL %s  \n Error :  %s \n", qrySQL.c_str(), sqlite3_errmsg(db)) ;
						logFile<<"Geometry field MakePoint update failed!" <<endl<<" Sql command : "<<qrySQL<<endl<<"Error : "<<to_string<const char *>(sqlite3_errmsg(db))<<endl;
					}
				}
			}
			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold

		}	// for each on stop
		// update the geometry column


		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}

		sqlite3_finalize(stmti );			
	return true;

}

template <typename q,  typename d, typename o,typename m,typename t, typename k, typename s, typename f >
bool inSpaTblVx(q& qrySQL, d& db,o& o1, m& m0,m& m1,t& tblNm,k& keyFld,s& srid, f& logFile)
{
	int ret=0, col=0;
	bool blnIns=false;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	op1=&o1;
	o o2;
	o o3;
	m::iterator it0;
	m::iterator it1;
	typedef m::iterator mit;
	mit mit1;
	mit mit2;

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		int id = 0;
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			o1 = it0->second;
			id = it0->second.get_id();
			it1 = m1.find(id);
			if(it1 !=m1.end()) {
				o2 = it1->second;
			} else {
				o2=o3;
			}
			recount++;
			checount++;
			blnIns = exInSpaTblVx(db,stmti,o1,o2,logFile);
				if (recount==100000)
				{
					int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
					if ( ret != SQLITE_OK )
					{
					// some error occurred
						logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
						sqlite3_free(zErrMsg);
						return false;
					}

					fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
					logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

					ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
					if ( ret != SQLITE_OK )
					{
					// some error occurred
						logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
						sqlite3_free(zErrMsg);
						return false;
					}
					recount = 0;
				}// if insert count is greater than threshold
		}
		// update the geometry column

		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}

		sqlite3_finalize(stmti );			
	return true;

}


template <typename d,  typename s, typename o, typename f >
bool exInSpaTblVx(d& db, s& stmti, o& o1, o& o2, f& logFile)
{
	int ret;
	int col=1, id =0;
	char *zErrMsg = NULL;
	o* op1;
	op1=&o1;
	o* op2;
	op2=&o2;

	sqlite3_reset( stmti );
	sqlite3_clear_bindings( stmti );
	col=1;
	id = o1.get_id();
	sqlite3_bind_int(stmti,col,id);
	col=2;
	sqlite3_bind_int(stmti,col,o1.get_idp());
	col=3;
	sqlite3_bind_int(stmti,col,o1.get_orig());
	col=4;
	sqlite3_bind_double(stmti,col,o1.get_cost());
	col=5;
	sqlite3_bind_int(stmti,col,o1.get_index());
	col=6;
	sqlite3_bind_int(stmti,col,o1.get_lowlink());
	col=7;
	sqlite3_bind_int(stmti,col,o2.get_orig());
	col=8;
	sqlite3_bind_double(stmti,col,o2.get_cost());
	col=9;
	sqlite3_bind_int(stmti,col,o2.get_index());
	col=10;
	sqlite3_bind_int(stmti,col,o2.get_lowlink());
//			sqlite3_bind_blob(stmti,11,geotxt)
	ret = sqlite3_step(stmti);
	if (ret != SQLITE_DONE)
	{
		logFile<<"Database Insert failed for id "<<id<< " !"<< endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
		sqlite3_free(zErrMsg);
		return false;

	} // if insert is unsuccessful

	return true;

} // update the statement and execute the vertex insert


template <typename q,  typename d, typename o,typename m,typename t, typename k, typename s, typename f >
bool inSpaTblEdge(q& qrySQL, d& db,o& o1, m& m0,m& m1,t& tblNm,k& keyFld,s& srid, f& logFile)
{
	int ret=0,col=1;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	o* op2;
	o* op3;
	op1=&o1;
	o o2;
	op2=&o2;
	o o3;
	op3=&o3;
	m::iterator it0;
	m::iterator it1;

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		int id = 0, esid = 0;
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			o1 = it0->second;
			esid = it0->second.get_esid();
			it1 = m1.find(esid);
			if(it1 !=m1.end()) {
				o2 = it1->second;
			} else {
				o2 = o3;
			}
			if( (o1.get_orig()!=o2.get_orig()) && (o1.get_tcost()==o2.get_tcost()) && (o1.get_orig() != -2 || o2.get_orig() != -2) ) {
				// some error occurred
				logFile<<"\n\tAccess Time Error, Edge 1: " <<to_string<long>(o1.get_id())<<"\t Stop : "<<to_string<long>(o1.get_orig());
				logFile<<"\t Stop On : "<<to_string<double>(o1.get_stopon())<<"\t Off : "<<to_string<double>(o1.get_stopoff())<<"\t Dir : "<<to_string<double>(o1.get_dirn());
				logFile<<"\t Vx From : "<<to_string<double>(o1.get_frid())<<"\t to : "<<to_string<double>(o1.get_toid());
				logFile<<"\t Cost : "<<to_string<double>(o1.get_tcost())<<"\t Edge 2 : "<<to_string<long>(o2.get_id());
				logFile<<"\t Stop On : "<<to_string<double>(o2.get_stopon())<<"\t Off : "<<to_string<double>(o2.get_stopoff())<<"\t Dir : "<<to_string<double>(o2.get_dirn());
				logFile<<"\t Vx From : "<<to_string<double>(o2.get_frid())<<"\t to : "<<to_string<double>(o2.get_toid());
				logFile<<"\t Stop : "<<to_string<long>(o2.get_orig())<<"\t Cost : "<<to_string<double>(o2.get_tcost())<<endl;
				fprintf(stderr,"\n\t Access Time error, Edge 1: %d , Stop : %d , Stop On : %d , off : %d , Dirn : %d , Vx From : %d , To : %d , Cost : %12.3f , Edge 2: %d , Stop : %d , Stop On : %d , off : %d , Dirn : %d , Vx From : %d , To : %d , Cost : %12.3f \n", o1.get_id(),o1.get_orig(),o1.get_stopon(),o1.get_stopoff(),o1.get_dirn(),o1.get_frid(),o1.get_toid(),o1.get_tcost(),o2.get_id(),o2.get_orig(),o2.get_stopon(),o2.get_stopoff(),o2.get_dirn(),o2.get_frid(),o2.get_toid(),o2.get_tcost());
			}
			recount++;
			checount++;
			sqlite3_reset( stmti );
			sqlite3_clear_bindings( stmti );
			id = o1.get_id();
			col=1;
			sqlite3_bind_int(stmti,col,id);
			col=2;
			sqlite3_bind_int(stmti,col,o1.get_esid());
			col=3;
			sqlite3_bind_int(stmti,col,o1.get_orig());
			col=4;
			sqlite3_bind_double(stmti,col,o1.get_tcost());
			col=5;
			sqlite3_bind_int(stmti,col,o1.get_lbl());  // edge labeled
			col=6;
			sqlite3_bind_int(stmti,col,o1.get_dirn()); // direction
			col=7;
			sqlite3_bind_int(stmti,col,o2.get_orig());
			col=8;
			sqlite3_bind_double(stmti,col,o2.get_tcost());
			col=9;
			sqlite3_bind_int(stmti,col,o2.get_lbl());
			col=10;
			sqlite3_bind_int(stmti,col,o2.get_dirn());
			int ret = sqlite3_step(stmti);
			if (ret != SQLITE_DONE)
			{
				logFile<<"SQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
				sqlite3_free(zErrMsg);
				ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				return false;

			} // if insert is successful
			// update the geometry column
			//string strGeo = "Update " + tblNm + " set geometry = ( Select Geometry " 
			//	" from " + dbATable + " where tblNm.EdgeId = dbATable.EdgeId); "  ;
			//if ( sqlite3_exec( db, strGeo.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
			//{
			//		logFile<<"Geometry field update failed for table "<<tblNm<<" !"<<endl;
			//}

			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold
		}
		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}
	sqlite3_finalize(stmti );			
	return true;

}

template <typename q,  typename d, typename o,typename m,typename t, typename k, typename s, typename f, typename h>
bool inSpaTblParcel(q& qrySQL, d& db,o& o1, m& m0,t& tblNm,k& keyFld,s& srid, f& logFile,h& blnHi)
{
	int ret=0,col=1;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	op1=&o1;
	o o2;
	m::iterator it0;
	m::iterator it1;

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		int id = 0;
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			o1 = it0->second;
			//o1.set_hist(blnHi);
			id = it0->second.get_id();
			recount++;
			checount++;
			sqlite3_reset( stmti );
			sqlite3_clear_bindings( stmti );
			id = o1.get_id();
			col=1;
			sqlite3_bind_int(stmti,col,id);
			string pacid = o1.get_pacid();
			col=2;
			sqlite3_bind_text(stmti,col,pacid.c_str(),pacid.length(),NULL);
			col=3;
			sqlite3_bind_int(stmti,col,o1.get_origon()); // stop for ons
			col=4;
			sqlite3_bind_double(stmti,col,o1.get_ons());
			col=5;
			if (o1.get_hist()) {
				sqlite3_bind_double(stmti,col,o1.get_hwkoncost());
			} else {
				sqlite3_bind_double(stmti,col,o1.get_wkoncost());
			}

			col=6;
			sqlite3_bind_int(stmti,col,o1.get_lbl());  // parcel labeled
			col=7;
			sqlite3_bind_int(stmti,col,o1.get_origoff()); // stop for offs
			col=8;
			sqlite3_bind_double(stmti,col,o1.get_offs());
			col=9;
			if (o1.get_hist()) {
				sqlite3_bind_double(stmti,col,o1.get_hwkoffcost());
			} else {
				sqlite3_bind_double(stmti,col,o1.get_wkoffcost());
			}
			col=10;
			sqlite3_bind_int(stmti,col,o1.get_lbl());  // parcel labeled
//			sqlite3_bind_blob(stmti,7,geotxt)
			int ret = sqlite3_step(stmti);
			if (ret != SQLITE_DONE)
			{
				logFile<<"SQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
				sqlite3_free(zErrMsg);
				return false;
				ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

			} // if insert is successful
			// update the geometry column
			//string strGeo = "Update " + tblNm + " set geometry = MakePoint(" 
			//	" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
			//	" " + to_string<s> (srid ) + " ) where " + keyFld + " = " + to_string<int> (id ) + " ;";
			//if ( sqlite3_exec( db, strGeo.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
			//{
			//		logFile<<"Geometry field update failed for table "<<tblNm<<" !"<<endl;
			//}

			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold
		} // loop over all objects
		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}
	sqlite3_finalize(stmti );			
	return true;

}

template <typename q,  typename d, typename o,typename m,typename t, typename k, typename s, typename f >
bool inSpaTblStop(q& qrySQL, d& db,o& o1, m& m0,t& tblNm,k& keyFld,s& srid, f& logFile)
{
	int ret=0, col=0;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	//tstop o3;
	op1=&o1;
	m::iterator it0;
	m::iterator it1;
	double dObjVal=0;
	int iObjVal=0;
	long lObjVal=0;
	string sObjVal="",bObjVal="";

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			cout<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		col =	sqlite3_column_count( stmti );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			o1 = it0->second;
			lObjVal = it0->second.get_id();
			recount++;
			checount++;
			sqlite3_reset( stmti );
			sqlite3_clear_bindings( stmti );
			col = 1; // id
			lObjVal = op1->get_id();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 2; //tripId
			lObjVal = op1->tripId();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 3; //stopId
			sObjVal = (op1->get_StopLbl());
			if (sObjVal.length()>0) {
				sqlite3_bind_int(stmti,col, from_string<long>( sObjVal));
			}
			//if (sObjVal.length()>0) {
			//	sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			//}
			col = 4;  //StOrdr
			iObjVal = op1->get_StOrdr();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 5;  //StopName
			sObjVal = op1->get_StopName();
			if (sObjVal.length()>0) {
				sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			}
			col = 6; // Predecessor id
			lObjVal = op1->get_Stopidp();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 7;  // Successor id
			lObjVal = op1->get_Stopids();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 8;  // Edeg Id
			lObjVal = op1->get_Edgeid();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 9; // Position Along
			dObjVal = op1->get_posalong();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 10; // if this is labeled
			iObjVal = op1->get_lbl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 11; // if historic stop
			iObjVal = op1->get_blnHist();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 12;  // if inbound
			iObjVal = op1->get_blnInbd();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 13; // external stop
			iObjVal = op1->get_blnExtr();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 14;  // included in this analysis
			iObjVal = op1->get_blnIncl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 15;  // is eliminatable
			iObjVal = op1->get_blnElim();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 16;  // Cumulative Distance 
			dObjVal = op1->get_CumDist();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 17;  // set cumulative ride time 
			dObjVal = op1->get_CRdTm();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 18;
			dObjVal = op1->get_undCRdTm();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 19;
			dObjVal = op1->get_CRdTmC();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 20;
			dObjVal = op1->get_HistOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 21;
			dObjVal = op1->get_HistOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 22;
			dObjVal = op1->get_Ons();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 23;
			dObjVal = op1->get_Offs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 24;
			dObjVal = op1->get_DepVol();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 25;
			dObjVal = op1->get_probStoph();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 26;
			dObjVal = op1->get_probStop();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 27;
			dObjVal = op1->get_depDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 28;
			dObjVal = op1->get_arrDelay();
			sqlite3_bind_double(stmti,col,dObjVal);


			col = 29;
			dObjVal = ((op1->get_depDelay() + op1->get_arrDelay())*op1->get_probStop());
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 30;
			dObjVal = op1->get_dwellDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 31;
			dObjVal = op1->get_rideDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 32;
			dObjVal = op1->get_PVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 33;
			dObjVal = op1->get_AVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 34;
			dObjVal = op1->get_CRdTmE();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 35;
			dObjVal = op1->get_hWkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 36;
			dObjVal = op1->get_hWkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 37;
			dObjVal = op1->get_WkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 38;
			dObjVal = op1->get_WkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 39;
			dObjVal = op1->get_WalkCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 40;
			dObjVal = op1->get_RideCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 41;
			dObjVal = op1->get_OperCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 42;
			dObjVal = op1->get_TCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 43;
			dObjVal = op1->get_xc();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 44;
			dObjVal = op1->get_yc();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 45;
			dObjVal = 0;
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 46;
			bObjVal = "MakePoint(" + to_string<double>(op1->get_xc()) + "," + to_string<double>(op1->get_yc()) + "," + to_string<int>(srid) + ")"; //op1->get_zc();
			sqlite3_bind_blob(stmti,col,bObjVal.c_str(),bObjVal.length(),NULL);

//			sqlite3_bind_blob(stmti,7,geotxt)
			int ret = sqlite3_step(stmti);
			if (ret != SQLITE_DONE)
			{
				logFile<<"SQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
				fprintf(stderr,"\n\tSQL %s \n error: %s \n", qrySQL.c_str(), sqlite3_errmsg(db));
				ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				sqlite3_free(zErrMsg);
				return false;

			} // if insert is successful
			// update the geometry column
			//string strGeo = "Update " + tblNm + " set geometry = MakePoint(" 
			//	" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
			//	" " + to_string<s> (srid ) + " ) where " + keyFld + " = " + to_string<int> (id ) + " ;";
			//if ( sqlite3_exec( db, strGeo.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
			//{
			//		logFile<<"Geometry field update failed for table "<<tblNm<<" !"<<endl;
			//}

			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold
		}
		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}
	sqlite3_finalize(stmti );
	return true;

}



template <typename q,  typename d, typename m,typename o,typename p,typename g, typename f >
bool inSpaTblDPStop(q& qrySQL, d& db,m& m0,o& o1,p& p1,g& g1, f& logFile)
{
	int ret=0, col=0;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	p* dp1;
	g* pg1;
	//tstop o3;
	op1=&o1;
	dp1 =&p1;
	//pc1 = c1;
	m::iterator it0;
	m::iterator it1;
	double dObjVal=0;
	int iObjVal=0;
	long lObjVal=0;
	string sObjVal="",siObjVal="",skObjVal="",snObjVal="";

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		col =	sqlite3_column_count( stmti );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			o1 = it0->second;
			recount++;
			checount++;
			sqlite3_reset( stmti );
			sqlite3_clear_bindings( stmti );
			// DP Key 
			col = 1;
			skObjVal = o1.get_dpkey();
			if (skObjVal.length()>0) {
				sqlite3_bind_text(stmti,col,skObjVal.c_str(),skObjVal.length(),NULL);
			}
			// StopId 
			col = 2;
			lObjVal = from_String<long>(o1.get_tstop().get_StopLbl());
			sqlite3_bind_int(stmti,col,lObjVal);
//			if (siObjVal.length()>0) {
//				sqlite3_bind_text(stmti,col,siObjVal.c_str(),siObjVal.length(),NULL);
//			}
			// StopName 
			col = 3;
			snObjVal = o1.get_tstop().get_StopName();
			if (snObjVal.length()>0) {
				sqlite3_bind_text(stmti,col,snObjVal.c_str(),snObjVal.length(),NULL);
			}
			// walk time coefficient
			col = 4;
			dObjVal = g1->get_walkcost();
			sqlite3_bind_double(stmti,col,dObjVal);

			// ride time coefficient
			col = 5;
			dObjVal = g1->get_ridecost();
			sqlite3_bind_double(stmti,col,dObjVal);

			// scenario id q
			col = 6;
			lObjVal = o1.get_q();
			sqlite3_bind_int(stmti,col,lObjVal);
			// i of dp
			col = 7;
			lObjVal = o1.get_i();
			sqlite3_bind_int(stmti,col,lObjVal);

			// j of dp
			col = 8;
			lObjVal = o1.get_j();
			sqlite3_bind_double(stmti,col,lObjVal);

			// k of dp
			col = 9;
			lObjVal = o1.get_k();
			sqlite3_bind_double(stmti,col,lObjVal);

			// l of dp
			col = 10;
			lObjVal = o1.get_l();
			sqlite3_bind_double(stmti,col,lObjVal);

			// m of dp
			col = 11;
			lObjVal = o1.get_m();
			sqlite3_bind_double(stmti,col,lObjVal);

			// Cum Ride Time
			col = 12;
			dObjVal = o1.get_tstop().get_CRdTm();
			sqlite3_bind_double(stmti,col,dObjVal);
			// und Cum Ride Time 
			col = 13;
			dObjVal = o1.get_tstop().get_undCRdTm();
			sqlite3_bind_double(stmti,col,  dObjVal);

			col = 14;
			dObjVal = o1.get_tstop().get_Ons();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 15;
			dObjVal = o1.get_tstop().get_Offs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 16;
			dObjVal = o1.get_tstop().get_DepVol();
			sqlite3_bind_double(stmti,col,dObjVal);


			col = 17;
			dObjVal = o1.get_tstop().get_probStop();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 18;
			dObjVal = o1.get_tstop().get_depDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 19;
			dObjVal = o1.get_tstop().get_arrDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 20;
			dObjVal = o1.get_tstop().get_dwellDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 21;
			dObjVal = o1.get_tstop().get_rideDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 22;
			dObjVal = o1.get_tstop().get_PVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 23;
			dObjVal = o1.get_tstop().get_AVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 24;
			dObjVal = o1.get_tstop().get_WkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 25;
			dObjVal = o1.get_tstop().get_WkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 26;
			dObjVal = o1.get_tstop().get_WalkCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 27;
			dObjVal = o1.get_tstop().get_RideCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 28;
			dObjVal = o1.get_tstop().get_OperCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 29;
			dObjVal = o1.get_tstop().get_TCost();
			sqlite3_bind_double(stmti,col,dObjVal);
			col = 30; // dpcount if multi period run otherwise = 1
			lObjVal = 1;//o1.get_dpcnt();
			sqlite3_bind_int(stmti,col,lObjVal);
			col = 31;
			dObjVal = (o1.get_tstop().get_arrDelay()+o1.get_tstop().get_depDelay())*o1.get_tstop().get_probStop();
			sqlite3_bind_double(stmti,col,dObjVal);
			col = 32;  // StOrdr
			lObjVal = o1.get_tstop().get_StOrdr();
			sqlite3_bind_int(stmti,col,lObjVal);

			// Computed Cum Ride Time 
			col = 33;
			dObjVal = o1.get_tstop().get_CRdTmC();
			sqlite3_bind_double(stmti,col,dObjVal);


			// Stop Spacing Cum Distance Time 
			col = 34;
			dObjVal = o1.get_tstop().get_CumDist();
			sqlite3_bind_double(stmti,col,dObjVal);

//			sqlite3_bind_blob(stmti,7,geotxt)
			ret = sqlite3_step(stmti);
			if (ret != SQLITE_DONE)
			{
				logFile<<"SQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
				fprintf(stderr,"\n\tSQL %s \n error: %s \n", qrySQL.c_str(), sqlite3_errmsg(db));
				ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				sqlite3_free(zErrMsg);
				return false;

			} // if insert is successful
			// update the geometry column
			//string strGeo = "Update " + tblNm + " set geometry = MakePoint(" 
			//	" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
			//	" " + to_string<s> (srid ) + " ) where " + keyFld + " = " + to_string<int> (id ) + " ;";
			//if ( sqlite3_exec( db, strGeo.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
			//{
			//		logFile<<"Geometry field update failed for table "<<tblNm<<" !"<<endl;
			//}

			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold
		}
		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}
	sqlite3_finalize(stmti );
	return true;

}







template <typename q,  typename d, typename o,typename p,typename g,typename m,typename c,typename i, typename s, typename f >
bool inSpaTblTripDPStop(q& qrySQL, d& db,o& o1,p& p1, g& g1, m& m0,c c1, i& i1,s& s1, f& logFile)
{
	int ret=0, col=0;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;
	p* dp1;
	g* pg1;
	c* pc1;
	//tstop o3;
	op1=&o1;
	dp1 =&p1;
	//pc1 = c1;
	m::iterator it0;
	m::iterator it1;
	double dObjVal=0;
	int iObjVal=0;
	long lObjVal=0;
	string sObjVal="";

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		col =	sqlite3_column_count( stmti );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			p1 = it0->second;
			recount++;
			checount++;
			sqlite3_reset( stmti );
			sqlite3_clear_bindings( stmti );
			// tripId 
			col = 1;
			lObjVal = i1; // op1->get_tripId();
			sqlite3_bind_int(stmti,col,lObjVal);
			// walk time coefficient
			col = 2;
			dObjVal = c1->get_walkcost();
			sqlite3_bind_double(stmti,col,dObjVal);

			// ride time coefficient
			col = 3;
			dObjVal = c1->get_ridecost();
			sqlite3_bind_double(stmti,col,dObjVal);
			//// dpKey
			//col = 4;
			//sObjVal = op1->get_tDPkey();
			//if (sObjVal.length()>0) {
			//	sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			//}
			// use the dp object data to populate the table
			//p1 = op1->get_dptStop();

			// scenario id q
			col = 4;
			lObjVal = p1.get_q();
			sqlite3_bind_int(stmti,col,lObjVal);
			// i of dp
			col = 5;
			lObjVal = p1.get_i();
			sqlite3_bind_int(stmti,col,lObjVal);

			// j of dp
			col = 6;
			lObjVal = p1.get_j();
			sqlite3_bind_double(stmti,col,lObjVal);

			// k of dp
			col = 7;
			lObjVal = p1.get_k();
			sqlite3_bind_double(stmti,col,lObjVal);

			// l of dp
			col = 8;
			lObjVal = p1.get_l();
			sqlite3_bind_double(stmti,col,lObjVal);

			// m of dp
			col = 9;
			lObjVal = p1.get_m();
			sqlite3_bind_double(stmti,col,lObjVal);

			// get the stop data from the dp result object
			g1 = p1.get_tstop();
			// id
			col = 10;
			lObjVal = g1.get_id();
			sqlite3_bind_int(stmti,col,lObjVal);
			// StopLbl 
			col = 11;
			sObjVal = g1.get_StopLbl();
			if (sObjVal.length()>0) {
				sqlite3_bind_int(stmti,col, from_string<long>( sObjVal));
			}
			//if (sObjVal.length()>0) {
			//	sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			//}

			col = 12;
			iObjVal = g1.get_StOrdr();
			sqlite3_bind_int(stmti,col,iObjVal);


			// StopName 
			col = 13;
			sObjVal = g1.get_StopName();
			if (sObjVal.length()>0) {
				sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			}
			// Stop id Predecessor
			col = 14;
			lObjVal = g1.get_Stopidp();
			sqlite3_bind_int(stmti,col,lObjVal);
			// stop id successor
			col = 15;
			lObjVal = g1.get_Stopids();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 16;
			lObjVal = g1.get_Edgeid();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 17;
			dObjVal = g1.get_posalong();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 18;
			iObjVal = g1.get_lbl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 19;
			iObjVal = g1.get_lbl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 20;
			iObjVal = g1.get_blnInbd();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 21;
			iObjVal = g1.get_blnExtr();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 22;
			iObjVal = g1.get_blnIncl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 23;
			iObjVal = g1.get_blnElim();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 24;
			dObjVal = g1.get_CumDist();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 25;
			dObjVal = g1.get_CRdTm();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 26;
			dObjVal = g1.get_undCRdTm();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 27;
			dObjVal = g1.get_CRdTmC();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 28;
			dObjVal = g1.get_HistOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 29;
			dObjVal = g1.get_HistOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 30;
			dObjVal = g1.get_Ons();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 31;
			dObjVal = g1.get_Offs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 32;
			dObjVal = g1.get_DepVol();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 33;
			dObjVal = g1.get_probStoph();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 34;
			dObjVal = g1.get_probStop();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 35;
			dObjVal = g1.get_depDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 36;
			dObjVal = g1.get_arrDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 37;
			dObjVal = g1.get_dwellDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 38;
			dObjVal = g1.get_rideDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 39;
			dObjVal = g1.get_PVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 40;
			dObjVal = g1.get_AVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 41;
			dObjVal = g1.get_CRdTmE();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 42;
			dObjVal = g1.get_hWkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 43;
			dObjVal = g1.get_hWkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 44;
			dObjVal = g1.get_WkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 45;
			dObjVal = g1.get_WkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 46;
			dObjVal = g1.get_WalkCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 47;
			dObjVal = g1.get_RideCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 48;
			dObjVal = g1.get_OperCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 49;
			dObjVal = g1.get_TCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 50;
			dObjVal = g1.get_xc();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 51;
			dObjVal = g1.get_yc();
			sqlite3_bind_double(stmti,col,dObjVal);

			//col = 52;
			//dObjVal = g1.get_zc();
			//sqlite3_bind_double(stmti,col,dObjVal);
			col = 52; // stop delay computed as (arrDelay+depDelay)*probStop
			dObjVal = g1.cstopDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

//			sqlite3_bind_blob(stmti,7,geotxt)
			int ret = sqlite3_step(stmti);
			if (ret != SQLITE_DONE)
			{
				logFile<<"SQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
				fprintf(stderr,"\n\tSQL %s \n error: %s \n", qrySQL.c_str(), sqlite3_errmsg(db));
				ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				sqlite3_free(zErrMsg);
				return false;

			} // if insert is successful
			// update the geometry column
			//string strGeo = "Update " + tblNm + " set geometry = MakePoint(" 
			//	" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
			//	" " + to_string<s> (srid ) + " ) where " + keyFld + " = " + to_string<int> (id ) + " ;";
			//if ( sqlite3_exec( db, strGeo.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
			//{
			//		logFile<<"Geometry field update failed for table "<<tblNm<<" !"<<endl;
			//}

			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold
		}
		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}
	sqlite3_finalize(stmti );
	return true;

}


template <typename q,  typename d, typename o,typename p,typename m,typename c,typename i, typename s, typename f >
bool inSpaTblDPStop2(q& qrySQL, d& db,o& o1,p& p1,  m& m0,c& c1, i& i1,s& s1, f& logFile)
{
	int ret=0, col=0;
	char *zErrMsg = NULL;
	sqlite3_stmt *stmti = NULL;
	o* op1;  // dp object
	p* dp1;  // stop object
	//tstop o3;
	op1=&o1;
	dp1 =&p1; 
	//pc1 = c1;
	m::iterator it0;
	m::iterator it1;
	double dObjVal=0;
	int iObjVal=0;
	long lObjVal=0;
	string sObjVal="";

	if ( qrySQL.length()==0 )
    return false;


		if ( sqlite3_prepare_v2( db, qrySQL.c_str(), -1, &stmti, NULL ) != SQLITE_OK )
		{
			// some error occurred
			logFile<<"SQLite prepare error: \tSQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
			fprintf(stderr,"\nSQLite prepare error: \tSQL %s \n error: %s .\n", qrySQL.c_str(), sqlite3_errmsg( db));
			return false;
		}


		int recount=0; // record count 
		int checount=0; // modulo
		ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg );
		col =	sqlite3_column_count( stmti );
		if ( ret != SQLITE_OK )
		{
		// some error occurred
			logFile<<"\n\tSQL  "<<endl<< "BEGIN TRANSACTION"<<"\n error: "<<to_string<const char *>( sqlite3_errmsg(db));
			fprintf(stderr,"\n\tSQL %s \n error: %s \n", "BEGIN TRANSACTION", sqlite3_errmsg(db));
			return false;
		}
		
		for(it0=m0.begin();it0!=m0.end();++it0)
		{
			p1 = it0->second;
			recount++;
			checount++;
			sqlite3_reset( stmti );
			sqlite3_clear_bindings( stmti );
			// get the stop data from the dp result object
			p1 = o1.get_tstop();
			// tripId 
			col = 1;
			lObjVal = i1; // op1->get_tripId();
			sqlite3_bind_int(stmti,col,lObjVal);
			// walk time coefficient
			col = 2;
			dObjVal = c1->get_walkcost();
			sqlite3_bind_double(stmti,col,dObjVal);

			// ride time coefficient
			col = 3;
			dObjVal = c1->get_ridecost();
			sqlite3_bind_double(stmti,col,dObjVal);
			//// dpKey
			//col = 4;
			//sObjVal = op1->get_tDPkey();
			//if (sObjVal.length()>0) {
			//	sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			//}
			// use the dp object data to populate the table
			//p1 = op1->get_dptStop();

			// scenario id q
			col = 4;
			lObjVal = o1.get_q();
			sqlite3_bind_int(stmti,col,lObjVal);
			// i of dp
			col = 5;
			lObjVal = o1.get_i();
			sqlite3_bind_int(stmti,col,lObjVal);

			// j of dp
			col = 6;
			lObjVal = o1.get_j();
			sqlite3_bind_double(stmti,col,lObjVal);

			// k of dp
			col = 7;
			lObjVal = o1.get_k();
			sqlite3_bind_double(stmti,col,lObjVal);

			// l of dp
			col = 8;
			lObjVal = o1.get_l();
			sqlite3_bind_double(stmti,col,lObjVal);

			// m of dp
			col = 9;
			lObjVal = o1.get_m();
			sqlite3_bind_double(stmti,col,lObjVal);

			col = 10;
			lObjVal = p1.get_id();
			sqlite3_bind_int(stmti,col,lObjVal);
			// StopName 
			col = 11;
			sObjVal = p1.get_StopName();
			if (sObjVal.length()>0) {
				sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			}
			// Stop id Predecessor
			col = 12;
			lObjVal = p1.get_Stopidp();
			sqlite3_bind_int(stmti,col,lObjVal);
			// stop id successor
			col = 13;
			lObjVal = p1.get_Stopids();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 14;
			lObjVal = p1.get_Edgeid();
			sqlite3_bind_int(stmti,col,lObjVal);

			col = 15;
			dObjVal = p1.get_posalong();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 16;
			iObjVal = p1.get_StOrdr();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 17;
			sObjVal = p1.get_StopLbl();
			if (sObjVal.length()>0) {
				sqlite3_bind_text(stmti,col,sObjVal.c_str(),sObjVal.length(),NULL);
			}
			col = 18;
			iObjVal = p1.get_lbl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 19;
			iObjVal = p1.get_lbl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 20;
			iObjVal = p1.get_blnInbd();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 21;
			iObjVal = p1.get_blnExtr();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 22;
			iObjVal = p1.get_blnIncl();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 23;
			iObjVal = p1.get_blnElim();
			sqlite3_bind_int(stmti,col,iObjVal);

			col = 24;
			dObjVal = p1.get_CumDist();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 25;
			dObjVal = p1.get_CRdTm();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 26;
			dObjVal = p1.get_CRdTmC();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 27;
			dObjVal = p1.get_HistOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 28;
			dObjVal = p1.get_HistOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 29;
			dObjVal = p1.get_Ons();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 30;
			dObjVal = p1.get_Offs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 31;
			dObjVal = p1.get_DepVol();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 32;
			dObjVal = p1.get_probStoph();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 33;
			dObjVal = p1.get_probStop();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 34;
			dObjVal = p1.get_depDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 35;
			dObjVal = p1.get_arrDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 36;
			dObjVal = p1.get_dwellDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 37;
			dObjVal = p1.get_rideDelay();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 38;
			dObjVal = p1.get_PVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 39;
			dObjVal = p1.get_AVal();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 40;
			dObjVal = p1.get_CRdTmE();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 41;
			dObjVal = p1.get_hWkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);
			col = 42;
			dObjVal = p1.get_hWkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 43;
			dObjVal = p1.get_WkTmOns();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 44;
			dObjVal = p1.get_WkTmOffs();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 45;
			dObjVal = p1.get_WalkCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 46;
			dObjVal = p1.get_RideCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 47;
			dObjVal = p1.get_OperCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 48;
			dObjVal = p1.get_TCost();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 49;
			dObjVal = p1.get_xc();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 50;
			dObjVal = p1.get_yc();
			sqlite3_bind_double(stmti,col,dObjVal);

			col = 51;
			dObjVal = p1.get_zc();
			sqlite3_bind_double(stmti,col,dObjVal);

//			sqlite3_bind_blob(stmti,7,geotxt)
			int ret = sqlite3_step(stmti);
			if (ret != SQLITE_DONE)
			{
				logFile<<"SQL  "<< qrySQL<<endl<<" error: "<<to_string<const char *>(sqlite3_errmsg( db))<<endl;
				fprintf(stderr,"\n\tSQL %s \n error: %s \n", qrySQL.c_str(), sqlite3_errmsg(db));
				ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				sqlite3_free(zErrMsg);
				return false;

			} // if insert is successful
			// update the geometry column
			//string strGeo = "Update " + tblNm + " set geometry = MakePoint(" 
			//	" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
			//	" " + to_string<s> (srid ) + " ) where " + keyFld + " = " + to_string<int> (id ) + " ;";
			//if ( sqlite3_exec( db, strGeo.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
			//{
			//		logFile<<"Geometry field update failed for table "<<tblNm<<" !"<<endl;
			//}

			if (recount==100000)
			{
				int ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}

				fprintf(stdout,"\n\t %d records updated using SQL %s ", checount,qrySQL.c_str() ) ;
				logFile<<"\n\t"<<checount << " records updated  " <<endl<<" using SQL "<<qrySQL <<endl;

				ret = sqlite3_exec( db, "BEGIN", NULL, NULL, &zErrMsg);
				if ( ret != SQLITE_OK )
				{
				// some error occurred
					logFile<<"\n\tSQL BEGIN "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
					sqlite3_free(zErrMsg);
					return false;
				}
				recount = 0;
			}// if insert count is greater than threshold
		}
		ret = sqlite3_exec( db, "COMMIT", NULL, NULL, &zErrMsg );
		if ( ret != SQLITE_OK )
		{	// some error occurred
			logFile<<"\n\tSQL COMMIT " <<qrySQL<<endl<<" Error: "<<to_string<const char *>( sqlite3_errmsg(db));
			sqlite3_free(zErrMsg);
			return false;
		}
	sqlite3_finalize(stmti );
	return true;

}

// attach a database to another opened sqlite database for query 
// provide main opened database  pointer , 
template <typename a,  typename b, typename c, typename f >
bool sqliteAttachDB(a& mDb, b& fNm,c& aNm, f& logFile)
{
	bool ret=true;
	b strSQL;
	int island = fNm.find("\\");
	while(island >= 0 ) {
		fNm = fNm.substr(0,island) + "/"+ fNm.substr(island+1,fNm.length()-island);
		island = fNm.find("\\");
	}
	strSQL="Attach  Database '" + fNm + "' As " + aNm + " ; "; 
	if ( sqlite3_exec( mDb, strSQL.c_str(), NULL, NULL, NULL ) == SQLITE_OK )
	{
		fprintf(logFile,"Table %s is attached as %s !\n SQLite Message : %s ", fNm.c_str(), aNm.c_str(), sqlite3_errmsg( mDb));
	} else {
		fprintf(logFile,"Table %s is not attached !\n SQLite Message : %s \n program exiting ", fNm.c_str(), sqlite3_errmsg( mDb));
		ret = false;
	}
return ret;
}

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.is_open();
}
bool checkExistence(const char *filename)
{
    ifstream f;
    f.open(filename);

    return f.is_open();
}
	template <typename a, typename q,typename o,  typename c, typename d,typename k> 
	c& dbkStOrdr(a& netdb , q& q1,o& o0, c& m1,d& logFile,k& key)
	{
//FLDPERD,FLDBEGT,FLDHDWAY,FLDPDLEN,COSTOPER
//SELECT FLDPERD, FLDBEGT, FLDHDWAY, FLDPDLEN, COSTOPER FROM pdhdwayIB ORDER BY FLDPERD, FLDBEGT
		int i=0, j=0;
		sqlite3_stmt *stmt = NULL;
		int rows=0;
		char *zErrMsg = NULL;
		string pkName="", sstopi="";
		int pkCount = 0;
		int fldNo = 0, col=0;
		int iFldVal=0, stopid=0,stopi=0;
		long lFldVal=0;
		double dFldVal=0,mcost=0;
		string svertid="", sstopid="", sFldVal="";
		typedef pair <k,o> pairKO; 
		m1.clear();

		if ( sqlite3_prepare_v2( netdb, q1.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );
		  fprintf(logFile,"\nSQLite error: %s\n\nSQL: %s \n", q1.c_str() ,  sqlite3_errmsg( netdb) );  //"SQLite error:  "<<q1<< " !"<<endl<<to_string<const char *>(sqlite3_errmsg( netdb)) <<endl;
		  return m1;
		}
		// query the stop detail data from the stop table

		while ( sqlite3_step( stmt ) == SQLITE_ROW )
		{
			rows++;
			// query the period detail data from the period table
			o o1;
			// (j==0) // k or as per query passed , must be number ( long )
			col=0;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal = sqlite3_column_int(stmt, col);
				o1.set_i(lFldVal);
			}
			// if (j==1) // StopId or per query passed,  must be number ( long)
			col=1;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal = sqlite3_column_int(stmt, col);
				o1.set_j(lFldVal);
			}
			// (j==2) // StOrdr or as per query passed, must be number(long)   
			col=2;
			if (sqlite3_column_bytes(stmt, col) !=0) {
				lFldVal = sqlite3_column_int(stmt, col);
				o1.set_k(lFldVal);
			} 
			if (key==1) {
				m1.insert(pairKO(o1.get_i(),o1)); // Stop id as map key (based on query)
			} else if (key==2) {
				m1.insert(pairKO(o1.get_j(),o1)); // Stop Order as map key (based on query)  
			} else if (key==3) {
				m1.insert(pairKO(o1.get_k(),o1)); // stop index K  as map key (based on query)
			}
		} // loop over rows
		sqlite3_finalize(stmt );
		return m1;
	}

	static inline void ReplaceAll2(std::string &str, const std::string& from, const std::string& to)
	{    size_t start_pos = 0;
		while((start_pos = str.find(from, start_pos)) != std::string::npos) {
			str.replace(start_pos, from.length(), to);
			start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
		}
	//	return str;
	}



