/* 
Stop Spacing Analysis and Optimization program that performs the Dynamic Programming Result Set for the 
periods in a Period Table that are included. 
If no scene in the form of a text or database file is given, a DP process is assumed and a full scene generation is performed.
Optimization procedure follows for each period in the included period list with individual result tables 
being created for each period.
If a scene file is given, then it is checked if it is a database. If so then it is checked for "%DPTrace" 
tables out of which are generated scene sets. Otherwise a set of lines with stops sets scene lists are assumed. 
The program then creates a single output database both with historic and alternative scene 
results and write them to a single database output file for all the scenario results regardless of the input 
types ( a SQLite database or a list of stops in text file is provided as input).  The output Result name is 
provided in the parameter file as a "ResultDB" variable. 
Added ability to distinguish fixed and eliminatable stops in the DP alternative
*/

#include <iostream>
#include <assert.h>
#include <string>
#include <string.h>
#include <strstream>
#include <stdlib.h>
#include <fstream>
#include <list> 
#include <iomanip>   // format manipulation
#include <algorithm>
#include <sstream>
#include <math.h>
#include <set>


#include "bustopsqlitev9sme.h"
#include "require.h"
#include "printSeq.h"
#include "StringConv.h"
#include "p2dp3d.h"
#include "vertexp.h"
#include "strip.h"
#include "stripkey.h"

#ifdef _MSC_VER
#define strcasecmp(a,b) stricmp(a,b)
#endif

	// container definition
	// mmapV1V2EgId - Vertex 1-2 -> Edge Id multi-map 
	// mmapEOIdEgId - Edge Object Id -> Edge Id multi-map 
	// mmaparced - Parcel ID -> Edge Id Multi-Map 
	// mapvert - Vertex id -> Vertex Object Map 
	// mmaped - edge id -> Edge Object multi-map
	// mmapStopEdge - Stop Id -> edge id multi-map
	// mmaparStop - Stop Assigned -> Parcel Object multi-map 
	// mmapved - Vertex id -> Edge Id - Multi-map
	 // input data field definition
	typedef map<int, string, less<int> > mapflds;
	mapflds :: iterator fld1map_Iter, fld2map_Iter;
	typedef pair < int, string >  fld_pair;
	mapflds mapfldshdr;
	mapflds mapfldsdet;
	mapflds::iterator mapfldsit;
	 
// global cost input data map
	 typedef map <long,globcost,less<long>> gcostmap;
     typedef pair <long, globcost> gcpair;
	 gcostmap gcmap;
	 gcostmap::iterator gcostmit;
// Period Headway input data map
	 typedef map <long,pdhdway,less<long>> phwmap;
     typedef pair <long, pdhdway> phwpair;
	 phwmap phmap;
	 phwmap::iterator phwmit;

// Period Headway input data map with period Key
	 typedef map <string,pdhdway,less<string>> phwsmap;
     typedef pair <string, pdhdway> phwspair;
	 phwsmap phsmap;
	 phwsmap::iterator phwsmit;
	 bool blnTripRun = false;
// Land use codes input data map
	 typedef map <string,lucodes> lucmap;
     typedef pair <string, lucodes> lucpair;
	 lucmap lumap;
	 lucmap::iterator lucmit;

// transit trip summary map sorted by Stop sequence
	 typedef map <long, stopKey, less<long>> maplngstopKey;
     typedef pair <long, stopKey> tskidpair;
	 maplngstopKey tskidmap, tskidpmap;
	 maplngstopKey::iterator tskidmit,tskidmit2;

// transit trip summary map sorted by Number of Stops 
	 typedef map <long, stopKey, greater<long>> maplngrkstop;
     typedef pair <long, stopKey> tseqidpair;
	 maplngrkstop tskgridmap, tskgridpmap;
	 maplngrkstop::iterator tskgridmit,tskgridmit2;

// transit trip summary map sorted by increasing trip start time sequence
	 typedef map <double, stopKey, less<double>> mapdblstopSeq;
     typedef pair <double, stopKey> tsKdblpair;
	 mapdblstopSeq tskdblmap, tskdblpmap;
	 mapdblstopSeq::iterator tskdblmit,tskdblmit2;

	 
	// transit stop data map sorted by Stop ID
	 typedef map <long, tstop, less<long>> maplngstop;
     typedef pair <long, tstop> tsidpair;
	 maplngstop tsidmap, tsidpmap;
	 maplngstop::iterator tsidmit,tsidmit2;
// transit stop data map sorted by increasing ride time
	 typedef map <double,tstop,less<double>> mapdblstop;
     typedef pair <double, tstop> tsrtpair;
	 mapdblstop tsrtmap;
	 mapdblstop::iterator tsrtmit,tsrtmit2;

// dp transit stop data map collection sorted by scenario ID
	 typedef map <long, dptStop, less<long>> maplngdpstop;
     typedef pair <long, dptStop> dptsidpair;
	 maplngdpstop dptsidmap, dptsidpmap;
	 maplngdpstop::iterator dptsidmit, dptsidmit2;

// trip dp transit stop data map collection sorted by trip ID
	 typedef map <long, tripDPs, less<long>> maplngtDPs;
     typedef pair <long, tripDPs> tripDPspair;
	 maplngtDPs tripDPsmap, tripDPsrmap;
	 maplngtDPs::iterator tDPsmit, tDPsmit2;

// dp transit stop data map collection sorted by dp key ID
	 typedef multimap <string, dptStop, less<string>> mmapstrdpstop;
     typedef pair <string, dptStop> strdptidpair;
	 mmapstrdpstop sdptsidpmmap, sdptsidprmmap;
	 mmapstrdpstop::iterator sdptsidmmit, sdptsidmmit2;

	 typedef map <string, dptStop, less<string>> mapstrdpstop;
	 mapstrdpstop sdptsidpmap, sdptsidprmap;
	 mapstrdpstop::iterator sdptsidmit, sdptsidmit2;


	 typedef multimap <string, tstop> mapsstop;
     typedef pair <string, tstop> strstopair;
	 mapsstop mapstrstop;
	 mapsstop::iterator msstopidmit, msstopidmit2;

//vertex map  object
	 typedef map<long,vertexp,less<long>> maplngvx;
// vertex map iterator objects
   maplngvx :: iterator vxmap_Iter, vomap_Iter;
   typedef pair < long, vertexp >  vx_pair;
   maplngvx mapvert,mapvert0,mapvert1,mapvert2,mapvert3;
   maplngvx mapvertor,mapvertor0,mapvertor1; // map vertex origin
   maplngvx scc; // strongly connected component
   maplngvx::iterator mapverit;

   typedef multimap<long,vertexp,less<long>> mmaplngvx;
// vertex multi map iterator objects
   mmaplngvx :: iterator vxmmap_Iter, vommap_Iter;
   mmaplngvx mmapvert,mmapvert0,mmapvert1,mmapvert2,mmapvert3;
   mmaplngvx::iterator mmapverit;

// Edge multi map object
   typedef multimap <long,edgev,less<long>> mmaplnged;
// edge multi-map iterator objects

   typedef   mmaplnged :: iterator edfmap_Iter, edrmap_Iter;
   typedef pair < long, edgev >  ed_Pair;
   mmaplnged mmaped,mmaped0,mmaped1,mapegor0,mapegor1;
   mmaplnged::iterator mmapedit;
   edfmap_Iter mmapedi;   //edge multimap iterator

// Edge map object
typedef map <long,edgev,less<long>> maplnged;
// edge map iterator objects
typedef   maplnged :: iterator efmap_Iter, ermap_Iter;
maplnged maped,maped0,maped1,maped2;
maplnged::iterator mapedit;
efmap_Iter mapedi;   //edge map iterator


//typedef deque<edgev> edgedeque;
// vertex origin map object

// MultiMap Edge-vertex map iterators object container definition
 typedef multimap <long, long> mmaplng;
 mmaplng mmapved;
 mmaplng mmapStopVx,mmapStopVx0,mmapVxStop,mmapVxStop0,mmapeidvid,mmapV1V2EOId; // Edge vertex map, dp stop vertex, dp vertex stop
 mmaplng mmapStopEdge,mmapEdgeStop, mmapV1V2, mmapV1V2EgId, mapEidOid,mmapEOIdEgId;
 typedef mmaplng :: iterator mumved_AIter, mumved_RIter;
 mumved_AIter mumavedit,mumESit,mumVSit;
 mumved_RIter mumrvedit;
 typedef pair < long, long >  lng_Pair;
// sorted multi map vertices for the voronoi algorithm
 typedef multimap <double, long> mmapdblng;
 mmapdblng mmapSver,mmapSver2,mmapVx0,mmapVx1;	 
 mmapdblng :: iterator mmapSverit;
 typedef pair < double, long >  dblng_Pair;

// multi map vertex ids and dbl values for the voronoi algorithm
 typedef multimap <long,double > mmaplngdbl;
 mmaplngdbl mmapVerlng, mmapVerlng2, mmapVx0R, mmapVx1R;	 
 mmaplngdbl :: iterator mmapSveRit;
 typedef pair < long,double >  lngdbl_Pair;
 lngdbl_Pair lngdblPair;

 // map dp ids and dbl values for the DP algorithm
 typedef map <long,double, less<long>> maplngdbl;
 maplngdbl mapDPRes;	 
 maplngdbl :: iterator mapDPRit;

  // sorted multi map predecessor & successor for j & k of the DP algorithm
 mmaplng mmapDPredK, mmapDPredJ;

 // sorted multi map stop cost with ID for the DP algorithm
 mmapdblng mmapDPIC,mmapDPRes;	 
 mmapdblng :: iterator mmapDPICit,mmapDPRit;
 typedef pair < double, long >  dblng_Pair;

 // multimap dp ij (stage j state i) values with decision variable k for the DP algorithm
 typedef multimap < string, long> mmapstrlng;
 mmapstrlng mmapIJDPK,mmapIJDPResICost;	 
 mmapstrlng :: iterator mmapIJDPRit;
mmapstrlng :: iterator mmapIJDPKit;

 typedef pair < string,long >  strlng_Pair;

 // map dp string keys and dbl values for the DP algorithm
 typedef multimap <string,double> mmapstrdbl;
 mmapstrdbl mmapIJDPRes,mmapIJKDPICost;	 
 mmapstrdbl :: iterator mmapIJKDPRit;
 mmapstrdbl :: iterator mmapKLMDPRit;
 typedef pair < string,double >  strdbl_Pair;
 
 // map dp Immediate Cost to string key values for the DP algorithm
 typedef multimap <double,string, less<double>> mmapdblstr;
 mmapdblstr mmapICostDPRes,mmapICostDPIJ;	 
 mmapdblstr :: iterator mmapICDPIJRit;
 typedef pair < double,string >  dblstr_Pair;

 // map dp decision k and (stage j + state i) string key for the DP algorithm optimal path
 typedef multimap <long, string, less<long>> mmaplngstr;
 mmaplngstr mmapStopDPRes;	 
 mmaplngstr :: iterator mmapStopDPRit;
 typedef pair < long, string >  lngstr_Pair;

 // map dp string keys and dbl values for the DP algorithm
 typedef map <string,double> mapstrdbl;
 mapstrdbl mapStrDPRes;	 
 mapstrdbl :: iterator mapStrDPRit;
 typedef pair < string,double >  strdbl_Pair;

 // map dp stop id and string key values for the DP algorithm optimal path
 typedef map <long, string, less<long>> maplngstr;
 maplngstr mapStopDPRes, mapMsg;	 
 maplngstr :: iterator mapStopDPRit, mapMsgIt;
 typedef pair < long, string >  lngstr_Pair;


 // parcel multi map for parcel and edge  and Parcel and Stop connection
 typedef multimap <long, parcel, less<long>> mmaplngpar;
 mmaplngpar mmaparced,mmaparcel,mmaparced0,mmaparced1;	 
 mmaplngpar mmaparStop,mmaparStop0,mmaparStop1;	 
 typedef mmaplngpar :: iterator mmapP_AIter;
 mmapP_AIter mmAPEi; 
 typedef pair < long, parcel >  PE_Pair;
// map parcel
  typedef map <long, parcel, less<long>> maplngpar;
 maplngpar maparcid,maparcel,maparcid0,maparcid1;	 
 typedef maplngpar :: iterator mapParIt;
 mapParIt mPi; 

//long globcost::create = 0, globcost::assign = 0,
//     globcost::copycons = 0, globcost::destroy = 0;
long globcost::id = 0;
// output file object for write 
ofstream outfile; // output file stream
ifstream infile; // input file stream
ofstream outedgefile;
ofstream outedgebdyfile;
ofstream outedgeobjectfile;
ofstream outvertfile;
ofstream outparcfile;


void fileName (char  [],char [],char  [], int onoff=-1,int alt=-1);
void fileName (char  [], string ,char  [], int onoff=-1,int alt=-1);

map<int , string> readatahdr(string& rec1, char *seps = "\t"); 
map<int , string> recoread(string& rec1, char *seps = "\t"); 
vertexp& readvertex1(string& rec1, vertexp& vx, map<int , string>::iterator maphdrit, char *seps = "\t");
edgev& readedge1(string& rec1, edgev& evx, map<int , string>::iterator maphdrit,globcost& gc , char *seps = "\t"); 
parcel& readparc1(string& rec1, parcel& parc, map<int , string>::iterator maphdrit, char *seps = "\t"); 
//tstop& readtstop1(string& rec1, tstop& pstop, map<int , string>::iterator maphdrit, char *seps = "\t"); 
globcost& readgcost1(string& rec1, globcost& gc1, map<int , string>::iterator maphdrit, char *seps = "\t"); 
pdhdway& readpdhdway1(string& rec1, pdhdway& pdh1, map<int , string>::iterator maphdrit, char *seps = "\t"); 
lucodes& readlucodes1x(string& rec1, lucodes& luc1, map<int , string>::iterator maphdrit, char *seps = "\t");
template <typename s> s& rideTimetoEnd(s& tsrtmap);
list<string> loadFields(sqlite3* netdb ,string mTableName,bool isQuery=true);

class ParcelAccum;
class ParcelOnAccum;
class ParcelOffAccum;
class VertAccum;
class VxVorAccum;
class EdgeAccum;
maplngvx& vxv (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,
		  mmaplng& mmapved,mmaplng& mmapEdgeStop,short onoff); //vertex voronoi
void edv (maplngvx& mapvert, mmaplnged& mmaped,mmaplng& mmapved,
		  mmaplng& mmapEdgeStop,short onoff); //edge voronoi
void parcstopwalk(mmaplng& mmapEOIdEgId,mmaplngpar& mmaparced,maplngvx& mapvert, 
			  mmaplnged& mmaped,mmaplngpar& mmaparStop,short onoff,bool blnHist=true);

template <typename l, typename p, typename v, typename e, typename o>
p& parcstopwalknew(l& mmapEOIdEgId, p& mmaparced, v& mapvert, e& mmaped,
				  p& mmaparStop,short onoff,bool blnHist, o& logfile );
//template <typename o>
//void parcstopwalknew(mmaplng& mmapEOIdEgId,mmaplngpar& mmaparced,maplngvx& mapvert, 
//			  mmaplnged& mmaped,mmaplngpar& mmaparStop,short onoff,bool blnHist, o& logfile);
void parcstopwalkx(mmaplng& mmapEOIdEgId,mmaplngpar& mmaparced,maplngvx& mapvert, 
			  mmaplnged& mmaped,mmaplngpar& mmaparStop,bool blnHist=true);

maplngstop&  stopspacing_analysis (mmaplng& mmapV1V2EOId,mmaplng& mmapV1V2EgId,
					mmaplngpar& mmaparced,maplngvx& mapvert,maplnged& maped,
					mmaplngpar& mmaparStop, mmapdblng& mmapSver,mmaplng& mmapVxStop,
					mmaplng& mmapStopVx,mmapdblng& mmapVx0, mmapdblng& mmapVx1,
					mmaplngdbl& mmapVx0R,mmaplngdbl& mmapVx1R, bool blnHist,bool blnEuclid,
					string strEdgeBaseName, string strVertexBaseName, string strParcelBaseName,string strStopBaseName,
					globcost& gc1, pdhdway phway, lucodes luci, maplngstop& stops, int q, ofstream& logfile,
					inputfilelist& inpList, stopKey& kStop, sqlite3* db=NULL, string sADB="NetDB", int xi=9);

template <typename T>
deque<T*> pagwalkstop (maplngstop& tsidmap, mmaplngpar& mmaparStop,T* pa);

deque<ParcelOnAccum> pagwalkstopon (maplngstop& tsidmap, mmaplngpar& mmaparStop,lucodes* luci,
								globcost* gc1,short onoff,ParcelOnAccum& pa,short blnHist=true, ostream& logfile=cout);
deque<ParcelOffAccum> pagwalkstopoff (maplngstop& tsidmap, mmaplngpar& mmaparStop,lucodes* luci,
								globcost* gc1,short onoff,ParcelOffAccum& pa,short blnHist=true, ostream& logfile=cout );


// Stop Index and StOrdr Map values 
	typedef map <long,long,less<long>> mkOrdr;
	typedef pair <long, long> pairkOrdr;
	mkOrdr mapkOrdr;
	mkOrdr::iterator mapkOrIt;

// Stop Index and StOrdr Map values 
	typedef map <long,ltStop,less<long>> mkLtStop;
	typedef pair <long, ltStop> pairkLtStp;
	mkLtStop mapxLtStp, mapkLtStp, mapstOrLtStp, mapstIdLtStp;
	mkLtStop::iterator mapkLtIt,mapoLtIt,mapiLtIt;

// Stop Index and StOrdr Map values 
	typedef map <string,mkLtStop,less<string>> mkLtString;
	typedef pair <string, mkLtStop> pairmkLtString;
//	typedef map <string,std::map <long,ltStop>> mkLtString;
//	typedef pair <string, std::map <long,ltStop>> pairmkLtString;
//	std::map <string,std::map <long,ltStop>>  mapksTimeLt;
	mkLtString mapTripLt , mapksTimeLt;
	mkLtString::iterator mapkLSIt,mapstrLSIt,mapLtIt;

// Stop Index and StOrdr Map values 
	typedef map <long ,mkLtStop,less<long>> mkLtLg;
	typedef pair <long, mkLtStop> pairmkLtLg;
	mkLtLg mapTripLtLg , mapksTimeLtLg;
	mkLtLg::iterator mapkLSItLg, mapstrLSItLg, mapLtItLg;

// Geometry function headers
double ds( Point P, Point Q);
double d2( Point P, Point Q);
int closest2D_Point_to_Line ( Point P[], int n, Line L);
float dist_Point_to_Line( Point P, Line L);
float dist_Point_to_Segment( Point P, Line S);



void main(int argc , char* argv[] )
{

	// SQLite variables
	sqlite3 *pTransInpDb;
	sqlite3 *pSceneDb;
	sqlite3 *pTransOutDb;
	FILE* txtStream;
	string fileName1="", tblStopBuf="", tblStopBufx2="";
	string strSQLInserTbl="" ,strDPSQLCreateTbl="" ,strDPSQLInserTblDef= "", strDPSQLInserTbl="" , strDPTblName="";
	char *zErrMsg = 0;
	int rc=0;
	int trafficStress=0, xi = 9;

	sqlite3_stmt *stmtStopTablei = NULL;
	sqlite3_stmt *stmtTripTablei = NULL;
	sqlite3_stmt *stmtDmndTablei = NULL;
	sqlite3_stmt *stmtstrSQLTablei= NULL;
	sqlite3_stmt *stmtDPTables=NULL;
	char *errMsg = NULL;


	//initialize the edges
    edgev ev;
    edgev* evp;
    evp = &ev;
  //initialize the vertices
    vertexp vx;
    vertexp* pVx1=&vx;
	vertexp* pVx2=&vx;
    vertexp* vxp=&vx;

	//initialize the parcels
    parcel par1;
    parcel* parp;
    parp = &par1;

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
	// dpstop initialization
	dptStop dpstop;
    dptStop* pdpstop;
	pdpstop = &dpstop;

	// trip dpstop initialization that includes walk / ride cost coefficients, trip id, dp key
	tripDPs tDPstp;
	tripDPs* ptDPstp;
	ptDPstp = &tDPstp;

	//initialize the Land use code data object
    lucodes luc1;
    lucodes* luci;
    luci = &luc1;

	//initialize the Period Headway data object
    pdhdway phw1;
    pdhdway* phway;
    phway = &phw1;

	//initialize the Global Cost data Object
    globcost gc1;
    globcost* gcost;
    gcost = &gc1;

	//initialize the trip key object
    kTrip kTrip1;
    kTrip* pkTrip;
    pkTrip = &kTrip1;

	//initialize the Key Stop data object that includes route, schlName,  timePeriod, tripStartTime, 
	//tripTime, tripKey , tripNumber, stopCnt - number of stops, timeNexTrip, histOnSum, histOffSum 
    stopKey kStop1;
    stopKey* pkStop;
    pkStop = &kStop1;

	//initialize the input file list data Object
    inputfilelist InpFiles;
    inputfilelist* pInFiles;
    pInFiles = &InpFiles;

    ltStop ltSt,ltiSt,ltoSt,ltkSt;
    ltStop* pltSt;
	ltStop* pltiSt;
	ltStop* pltoSt;
	ltStop* pltkSt;
    pltSt = &ltSt;
    pltiSt = &ltiSt;
    pltoSt = &ltoSt;
    pltkSt = &ltkSt;
	list<string> lstPdNm;
	list<string> lstblNm;
	list<string> lstTripId;
	list<string>::iterator lstit;
	set<int> sElim;
	set<int>::iterator sElit;



// input file object and read 
//ifstream einfile; // edge file
//ifstream vinfile; // vertex file 
	 ifstream parfile; // parcel file
	 ifstream tstopfile; // transit stop file 
	 ifstream gcostfile; // Global Parameters file 
	 ifstream pdhdwayfile; // Period Parameters file 
	 ifstream lucifile; // parcel Land Use file
	 ifstream inputlistfile; // input file with list of tables , parameters etc
	 ifstream sceneifile; // scenario input file
	 ofstream logfile; // log file
	 ofstream stopout;

	 stringstream strout; 

    string rec1="", str1="", strSQL="", strDirn="", strTmPd="";
	double xcost=0,ycost=0,ecost=0, tt=0;
    int il=0, origin=0;
//    char seps[]   = ",\t\n";
//	char *seps = ",";
    char xs[]   = "\t";
    int ib=10, intSRID=0; //base for the string to number conversion
    int ne=0,nv=0, col=0; // ne - no of edges, nv - no of vertices
    long i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0, onoff=0,vid=0,v1=0,v2=0,eid=0,ip=0,pid=0,id=0,key=0; // vid - vertex id(v1 & v2 - vertex 1 & 2 of edge), eid - edge id, pid - parcel id
	int M=3, D=0;
	short lbl=0; // whether an object is labeled or not unlabeled=0, labeled=0.
	bool blnHistoric=true, qNet=true, blnS=true,blnCreate=false, blnInsert=false;
	bool blnPd=true, blnScene=false,blnSceneDB=false, blnPdNm=false,blnParcel=false;
	string strSQLstopTable="",strSQLTabList="", strPdHdwySQLTable = "", strSQLTable = "", keyFld="",tblDP="";	// sql string
	string txtStreamName="", logfileName="",txtStdoutName="",txtStderrName="", basefileName="";
	string strStopBaseName="", strEdgeBaseName="",strVertexBaseName="", strParcelBaseName="",tripNum="";
	string strStopFileName="", strEdgeFileName="",strVertexFileName="", strParcelFileName="";
	string strPdHdwayName="", strGlobCostName="",strLUCodeName="", networkDbName="", strcTransOutDbName="";
	string strInpListName="", strOutFileName="",strEdgeBinName="",strVertexBinName="",strParcelBinName="";
	string tblTrip="", strADB = "netDB", strsceneFileName="", strWR="",pdNm="";
	//txtStream.setRealNumberNotation(QTextStream::FixedNotation);	
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );


	i=0;
	mapMsg.insert(lngstr_Pair( i,"Walk from Stop to Parcel (offs) "));
	i=1;
	mapMsg.insert(lngstr_Pair(i,"Walk from Parcel to Stop (ons) "));
// open the input list of files 
	//strcpy(inputfilenames, argv[1]);
	if (!argv[1]) {
		cout << "Name of Input File List  (max length 200 chars ) : ";
		cin >> strInpListName;
	}
	strInpListName = string(argv[1]);
	

	inputlistfile.open(strInpListName.c_str(), ios::in );
	while (!inputlistfile) 
	{ // w1
    	cout << "Input List file not found. Please re-enter again : ";
	    cin >> strInpListName;
		inputlistfile.clear();
    	inputlistfile.open(strInpListName.c_str(), ios::in );

	   if (++j>3) 
	   { // i1
		   cout<<"More than three trials opening the input file name!"<<endl;
		   cout<<" Press any key to exit!"<<endl;
		   cin>>str1;
		   exit (0);
	   } // 1i
	} // 1w


	while(!inputlistfile.eof())
	{ // w1 call the reading routine
		getline(inputlistfile,rec1);
		if (rec1.length() >0)
		{ // i1
			InpFiles = readInputFiles(rec1,*pInFiles, xs);
			cout <<endl<<rec1<<endl;
		} // 1i
	} // 1w
	bool blnEuclid =InpFiles.get_euclid();

	if (inputlistfile.is_open())
	{ // i1
		inputlistfile.close();
	} // 1i
	inputlistfile.clear();

	if (InpFiles.get_tblstop().length()>1)
	{
		strStopBaseName = InpFiles.get_filestop().substr(0,InpFiles.get_filestop().find_last_of(".")) ;
		if (strStopBaseName.length() > 0 )
		{
			logfileName = (strStopBaseName + "_opt.log"); 
			txtStreamName = (strStopBaseName + "_opt.out");
			txtStdoutName = (strStopBaseName + "_std.out");
			txtStderrName = (strStopBaseName + "_err.out");
		} else
		{
			logfileName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_opt.log" );
			txtStreamName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_opt.out" );
			txtStdoutName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_std.out" );
			txtStderrName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_err.out" );
		}
	} else
	{
		logfileName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_opt.log" );
		txtStreamName =  ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_opt.out" );
		txtStdoutName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_std.out" );
		txtStderrName = ("c:/temp/StOptimize_" + to_string<int>(rand()) + "_err.out" );
	}
	fopen_s(&txtStream,txtStreamName.c_str(), "w");
	if (txtStream==NULL)
	{
		fprintf(stderr, "File %s is not a valid path \n", txtStreamName.c_str());
		exit(1);
	}

// reassign stdout to files
	FILE *outstream;

   rc = freopen_s( &outstream,txtStdoutName.c_str() , "w", stdout );

   if( rc != 0 )
      fprintf( stdout, "error on stdout assignment using s_freopen for file name %s .\n",txtStdoutName.c_str() );

// reassign stderr to files
	FILE *errstream;

   rc = freopen_s( &errstream,txtStderrName.c_str() , "w", stderr );

   if( rc != 0 )
      fprintf( stdout, "error on stdout assignment using s_freopen %s .\n",txtStderrName.c_str() );


	spatialite_init(0);
	intSRID = InpFiles.get_srid();
	blnParcel = InpFiles.get_parcel();
	networkDbName =	 InpFiles.get_inputDB() ;	 //to_string<char*> (networkDatabaseName);
	//const char* networkDatabaseName = new char[networkDbName.length() + 1];
	//networkDatabaseName = networkDbName.c_str();
	rc = sqlite3_open_v2( networkDbName.c_str(), &pTransInpDb, SQLITE_OPEN_READWRITE, NULL );
	if( rc != SQLITE_OK){  // i1
			fprintf(stderr, "Can't open network database: %s error %s .\n", networkDbName.c_str(), sqlite3_errmsg(pTransInpDb));
			fprintf(txtStream, "Can't open network database: %s error %s .\n", networkDbName.c_str(), sqlite3_errmsg(pTransInpDb));
			sqlite3_close(pTransInpDb);
			exit(1);
	} // 1i

	strSQLTabList = " select type , name , tbl_name , rootpage , sql from sqlite_master " 
			" where tbl_name like '%' and type like 'table' "
			" ORDER BY tbl_name; " ;

	if ( sqlite3_prepare_v2( pTransInpDb, strSQLTabList.c_str(), -1, &stmtstrSQLTablei, NULL ) != SQLITE_OK )
	{ // i1
		// some error occurred
		fprintf(txtStream,"SQLite prepare error: \tSQL %s \n error: %s .\n", strSQLTabList.c_str(), sqlite3_errmsg( pTransInpDb));
		fprintf(stderr,"SQLite prepare error: \tSQL %s \n error: %s .\n", strSQLTabList.c_str(), sqlite3_errmsg( pTransInpDb));
		exit(1);
	} // 1i
	// open/create the hull output sqlite/spatialite database


// globcostiname - SQL statement to query the global cost detail table data
//SELECT  COSTWALK, COSTRIDE, UNITONTM, UNITOFFTM, MAXWLKDIST, PROPENSITY, FILESTEM, NOPERIODS, 
//WALKSPD, FILEPATH, MAXSKIP, DPDIMENSION , unitConv FROM GlobCostPXIBD5
	//strcpy(globcostiname, InpFiles.get_tblgcost().c_str());

		// strSQLTable - SQL statement to query the global cost table data
	strGlobCostName = InpFiles.get_tblgcost();
	strSQLTable = "SELECT COSTWALK, COSTRIDE, UNITONTM, UNITOFFTM, MAXWLKDIST, PROPENSITY, "
		" FILESTEM, NOPERIODS, WALKSPD, FILEPATH, MAXSKIP, DPDIMENSION, unitConv "
		" FROM " + strGlobCostName;
		gcmap = globalCost(pTransInpDb, strSQLTable,*gcost,gcmap,txtStream);

	// open the log file 

	logfile.open(logfileName.c_str(), ios::out | ios::trunc);



	M = gcost->get_maxskip() + 1;

	D = gcost->get_dpdimension();


	gcost->show_globcosthdr(logfile);
	gcost->show_globcost(logfile);

	str1 = datetimeStamp(logfile);
	cout<<str1<<endl;
	cout <<endl<< D <<" Dimension Run started"<<endl;
	logfile << D <<" Dimension Run started "<<endl;
	logfile<<"Input SQLite Database "<<"\t"<<InpFiles.get_inputDB()<<endl;
	logfile<<"Result output SQLite Database "<<"\t"<<InpFiles.get_resultDB()<<endl;
	logfile<<"Global Cost Table "<<"\t"<<InpFiles.get_tblgcost()<<endl;
	logfile<<"Period Table "<<"\t"<<InpFiles.get_tblperiod()<<endl;
	logfile<<"Land Use Table "<<"\t"<<InpFiles.get_tbllanduse()<<endl;
	logfile<<"Scenario File "<<"\t"<<InpFiles.get_tblscene()<<endl;
	logfile<<"Stop Table "<<"\t"<<InpFiles.get_tblstop()<<endl;
	logfile<<"Parcel Table "<<"\t"<<InpFiles.get_tblparcel()<<endl;
	logfile<<"Edge Table "<<"\t"<<InpFiles.get_tbledge()<<endl;
	logfile<<"Vertex Table "<<"\t"<<InpFiles.get_tblvertex()<<endl;
	logfile<<"Global Cost File  "<<"\t"<<InpFiles.get_filegcost()<<endl;
	logfile<<"Period File  "<<"\t"<<InpFiles.get_fileperiod()<<endl;
	logfile<<"Land Use File  "<<"\t"<<InpFiles.get_filelanduse()<<endl;
	logfile<<"Scenario File "<<"\t"<<InpFiles.get_filescene()<<endl;
	logfile<<"Stop File  "<<"\t"<<InpFiles.get_filestop()<<endl;
	logfile<<"Parcel File  "<<"\t"<<InpFiles.get_fileparcel()<<endl;
	logfile<<"Edge File  "<<"\t"<<InpFiles.get_fileedge()<<endl;
	logfile<<"Vertex File  "<<"\t"<<InpFiles.get_filevertex()<<endl;
	logfile<<"Run Historic Scenario First "<<"\t"<<InpFiles.get_historic()<<endl;
	logfile<<"Use Euclidean Geometry "<<"\t"<<blnEuclid<<endl;

	// check to see if this is a scene data and if so find if the list is from a DP result datbase 
	strsceneFileName = InpFiles.get_filescene();
	rc = sqlite3_open_v2( strsceneFileName.c_str(), &pSceneDb, SQLITE_OPEN_READONLY , NULL );
	if (rc == SQLITE_OK) {  blnSceneDB = blnScene = true; }
	if( !blnScene ){
		blnScene = checkExistence(strsceneFileName.c_str());
	}

	if (blnScene) {
		// The file exists, close it and check to see if it is possible to open it as a sqlite database
		//if (sceneifile.is_open()) {
		//	sceneifile.close();
		//} else {
		//	sqlite3_close(pSceneDb);
		//}
		// this maybe a list of scenarios as a text or a list of result tables from the  DBS run of DP sqlite DB
		// try opening it as a sqlite table

			// open Scene DB
			rc = sqlite3_open_v2( strsceneFileName.c_str(), &pSceneDb, SQLITE_OPEN_READONLY , NULL );
			if( rc != SQLITE_OK){
				fprintf( stderr,"\tSQLite error opening spatialite transit optimization output database : %s ; \n error %s \n", strsceneFileName.c_str(), sqlite3_errmsg( pSceneDb));
				fprintf( txtStream,"\tSQLite error opening spatialite transit optimization output database : %s ;\n error %s \n", strsceneFileName.c_str(), sqlite3_errmsg( pSceneDb));
			} 

			// strWR and rec1 are to be used to find trip and period id 
			//strWR = "D" + to_string<short>(D);
			strWR = ("W" + to_string<int> ((int) (gcost->get_walkcost()*10)));
			strWR.append("R" + to_string<int> ((int) (gcost->get_ridecost()*10)));
			// rec1 to be used to find trip and period id 
			rec1 = InpFiles.get_tblstop() + "D" + to_string<short> (D) + "t" ;

			//strDPTblName = inParams.get_tblstop() + "D" + to_string<short> (D) + "t" + to_string<long> (kStop1.tripId())+ kStop1.tripPeriod() + "W" + to_string<int> ((int) (gcost->get_walkcost()*10))+ "R" + to_string<int> ((int) (gcost->get_ridecost()*10)); 
			strDPTblName = rec1 + "%" + strWR + "%DPTrace%"; 
			replace(strDPTblName.begin(),strDPTblName.end(),' ','_');
			ReplaceAll2(strDPTblName,"__","_");

			// loop over all the tables in the DPResult database and read the individual trip/Period dp data 
			str1= " select tbl_name ,  name  , rootpage  from sqlite_master " 
			" where tbl_name like '" + strDPTblName +"' and type like 'table' "
			" ORDER BY tbl_name; " ;
			if ( sqlite3_prepare_v2( pSceneDb, str1.c_str(), -1, &stmtDPTables, NULL ) != SQLITE_OK )
			{
			  // some error occurred
			  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", str1.c_str() ,  sqlite3_errmsg(pSceneDb) );
			  fprintf(txtStream,"\nSQLite error: %s\n\nSQL: %s \n", str1.c_str() ,  sqlite3_errmsg(pSceneDb) );  //"SQLite error:  "<<q1<< " !"<<endl<<to_string<const char *>(sqlite3_errmsg( netdb)) <<endl;
			  //exit (-2);
			}
			j = 0;
			int jMx =0;
			// find the maximum number of stop sets among the dp result tables to use the basis for the optimization
			while ( sqlite3_step( stmtDPTables ) == SQLITE_ROW )
			{
				// get the table name and query the DP table 
				col=0;
				if (sqlite3_column_bytes(stmtDPTables, col) !=0) {
					tblDP  = to_string<const unsigned char * > (sqlite3_column_text(stmtDPTables, col));
					// extract the period name from "RTD_Trips_Rte20EBECPDXD5t8379662Midday_59W112R281_MPDPTrace" and store it in lstPdNm 
					ne = tblDP.find(strWR);
					keyFld = tblDP.substr((tblDP.find(rec1)+ rec1.size()),(tblDP.find(strWR)-tblDP.find(rec1)-rec1.size()));
					tripNum = keyFld.substr(0,7);
					pdNm = keyFld.substr(7,keyFld.size()-7);
					// find the DP K indices of the stOrdr key and store it into a map 
					str1 =	"SELECT distinct Stopid , StOrdr, K   "
						"FROM " + tblDP + " "   
						"ORDER BY K, StOrdr ";
					// make the multi-map for the index k and the StOrdr 
					key = 3;
					mapkLtStp = dbkStOrdr(pSceneDb, str1,*pltSt,mapkLtStp,txtStream,key);
					j = mapkLtStp.size();
					if (j>0) {
						//mapksTimeLtLg.insert(pairmkLtLg(from_string<long>(tripNum),mapkLtStp));
						if (j>jMx) {
							jMx = j;
							tblTrip = tblDP;
							lstblNm.push_front(tblDP);
							lstPdNm.push_front(pdNm);
							lstTripId.push_front(tripNum);
						} else {
							lstblNm.push_back(tblDP);
							lstPdNm.push_back(pdNm);
							lstTripId.push_back(tripNum);
						}
						blnSceneDB = true;
					} else {
							cout <<"No DP scenario data read from database "<<tblDP<<endl<<"Program ending!";
							logfile <<"No DP scenario data read from database "<<tblDP<<endl<<"Program ending!";
							fprintf(stderr, "\nNo DP scenario data read from file %s\n Exit Program! \n",  tblDP.c_str() );
							fprintf(txtStream, "\nNo DP scenario data read from file %s\n Exit Program! \n",  tblDP.c_str() );
							// exit (-1);
					}


					} // a table name is found
				} // read next row from tblDP
				// 

	}


 
	// tstopiname - intermediate results for stop detail table data
//	strcpy(tstopiname, InpFiles.get_filestop().c_str());
//	strcpy(edgeinfname,InpFiles.get_fileedge().c_str());
//	strcpy(vertinfname,InpFiles.get_filevertex().c_str());
//	strcpy(pariname,InpFiles.get_fileparcel().c_str());

	if (InpFiles.get_fileedge().length()>1)
	{
		if (InpFiles.get_fileedge().find_last_of(".")>0) {
			strEdgeBaseName = InpFiles.get_fileedge().substr(0,InpFiles.get_fileedge().find_last_of(".")) ;
		} else {
			strEdgeBaseName = InpFiles.get_fileedge() ;
		}
	} else
	{
		strEdgeBaseName = ("c:/temp/edgefile_" );
	}

	if (InpFiles.get_filevertex().length()>1)
	{
		if (InpFiles.get_filevertex().find_last_of(".")>0) {
			strVertexBaseName = InpFiles.get_filevertex().substr(0,InpFiles.get_filevertex().find_last_of(".")) ;
		} else {
			strVertexBaseName = InpFiles.get_filevertex() ;
		}
	} else
	{
		strVertexBaseName = ("c:/temp/vertexfile_" );
	}

	if (InpFiles.get_fileparcel().length()>1)
	{
		if (InpFiles.get_fileparcel().find_last_of(".")>0) {
			strParcelBaseName = InpFiles.get_fileparcel().substr(0,InpFiles.get_fileparcel().find_last_of(".")) ;
		} else {
			strParcelBaseName = InpFiles.get_fileparcel() ;
		}
	} else
	{
		strParcelBaseName = ("c:/temp/parcelfile_" );
	}

	// get if historic is to be run first before alternative analysis
	blnHistoric = InpFiles.get_historic() ;

	/*	cout << "Period Parameter table file name max length 200 chars (Enter to use default) : ";
		cin >> globcostiname;
	*/

	// pdhdwayiname - SQL statement to query the period table data
	//strcpy(pdhdwayiname,InpFiles.get_fileperiod().c_str());
	strPdHdwayName = InpFiles.get_fileperiod();
	/*	cout << "Parcel Land Use Code file name max length 200 chars (Enter to use default) : ";
		cin >> luciname;
	*/
	// luciname - sql query for landuse data
		//strcpy(luciname,InpFiles.get_tbllanduse().c_str());
	strLUCodeName = InpFiles.get_tbllanduse();

	// strSQLTable - SQL statement to query the land use table data
	//	strSQLTable =  " SELECT Ptype, Desc, LU , KeyFld, OffCoef, OnCoef  "
	//		" FROM  " + strLUCodeName + " "
	//		" Order by Ptype ";


		//lumap = readLUCodesTable(pTransInpDb, strSQLTable,str1,*luci,lumap,txtStream,strTmPd);
		blnPd = InpFiles.get_pdRun();

		// read period headway data table 
		// Find the maximmum number of stops per route and use that as the historic stop set
		strPdHdwayName=InpFiles.get_tblperiod();
		j=0;
		pdNm = "";
		if (blnPd) {
		// strSQLTable - SQL statement to query the period table data
			if (blnScene) {
				if (lstPdNm.size()>0) {
					blnPdNm = true;
					pdNm = "  TimePeriod in ( " ;
					for (lstit = lstPdNm.begin();lstit !=lstPdNm.end();++lstit)
					{
						j++;
						str1 = *lstit;
						replace(str1.begin(),str1.end(),'_',' ');
						if (j< lstPdNm.size()) {
						// get the period name to be used to get the stop set 
							pdNm  = pdNm + ( "\"" + str1 + "\",") ;
						} else {
							pdNm  = pdNm+ ( "\"" + str1 + "\")") ;
						}
					} // blnPd is true
				} else  {  // no time periods were detected the length of the Pd list is zero
					blnPdNm = false;
					pdNm = " TimePeriod like '%' " ;
				}
			} else { // no db scene file found 
				pdNm = " TimePeriod like '%' " ;
			}
			strPdHdwySQLTable = "SELECT FLDPERD, FLDBEGT, FLDHDWAY, FLDPDLEN, COSTOPER, TimePeriod, " 
						" numTrips, Include "
						" FROM  " + strPdHdwayName + " Where " +  pdNm + " "
						" ORDER BY FLDPERD, FLDBEGT";
			strSQLTable = "	select RteName,schlName,timePeriod , min(TripSTime) TripSTime,TripId,tripNumber, "
					" count(stop_Id) NumStops,  sum(Ons) SumHOns , sum(Offs) SumHOffs, dirName from " + InpFiles.get_tbltrip() + " "
					" where RteName like '" + InpFiles.get_rtename() +  "' and hist <> 0 and "
					" dirName like '%" + InpFiles.get_dirname() +  "%' and "  +  pdNm + " "  // " TimePeriod like '" + lstPdNm.front() + "' "
					" group by RteName,timePeriod "
					" order by tripNumber , RteName,count(stop_Id) , dirName ;";
		} else {
			strPdHdwySQLTable = "SELECT FLDPERD, FLDBEGT, FLDHDWAY, FLDPDLEN, COSTOPER, pdKey, numTrips, Include "
						" FROM  " + strPdHdwayName + " " 
						" ORDER BY FLDPERD, FLDBEGT";
			strSQLTable = "	select RteName,schlName, timePeriod,TripSTime,TripId,tripNumber, "
					" count(stop_Id) NumStops, sum(Ons) SumHOns , sum(Offs) SumHOffs, dirName from " + InpFiles.get_tbltrip() + " "
					" where RteName like '" + InpFiles.get_rtename() +  "' and hist >= 0  and "
					" dirName like '%" + InpFiles.get_dirname() +  "%'"
					" group by RteName,TripSTime,TripId,tripNumber "
					" order by tripNumber , RteName,count(stop_Id) , dirName ;";
		}
		// get the period headway list in the period table
		phmap = pdHdwy(pTransInpDb, strPdHdwySQLTable,*phway,phmap,txtStream);
		phsmap = remapdHdwy2pdKey(phmap,phsmap,*phway,str1,blnPd);
		// from the trip table  get the list of stops grouped by RteName,TripSTime,Tripid,tripNumber  
		tskidmap.clear();
		tskidmap = stopTableCount(pTransInpDb, strSQLTable,kStop1,tskidmap,txtStream);
		
		//tskdblmap = remap2TripTime(tskidmap,tskdblmap,kStop1,tt);
		fprintf( stderr,"\n\tStop Trip Summary Size : %d \t for Trip/Period %s and trip Number %s from table %s ; Route %s ;  and Dirn %s \n SQL : %s  \n",tskidmap.size(), kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str(),strSQLTable.c_str());
		fprintf( txtStream,"\n\tStop Trip Summary size : %d \t for  Trip/Period %s and trip Number %s from table %s ; Route %s ;  and Dirn %s \n SQL : %s \n", tskidmap.size(), kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str(),strSQLTable.c_str());
		if (tskidmap.size() == 0) {
			//logfile<<"No Stop Table summary data found for table "<<InpFiles.get_tbltrip()<< " and \n SQL : "<<"\t"<<strSQLTable<<endl;
			return;
		}
		kStop1.serializetexthdr(logfile);
		writeTextObjectData(tskidmap,kStop1,ip,logfile,""); 
		//tskgridmap = remap2TripCount(tskidmap,tskgridmap,kStop1,ip);
		for (tskidmit=tskidmap.begin();tskidmit!=tskidmap.end();tskidmit++) {

			kStop1 = tskidmit->second;
			kStop1.serializetext(logfile);
			blnHistoric = InpFiles.get_historic() ;	// get if historic is to be run first before alternative analysis

			if (blnPd) { // if it is a period run
				phwsmit  = phsmap.find(kStop1.tripPeriod()); //toString<long>(kStop1.tripId())); 
			} else {
				phwsmit  = phsmap.find(to_string<long> (kStop1.tripId())); //toString<long>(kStop1.tripId())); 
			}
			if (phwsmit!=phsmap.end()) {
				phw1 = phwsmit->second;
				blnTripRun=phw1.get_include(); 
				if (!blnTripRun) { // if the trip or period is included in the selecton 
					phwsmit  = phsmap.find(toString<long>(kStop1.tripId()));
					if (phwsmit !=phsmap.end()) { 
						blnTripRun = (phwsmit->second.get_include());
					}
				}
				if (blnTripRun) { // if this period or trip is included in the analysis

					pid = phw1.get_pdId();
					fprintf( stderr,"\n\tStarting Trip/Period %s and trip Number %s from table %s ; Route %s ;  and Dirn %s Stop Spacing Analysis ... \n", kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str());
					fprintf( txtStream,"\n\tTrip/Period %s and trip Number %s   table %s ; Route %s ;  and Dirn %s Stop Spacing Analysis ... \n", kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str());

					M = gcost->get_maxskip() + 1;
					D = gcost->get_dpdimension();
					strTmPd = kStop1.tripPeriod() ;				
				// strSQLTable - SQL statement to query the land use table data
					strSQLTable =  " SELECT Ptype, Desc, LU , KeyFld, OffCoef, OnCoef, PdOveRide  "
						" FROM  " + strLUCodeName + " "
						" Order by Ptype ";
					
					// get landuse factors based on time period (AM / PM) and zero out Values 
					// AM PD - the attracton for RES - and Producton for COM and Vice Versa for PM  
					lumap = readLUCodesTable(pTransInpDb, strSQLTable,str1,*luci,lumap,txtStream,strTmPd);
					
					if (blnScene) {
						strcTransOutDbName = InpFiles.get_resultDB(); //+ "_" + kStop1.tripPeriod() + "_t" + to_string<long>(kStop1.tripId()) + "_LTS" + to_string<int>(InpFiles.get_lts()) +".sqlite"; // networkDbName.substr(0,networkDbName.find_last_of(".")) ;
					} else {
						strcTransOutDbName = InpFiles.get_resultDB() + "_" + kStop1.tripPeriod() + "_t" + to_string<long>(kStop1.tripId()) + "_LTS" + to_string<int>(InpFiles.get_lts()) +".sqlite"; // networkDbName.substr(0,networkDbName.find_last_of(".")) ;
					}
					replace(strcTransOutDbName.begin(),strcTransOutDbName.end(),' ','_');

					//strcTransOutDbName = strcTransOutDbName.append("_Transit.db");
					//const char* cTransOutDbName = new char[strcTransOutDbName.length() + 1];
					//cTransOutDbName = strcTransOutDbName.c_str();
					rc = createDb(strcTransOutDbName,pTransOutDb,txtStream);
			/*		if (rc) {
						fprintf( stderr,"\tSQLite error opening/creating spatialite Convex Hull output database : %s ; error %s \n", cTransOutDbName, sqlite3_errmsg( networkDbase));
						fprintf( txtStream,"\tSQLite error opening/creating spatialite Convex Hull output database : %s ; Error %s \n", cTransOutDbName, sqlite3_errmsg( networkDbase));
						exit(1);
					}
			*/		
					rc = sqlite3_open_v2( strcTransOutDbName.c_str(), &pTransOutDb, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );
					if( rc != SQLITE_OK){
						fprintf( stderr,"\tSQLite error opening spatialite transit optimization output database : %s ; \n error %s \n", strcTransOutDbName.c_str(), sqlite3_errmsg( pTransOutDb));
						fprintf( txtStream,"\tSQLite error opening spatialite transit optimization output database : %s ;\n error %s \n",strcTransOutDbName.c_str(), sqlite3_errmsg( pTransOutDb));
						exit(1);
					}

					// attach the input and result databases for making the geoemtric query and writing to the output database
					strADB = "netDB";
					blnS = sqliteAttachDB(pTransOutDb, networkDbName,strADB,txtStream);
					if ( blnS )
					{
						fprintf(txtStream,"Table %s is attached as %s !\n SQLite Message : %s ", networkDbName.c_str(),strADB.c_str(), sqlite3_errmsg( pTransOutDb));
					} else {
							fprintf(txtStream,"Table %s could not be attached ! Routine ending \n SQLite Error : %s! ", networkDbName.c_str(), sqlite3_errmsg( pTransOutDb));
							exit(1);
					}

					// if historic scenario needs to run first then read the raw input files/tables else skip 
					//and read binary files from previous run

					if (blnHistoric) { // i0
				
						// perform the historic set run by using the largest stop set in the trips table 
						if (blnPd) {
							strSQLTable = "Select StOrdr, SchlName, RteName, DirName, TimePeriod, TripSTime, TripId," 
										"  TripNumber, SchlTime, AcTimeArr, AcTimeDep, SchlRunTm, Stop_Id, StrMainCrs, " 
										"  Ons, Offs, st_x(Geometry) x, st_y(Geometry) y, EdgeId, PosAlong, arrDelay, "
										" depDelay, NearDist,hist, eliminate, Geometry  from  " + InpFiles.get_tbltrip() + " "
										" where RteName like '" + InpFiles.get_rtename() + "' and "
										" dirName like '" + InpFiles.get_dirname() + "%' and " 
										" timePeriod = '" + kStop1.tripPeriod() + "' and hist <> 0 and "
										" tripNumber = " + to_string<long>(kStop1.tripNumber()) + "  "
										" Order by StOrdr ;";
						} else {
							//strSQLTable = "Select StOrdr, SchlName, RteName, DirName, TimePeriod, TripSTime, TripKey," 
							//		"  TripNumber, SchlTime, AcTimeArr, AcTimeDep, SchlRunTm, Stop_Id, StrMainCrs, " 
							//		"  Ons, Offs, st_x(Geometry) x, st_y(Geometry) y, EdgeId, PosAlong, arrDelay, "
							//		" depDelay, NearDist, Hist, Geometry from  " + InpFiles.get_tbltrip() + " "
							//		" where RteName like '" + InpFiles.get_rtename() + "' and "
							//		" dirName like '" + InpFiles.get_dirname() + "%' and  " 
							//		" tripNumber = '" + to_string<long>(kStop1.tripNumber()) + "' "
							//		" Order by StOrdr ;";
							strSQLTable = "Select StOrdr, SchlName, RteName, DirName, TimePeriod, TripSTime, TripId," 
									"  TripNumber, SchlTime, cRdTm, undCRdTm, cumDist,  Stop_Id, Stopname, " 
									"  Ons, Offs, st_x(Geometry) x, st_y(Geometry) y, EdgeId, PosAlong, arrDelay, "
									" depDelay, SchlRunTm, Hist, dwlDelay, rideDelay, eliminate, Geometry from  " + InpFiles.get_tbltrip() + " "
									" where RteName like '" + InpFiles.get_rtename() + "' and hist >= 0 and "
									" dirName like '" + InpFiles.get_dirname() + "%' and  " 
									" tripNumber = " + to_string<long>(kStop1.tripNumber()) + " "
									" Order by StOrdr ;";
						}
						tsidmap.clear();
						// stop table name kStop1.route() +kStop1.schlName()+ to_string<long>(kStop1.tripId())+kStop1.tripPeriod()
						// choose only the historic set of stops first 
						if (blnPd) {
							tsidmap = stopTableRTD(pTransInpDb, strSQLTable,*pstop,tsidmap,kStop1,gcost,phw1,blnPd,intSRID,outfile,sElim);
						} else {
							tsidmap = stopTableRTD2(pTransInpDb, strSQLTable,*pstop,tsidmap,kStop1,gcost,phw1,blnPd,intSRID,outfile,sElim);
						}
						//						assert(tsidmap.size());
						
						if (tsidmap.size()<=1 )
						{
							fprintf(stderr, "\nSQLite error: %s\nSQL: %s \n", strSQLTable.c_str() ,  sqlite3_errmsg( pTransInpDb) );
							fprintf(txtStream, "\nSQLite error: %s\nSQL: %s \n", strSQLTable.c_str() ,  sqlite3_errmsg( pTransInpDb) );
							exit(1);
						}

						if (outfile.is_open()) {
							outfile.close();
						}
						outfile.clear();

					// save the stop data to a binary file
						strStopFileName= strStopBaseName + ".bin";
						//fileName(tstopiname,".bin",outfilename);
						outfile.open(strStopFileName.c_str(), ios::out | ios::binary|ios::trunc);
						writeBinObjectData(tsidmap,stop0,i, outfile);
						strStopFileName= strStopBaseName ;

						if (outfile.is_open()) {
							outfile.close();
						}
						outfile.clear();
						logfile<<endl<<" Stop Table Input Data"<<endl;
						stop0.serializetexthdr(logfile);
						writeTextObjectData(tsidmap,stop0,ip,logfile,"");


					// three dimensional pattern generator used only for reporting
						if (D==3) {
							PatternGenerator3d(tsidmap,M,strStopBaseName,strStopFileName);
						} else if (D==5) {
					//readtstopMapData(tsidmap,*pstop,ip,tstopfile,xs);
					// five dimensional pattern generator used only for reporting 
							PatternGenerator5d(tsidmap,M,strStopBaseName,strStopFileName);
						}
					/* get the maximum ride time and calculate the ride time from boarding station 
					to the end of the line */

						tsrtmap.clear();
						tsidmit = tsidmap.begin();
						while (tsidmit!=tsidmap.end()) {
								tsrtmap.insert(tsrtpair(tsidmit->second.get_CRdTm(),tsidmit->second));
								tsidmit++;
						}
						rideTime2End(tsrtmap);


					// save the stop data to a binary file after the ridetime to end is calculated
						
						strStopFileName= strStopBaseName + ".bin";
						//fileName(tstopiname,".bin",outfilename);
						//outfile.open(outfilename, ios::out | ios::binary|ios::trunc);
						if (outfile.is_open()) {
							outfile.close();
						}
						outfile.clear();
						outfile.open(strStopFileName.c_str(), ios::out | ios::binary|ios::trunc);
						writeBinObjectData(tsrtmap,stop0,xcost, outfile);

						if (outfile.is_open()) { // i1
							outfile.close();
						} // 1i
						outfile.clear();

					// transfer the objects with the calculated runtime to the end into the main collection (tsidmap)
						tsidmap.clear();
						remapObj2Id(tsrtmap,tsidmap,stop0,ip);

						//Read Parcel data
						tblStopBuf = "s" + kStop1.route()+kStop1.schlName() + kStop1.dir()+to_string<long>(kStop1.tripId()) + kStop1.tripPeriod(); 
						replace(tblStopBuf.begin(),tblStopBuf.end(),' ','_');
						replace(tblStopBuf.begin(),tblStopBuf.end(),'(','_');
						replace(tblStopBuf.begin(),tblStopBuf.end(),')','_');
						replace(tblStopBuf.begin(),tblStopBuf.end(),':','_');
						replace(tblStopBuf.begin(),tblStopBuf.end(),'-','_');
						replace(tblStopBuf.begin(),tblStopBuf.end(),'/','x');
						ReplaceAll2(tblStopBuf,"__","_");
						tblStopBuf += "_Buf";
						//strSQLTable = " SELECT  t1.ID, t1.GEOID10, t1.STATE, t1.COUNTY, t1.PLACE, t1.TRACT, t1.BLOCK, "
						//	" t1.POP, t1.HU, t1.EdgeID, t1.Ptype, t1.PtypeUR, trim(Round(st_x(st_centroid(t1.Geometry)),3)) xC, "
						//	" trim(Round(st_y(st_centroid(t1.Geometry)),3)) yC, " 
						//	" trim(round(st_area(t1.Geometry),3)) Area, t1.PosAlong, t1.NearDist "
						if (blnParcel) {
							strSQLTable = " SELECT  t1.ID, t1.PID, t1.Parcel_id, "
							" t1.POPOn, t1.EdgeID, t1.LU, t1.LUDesc, trim(Round(st_x(st_centroid(t1.Geometry)),3)) xC, "
							" trim(Round(st_y(st_centroid(t1.Geometry)),3)) yC, " 
							" trim(round(st_area(t1.Geometry),3)) Area, t1.PosAlong, t1.NearDist, t1.PopOff "
							" FROM  " + InpFiles.get_tblparcel() + " t1 , " + tblStopBuf + " t2 "
							" where st_within(st_centroid(t1.Geometry),t2.Geometry) "
							" ORDER BY ID; ";
						} else {
							strSQLTable = " SELECT  t1.ID, t1.GEOID10, t1.BLOCK, "
							" t1.POP, t1.EdgeID, t1.Ptype, t1.LU, trim(Round(st_x(st_centroid(t1.Geometry)),3)) xC, "
							" trim(Round(st_y(st_centroid(t1.Geometry)),3)) yC, " 
							" trim(round(st_area(t1.Geometry),3)) Area, t1.PosAlong, t1.NearDist "
							" FROM  " + InpFiles.get_tblparcel() + " t1 , " + tblStopBuf + " t2 "
							" where st_within(st_centroid(t1.Geometry),t2.Geometry) "
							" ORDER BY ID; ";
						} // else {
						//strSQLTable = " SELECT  t1.ID, t1.GEOID10, t1.BLOCK, "
						//	" t1.POP, t1.EdgeID, t1.Ptype, t1.LU, trim(Round(st_x(st_centroid(t1.Geometry)),3)) xC, "
						//	" trim(Round(st_y(st_centroid(t1.Geometry)),3)) yC, " 
						//	" trim(round(st_area(t1.Geometry),3)) Area, t1.PosAlong, t1.NearDist "
						//	" FROM  " + InpFiles.get_tblparcel() + " t1 "
						//	" ORDER BY ID; ";
						//}

				// table reading routine that is statement based query
						mmaparced.clear();
						mmaparced = parcelMapTableState(pTransInpDb, strSQLTable,*parp,mmaparced,txtStream);
						if (mmaparced.size() == 0 ) {
							fprintf(stderr, "\nSQLite error: %s\nSQL: %s \n", strSQLTable.c_str() ,  sqlite3_errmsg( pTransInpDb) );
							fprintf(txtStream, "\nSQLite error: %s\nSQL: %s \n", strSQLTable.c_str() ,  sqlite3_errmsg( pTransInpDb) );
							exit(1);
						}
						// write a binary parcel object data
						if (outparcfile.is_open()) {
							outparcfile.close();
						}
						outparcfile.clear();
						//strParcelBaseName = InpFiles.get_fileparcel();
						strParcelFileName = strParcelBaseName + ".bin";
						//fileName(pariname,".bin",outfilename);
						outparcfile.open(strParcelFileName.c_str(), ios::out | ios::binary|ios::trunc);
						writeBinObjectData(mmaparced,par1,ip,outparcfile);

						if (outparcfile.is_open()) {
							outparcfile.close();
						}
						outparcfile.clear();

						if (blnEuclid) {
							mmapStopVx.clear();
							mmaped.clear();
							mmapVx1.clear();
							mapvert0.clear();
							mapvert1.clear();
							// start building the vertex and edge network using stop and parcel data
							vxveuclid(mmaparced,tsidmap,mmapStopVx,mapvert0,mapvert1,mmapVx1,gc1,blnHistoric);

						} // 1i end euclid historic run
						else 
						//if (!blnEuclid) 
						{ // 1i not an eucledian run

						  //	 openfile(outfile);
							//fext=ext; // file extension
							//strEdgeBaseName = InpFiles.get_fileedge() ;
							//fileName(edgeinfname,".edg",outfilename);
							strEdgeFileName = strEdgeBaseName + ".edg";
							outedgefile.clear();
							outedgefile.open(strEdgeFileName.c_str(), ios::out | ios::trunc);

							//fileName(edgeinfname,".bdy",outfilename);
							strEdgeFileName = strEdgeBaseName + ".bdy";
							outedgebdyfile.clear();
							outedgebdyfile.open(strEdgeFileName.c_str(), ios::out | ios::trunc);

							j=0, nv=0;
							tblStopBufx2 = tblStopBuf + "x2";
//						if (blnPd) {

							strSQLTable = " SELECT  t2.EDGEID, t2.EDGESID, t2.VX1, t2.VX2, st_length(t2.Geometry) LENGTH, "
							" t2.NAME, t2.COST,  t2.spd_lim, t2.COMMENT, t2.CLASS, t2.ROW_WIDTH, t2.FOC_WIDTH, "
							" t2.ONEWAY, t2.PRIVATE, t2.FCODE, t2.FACTYPE,  t2.PATHWIDTH, t2.PARKWIDTH, "
							" t2.NUMLANE, t2.ILLPARKING, t2.ADT, t2.RTLANE, t2.SIGNAL, t2.MEDWIDTH, " 
							" t2.SHARESPACE, t2.CL, t2.RTSTRESS, t2.OVERIDE "
							" FROM " + InpFiles.get_tbledge() + " t2 " + "  , " + tblStopBufx2 + " t1 "
							" where (t2.RTSTRESS <= " + to_string<int>(InpFiles.get_lts()) + " or t2.OVERIDE > 0 ) "
							" and st_within(st_centroid(t2.Geometry),t1.Geometry)  "
							" ORDER BY t2.edgeid ; ";
						//} else {
						//	strSQLTable = " SELECT  t2.EDGEID, t2.EDGESID, t2.VX1, t2.VX2, st_length(t2.Geometry) LENGTH, "
						//	" t2.NAME, t2.COST,  t2.spd_lim, t2.COMMENT, t2.CLASS, t2.ROW_WIDTH, t2.FOC_WIDTH, "
						//	" t2.ONEWAY, t2.PRIVATE, t2.FCODE, t2.FACTYPE,  t2.PATHWIDTH, t2.PARKWIDTH, "
						//	" t2.NUMLANE, t2.ILLPARKING, t2.ADT, t2.RTLANE, t2.SIGNAL, t2.MEDWIDTH, " 
						//	" t2.SHARESPACE, t2.CL, t2.RTSTRESS, t2.OVERIDE "
						//	" FROM " + InpFiles.get_tbledge() + " t2 "  
						//	" where (t2.RTSTRESS <= " + to_string<int>(InpFiles.get_lts()) + " or t2.OVERIDE > 0 ) "
						//	" ORDER BY t2.edgeid ; ";
						//}
							maped.clear();
							maped = readedgeTableMapData (pTransInpDb, strSQLTable,ip,*evp,maped,txtStream,*gcost);
							
							// write the edge object data
							//fileName(edgeinfname,".bin",outedgefilename);
							strEdgeFileName = strEdgeBaseName + ".bin";

							outedgeobjectfile.clear();

							outedgeobjectfile.open(strEdgeFileName.c_str(), ios::out | ios::binary|ios::trunc);
							writeBinObjectData(maped,ev,ip,outedgeobjectfile);
							if (outedgeobjectfile.is_open()) {
								outedgeobjectfile.close();
							}
							outedgeobjectfile.clear();

					  // open vertex data file  -  
							//strVertexBaseName = InpFiles.get_filevertex();
							//fileName(vertinfname,".vex",outfilename);
							strVertexFileName = strVertexBaseName + ".vex";
							if (outvertfile.is_open()) {
								outvertfile.close();
							}
							outvertfile.clear();
							outvertfile.open(strVertexFileName.c_str(), ios::out | ios::trunc);
							j=0, nv=0;
							while (!outvertfile) 
							{ //w1
								cout << "Vertex Output file could not be created at "<<strVertexFileName;
								cout << "Please reenter file name again : "<<endl; 
    							outvertfile.clear();
    							outvertfile.open(strVertexFileName.c_str(), ios::out | ios::binary|ios::trunc);
								if (++j>3) 
								{ //i1 
									cout << "More than three trials!"<<endl;
									cout << "Press enter to exit!"<<endl;
									cin>>strVertexFileName;
									exit (0);
								} //1i
							} //1w

							// read vertex data

							mapvert3.clear();

//						if (blnPd) {

							strSQLTable = "  SELECT t1.VERTID, t1.EDGEID, t1.IDP, t1.COST, t1.SIGNAL,st_x(t1.Geometry) X, "
							" st_y(t1.Geometry) Y, t1.Geometry FROM  " + InpFiles.get_tblvertex() + " t1 " 
							" where t1.Vertid in (select t3.vx1 from " + InpFiles.get_tbledge() + " t3 , " + tblStopBufx2 + " "
							" t4 where (t3.RTSTRESS <= " + to_string<int>(InpFiles.get_lts()) + " or t3.OVERIDE > 0 ) "
							" and st_within(st_centroid(t3.Geometry),t4.Geometry) ) Or t1.VertId in (select t5.vx2 from " 
							" " + InpFiles.get_tbledge() + " t5 , " + tblStopBufx2 + " t6 where (t5.RTSTRESS <= " +  " "
							" " + to_string<int>(InpFiles.get_lts()) + " or t5.OVERIDE > 0 ) " 
							" and st_within(st_centroid(t5.Geometry),t6.Geometry) ) ORDER BY t1.vertId ";
						//} else {
						//	cout << "Source Vertex Query : "<<strSQLTable<<endl;
						//	strSQLTable = "  SELECT t1.VERTID, t1.EDGEID, t1.IDP, t1.COST, t1.SIGNAL,st_x(t1.Geometry) X, "
						//	" st_y(t1.Geometry) Y, t1.Geometry FROM  " + InpFiles.get_tblvertex() + " t1 " 
						//	" where t1.Vertid in (select t3.vx1 from " + InpFiles.get_tbledge() + " t3 , " + tblStopBufx2 + " "
						//	" t4 where (t3.RTSTRESS <= " + to_string<int>(InpFiles.get_lts()) + " or t3.OVERIDE > 0 ) "
						//	" and st_within(st_centroid(t3.Geometry),t4.Geometry) ) Or t1.VertId in (select t5.vx2 from " 
						//	" " + InpFiles.get_tbledge() + " t5 where (t5.RTSTRESS <= " 
						//	" " + to_string<int>(InpFiles.get_lts()) + " or t5.OVERIDE > 0 ) " 
						//	"  ) ORDER BY t1.vertId ";

						//}
							outfile << "Source Vertex Query : "<<strSQLTable<<endl;

							mapvert3 = readvertexTableMapData (pTransInpDb, strSQLTable,ip,vx,mapvert3,txtStream);
							
						} // 1i not eucledean run therefore read edge and vertex data from table/files 

					// use only the vertex that is an end point of an edge and write to a binary file 
						// get a multimap of the end point vertices 
						mmapV1V2.clear();
						mmapV1V2 = funcV1V2(maped,mmapV1V2,*evp,i,j);
						mapvert.clear();
						vxmap_Iter = mapvert3.begin();
						for (vxmap_Iter=mapvert3.begin();vxmap_Iter!=mapvert3.end();++vxmap_Iter)
						{ // f1
							ip = vxmap_Iter->first;
							mumavedit = mmapV1V2.find(ip);
							if  (mumavedit != mmapV1V2.end()) { //i1
								vx = vxmap_Iter->second;
								mapvert.insert(vx_pair(ip,vx));
							}  // 1i
						} // 1f
						
					// write the vertex object data as a binary file 

						strVertexFileName = strVertexBaseName + ".bin";
						//fileName(vertinfname,".bin",outfilename);
						if (outfile.is_open()) {
							outfile.close();
						}
						outfile.clear();
						outfile.open(strVertexFileName.c_str(), ios::out | ios::binary|ios::trunc);
						writeBinObjectData(mapvert,vx,ip,outfile);
						if (outfile.is_open()) {
							outfile.close();
						}
						outfile.clear();

					// set the stop -> edge -> vertex relationship
					 // first make a multimap edgeid and edge object from the map of edge Object id & edge object
						if (!blnEuclid) {  // i1
							mapEidOid = funcEdgeId2ObjId(maped,mapEidOid,*evp,i,j);
							mmapV1V2EgId = funcV1V2EdgeId(maped,mmapV1V2EgId,*evp,i,j);
							mmapV1V2EOId = funcV1V2EOId(maped,mmapV1V2EOId,*evp,i,j);
							mmapeidvid = funcEdgeIdV1V2(maped,mmapeidvid,*evp,i,j);
							mmapEOIdEgId = funcEdgeObjId2Id(maped,mmapEOIdEgId,*evp,i,j);
						// assign the stop starting points to the vertices and edges
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

							maped =  funCoStopVxEdge(maped,mapEidOid,ev,tsidmap,*pstop,mapvert,*pVx1,mmapVxStop,
							mmapStopVx,mmapStopEdge,mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,mmapVx0R,mmapVx1R,
							mapvertor0,mapvertor1,outedgefile);
						// mmapStopEdge,mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,mmapVx0R,mmapVx1R,mapStopVx,mmapVxStop
							onoff=0;
							strVertexFileName = strVertexBaseName + "_" + to_string<long>(onoff) +"Vx.bin";
							
							//fileName(vertinfname,"Vx.bin",outfilename,onoff);
							if (outfile.is_open()) {  // i2
								outfile.close();
							} // 2i
							outfile.clear();

							outfile.open(strVertexFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindualvar(mmapVx0,xcost,i,outfile);
							if (outfile.is_open()) { // i2
								outfile.close();
							} // 2i
							outfile.clear(); 

						// this is the reverse of mmapVx0
							//fileName(vertinfname,"VxR.bin",outfilename,onoff);
							strVertexFileName = strVertexBaseName + "_" + to_string<long>(onoff) +"VxR.bin";
							if (outfile.is_open()) { // i2
								outfile.close();
							} // 2i
							outfile.clear();

							outfile.open(strVertexFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindualvar(mmapVx0R,i,xcost,outfile);
							if (outfile.is_open()) {  // i2
								outfile.close();
							}  // 2i
							outfile.clear(); 

							onoff=1;
							//fileName(vertinfname,"Vx.bin",outfilename,onoff);
							strVertexFileName = strVertexBaseName + "_" + to_string<long>(onoff) +"Vx.bin";
							if (outfile.is_open()) {  // i2
								outfile.close();
							} // 2i
							outfile.clear();

							outfile.open(strVertexFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindualvar(mmapVx1,xcost,i,outfile);
							if (outfile.is_open()) {  // i2
								outfile.close();
							}  // 2i
							outfile.clear(); 

						// this is the reverse of mmapVx1
							//fileName(vertinfname,"VxR.bin",outfilename,onoff);
							strVertexFileName = strVertexBaseName + "_" + to_string<long>(onoff) +"VxR.bin";
							if (outfile.is_open()) {  // i2
								outfile.close();
							}  // 2i
							outfile.clear();

							outfile.open(strVertexFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindualvar(mmapVx1R,i,xcost,outfile);
							if (outfile.is_open()) {  // i2
								outfile.close();
							}  // 2i
							outfile.clear(); 


						//write to origin vertex data output file - empty file if it already exists before writing new data
				 
							//fileName(tstopiname,"_Orig.bin",outfilename);
							strStopFileName = strStopBaseName + "_Orig.bin";
							ofstream originout(strStopFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindblong(mmapSver,originout);
							if (originout.is_open()) { // i2
								originout.close();
							}  // 2i
							originout.clear(); 

						//open vertex-stop map data output output file - empty file if it already exists before writing
							//fileName(tstopiname,"_VxStop.bin",outfilename);
							strStopFileName = strStopBaseName + "_VxStop.bin";
							originout.open(strStopFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindualvar(mmapVxStop,v1,v2,originout);
							if (originout.is_open()) {  // i2
								originout.close();
							}  // 2i
							originout.clear(); 

				 //open stop-vertex map data output output file - empty file if it already exists before writing
							//fileName(tstopiname,"_StopVx.bin",outfilename);
							strStopFileName = strStopBaseName + "_StopVx.bin";
							originout.open(strStopFileName.c_str(),  ios::out | ios::binary|ios::trunc);
							writeBindualvar(mmapStopVx,v1,v2,originout);
							if (originout.is_open()) {  // i2
								originout.close();
							}  // 2i
							originout.clear(); 
						} // i1 not euclidean run
					} // i0 set up of historic data input/ read if historic run is required  
			// start analysis here 
				if (!blnEuclid) {  // i0
					//fileName(tstopiname,"_Orig.bin",outfilename);
					strStopFileName = strStopBaseName + "_Orig.bin";
					mmapSver2.clear();
					infile.open(strStopFileName.c_str(),  ios::in | ios::binary);
					mmapSver2= readBindblong(mmapSver2,infile);
					if (infile.is_open()) { // i1
						infile.close();
					 } // 1i
					infile.clear();

			 // output headers for the edge and vertex files
					pVx1->show_verthdr(outvertfile);
					evp->show_edgehdr(outedgefile);
				} // 0i blnEuclid 

			// historic analysis
				int alt=0;

				if (blnHistoric) { // i0 begin historic analysis
			// assign the period default period is 3
					phwmit = phmap.find(pid);
					if (phwmit!=phmap.end()) { // i1
						*phway = phwmit->second;
					} // 1i
					mmaparStop.clear();
					tsidmap = stopspacing_analysis (mmapV1V2EOId,mmapV1V2EgId, mmaparced, mapvert, maped,mmaparStop,
					mmapSver2, mmapVxStop,mmapStopVx, mmapVx0, mmapVx1, mmapVx0R, mmapVx1R, blnHistoric,blnEuclid,
					strEdgeBaseName, strVertexBaseName, strParcelBaseName, strStopBaseName, gc1, *phway, *luci, tsidmap, 
					alt, logfile, InpFiles,kStop1, pTransOutDb, strADB,xi);

			//open binary output stop set file - empty file if it already exists before writing new data
					 //fileName(tstopiname,"alt.bin",outfilename,alt);
					strStopFileName = strStopBaseName + "_" + to_string<int>(alt) +"alt.bin";
					 stopout.open(strStopFileName.c_str(), ios::out | ios::binary|ios::trunc);
					 writeBinObjectData(tsidmap,stop0,ip,stopout);
					 if (stopout.is_open()) { // i1
						stopout.close();
					 } // 1i
					 stopout.clear();


					 cout <<endl<<"Historic Analysis is completed."<<endl;

				// end historic analysis here
				} // 0i
				// start the alternative analysis
				// maximum number of stops to be removed l = 4, q - serial id, k - pattern repeatition
					 i=0;m=0;n=0;o=0;p=0;q=0; onoff=0;
				// open and write the DP pattern into a text file for reference
				 q=-1;
				strStopFileName = strStopBaseName + "_dpattern.csv";
				 //fileName(tstopiname,"_dpattern.csv",outfilename,q);
				 //stopout.open(strStopFileName.c_str(), ios::out |ios::trunc);
				 // read the stop data into the alternative map data set
				tsidmap.clear();
				q=0;
				//fileName(tstopiname,"alt.bin",outfilename,q);
				if (tstopfile.is_open()) {  // i1
					tstopfile.close();
				}  // 1i
				tstopfile.clear();
				strStopFileName = strStopBaseName + "_" + to_string<long>(q) +"alt.bin";
				tstopfile.clear();
				tstopfile.open(strStopFileName.c_str(), ios::in | ios::binary);
				readBinStopObjectData(tsidmap,stop0,tstopfile);
				if (tstopfile.is_open()) {  // i1
					tstopfile.close();
				}  // 1i
				tstopfile.clear();

				if (!blnHistoric) { // i1
					if (!blnEuclid) { // i2 not eucledean run
						if (infile.is_open()) { // i3
							infile.close();
						}  // 3i
						infile.clear();
						strEdgeFileName = strEdgeBaseName + ".edg";
						outedgefile.clear();
						outedgefile.open(strEdgeFileName.c_str());
						strEdgeBinName = strEdgeBaseName + ".bin";
						infile.open(strEdgeBinName.c_str(), ios::in | ios::binary);
						maped.clear();
						readedgeMapData(maped,*evp,ip,infile,gc1,xs);
						// close edge data if open
						if (infile.is_open()) { // i3
							infile.close();
						} // 3i
						infile.clear();

						// write the edge object data
						//fileName(edgeinfname,".bin",outedgefilename);
						strEdgeFileName = strEdgeBaseName + ".bin";

						outedgeobjectfile.clear();

						outedgeobjectfile.open(strEdgeFileName.c_str(), ios::out | ios::binary|ios::trunc);
						writeBinObjectData(maped,ev,ip,outedgeobjectfile);
						if (outedgeobjectfile.is_open()) { // i3
							outedgeobjectfile.close();
						} // 3i
						outedgeobjectfile.clear();

						//fileName(vertinfname,".bin",outfilename);
						strVertexFileName = strVertexBaseName + ".bin";
						infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
						readBinVertexObjectData(mapvert,ip,infile);
						if (infile.is_open()) { // i3
							infile.close();
						}  // 3i
						infile.clear();


						mapEidOid = funcEdgeId2ObjId(maped,mapEidOid,*evp,i,j);
						mmapV1V2 = funcV1V2(maped,mmapV1V2,*evp,i,j);
						mmapV1V2EgId = funcV1V2EdgeId(maped,mmapV1V2EgId,*evp,i,j);
						mmapV1V2EOId = funcV1V2EOId(maped,mmapV1V2EOId,*evp,i,j);
						mmapeidvid = funcEdgeIdV1V2(maped,mmapeidvid,*evp,i,j);
						mmapEOIdEgId = funcEdgeObjId2Id(maped,mmapEOIdEgId,*evp,i,j);
						// assign the stop starting points to the vertices and edges

						if (qNet) { // i3
							// Using the stop set, fix the vertex and edge costs that are the entry points into the network.  
							maped =  fixVxEdgeCoStop(maped,mapEidOid,ev,tsidmap,*pstop,mapvert,*pVx1,mmapVxStop,
							mmapStopVx,mmapStopEdge,mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,mmapVx0R,mmapVx1R,
							mapvertor0,mapvertor1,outedgefile);
						} else { // i3
							maped =  funCoStopVxEdge(maped,mapEidOid,ev,tsidmap,*pstop,mapvert,*pVx1,mmapVxStop,
							mmapStopVx,mmapStopEdge,mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,mmapVx0R,mmapVx1R,
							mapvertor0,mapvertor1,outedgefile);
						} // 3i
					} // 2i  not eucledean run
				} // 1i Historic is not run 
				q=0;
				long dpi = 0 ;
				int z=M+1;


				n = (int) tsidmap.size();
				tsidpmap.clear();
			//	scenarioname
				//sceneifile.open(strsceneFileName.c_str());
				// build the DP table name 
				pdNm = kStop1.tripPeriod();
				tripNum = to_string<long>(kStop1.tripId());

				if ( blnScene ) 
				{ // i1 scenario file is provided
					// set q to a high value so that it does not overwrite the previous run scenarios and run alternatives
					q=100000;

				//Read scenario list data 
					//	scenarioname
					if (blnSceneDB) {  // this is a database file list
						// open the dp results table and get the stops selected
						strWR = ("W" + to_string<int> ((int) (gcost->get_walkcost()*10)));
						strWR.append("R" + to_string<int> ((int) (gcost->get_ridecost()*10)));
						// rec1 to be used to find trip and period id 
						rec1 = InpFiles.get_tblstop()   ;
						//strDPTblName = inParams.get_tblstop() + "D" + to_string<short> (D) + "t" + to_string<long> (kStop1.tripId())+ kStop1.tripPeriod() + "W" + to_string<int> ((int) (gcost->get_walkcost()*10))+ "R" + to_string<int> ((int) (gcost->get_ridecost()*10)); 
						if (blnPdNm) {
							tblDP = rec1 + "D" + to_string<short> (D) + "t" + tripNum + pdNm + strWR + "%DPTrace"; 
						} else {

							tblDP = rec1 + "%" + strWR + "%DPTrace"; 
						}
						replace(tblDP.begin(),tblDP.end(),' ','_');
						ReplaceAll2(tblDP,"__","_");
						strDPTblName = tblDP;
						// loop over all the tables in the DPResult database and with the filter string & run individual trip/Period scenario 
						str1= " select tbl_name ,  name  , rootpage  from sqlite_master " 
						" where tbl_name like '" + strDPTblName +"' and type like 'table' "
						" ORDER BY tbl_name; " ;
						if ( sqlite3_prepare_v2( pSceneDb, str1.c_str(), -1, &stmtDPTables, NULL ) != SQLITE_OK )
						{
						  // some error occurred
						  fprintf(stdout, "\nSQLite error: %s\n\nSQL: %s \n", str1.c_str() ,  sqlite3_errmsg(pSceneDb) );
						  fprintf(stderr, "\nSQLite error: %s\n\nSQL: %s \n", str1.c_str() ,  sqlite3_errmsg(pSceneDb) );
						  fprintf(txtStream,"\nSQLite error: %s\n\nSQL: %s \n", str1.c_str() ,  sqlite3_errmsg(pSceneDb) );  //"SQLite error:  "<<q1<< " !"<<endl<<to_string<const char *>(sqlite3_errmsg( netdb)) <<endl;
						}
						j = 0;
						// find the maximum number of stop sets among the dp result tables to use as the basis for the optimization
						while ( sqlite3_step( stmtDPTables ) == SQLITE_ROW )
						{
							// get the table name and query the DP table 
							col=0;
							if (sqlite3_column_bytes(stmtDPTables, col) !=0) {
								tblDP  = to_string<const unsigned char * > (sqlite3_column_text(stmtDPTables, col));
								// find the DP K indices of the stOrdr key and store it into a map 
								//str1 =	"SELECT distinct Stopid , StOrdr, K   "
								//	"FROM " + tblDP + " "   
								//	"ORDER BY K, StOrdr ";

								// find the DP K indices of the stOrdr key and store it into a map 
								str1 =	"SELECT distinct Stopid , StOrdr, K   "
									"FROM " + tblDP + " where ons > 0 or offs> 0 "   // only stops with demand
								//	"FROM " + tblDP + " where ons >= 0 or offs>= 0 "   // include all stops in the dp result 
									"ORDER BY K, StOrdr ";
								// make the multi-map for the index k and the StOrdr 
								key = 3;
								mapkLtStp.clear();
								mapkLtStp = dbkStOrdr(pSceneDb, str1,*pltSt,mapkLtStp,txtStream,key);
								j = mapkLtStp.size();
								tsidpmap.clear();
								if (j>0) {
									for (mapkLtIt = mapkLtStp.begin();mapkLtIt != mapkLtStp.end(); ++mapkLtIt)
									{
										i = mapkLtIt->first;
										tsidmit = tsidmap.find(i);
										if (tsidmit!=tsidmap.end()) { // i2
											tsidpmap.insert(tsidpair(i,tsidmit->second));
										} else { // i2
											outfile<<(q-100000)<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<i;
										} // 2i
									}

								// do the analysis using the stop set
									if (!blnEuclid) { // i2
										//fileName(tstopiname,"_Orig.bin",outfilename,-1);
										strStopFileName = strStopBaseName + "_Orig.bin";
										mmapSver2.clear();
										ifstream originin(strStopFileName.c_str(),  ios::in | ios::binary);
										readBindblong(mmapSver2,originin);
										if (originin.is_open()) {  // i3
											originin.close();
										} // 3i
										originin.clear();
										// get the vertex assignment data
										onoff = 0;
										mmapVx0.clear();	
										//fileName(vertinfname,"Vx.bin",outfilename,onoff);
										strVertexFileName = strVertexBaseName +"_" + to_string<long>(onoff) + "Vx.bin";
										originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
										readBindblong(mmapVx0,originin);
										if (originin.is_open()) { // i3
											originin.close();
										} // 3i
										originin.clear();  

										onoff = 1;
										mmapVx1.clear();	
										//fileName(vertinfname,"Vx.bin",outfilename,onoff);
										strVertexFileName = strVertexBaseName +"_" + to_string<long>(onoff) + "Vx.bin";
										originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
										readBindblong(mmapVx1,originin);
										if (originin.is_open()) { // i3
											originin.close();
										}  // 3i
										originin.clear();  

										// transfer the vertex network to the two on and off multi-maps
										// read the vertex object data as a binary file
										infile.clear();
										strVertexFileName = strVertexBaseName + ".bin";
										//fileName(vertinfname,".bin",outvertfilename);
										infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
										mapvert0.clear();
										readBinVertexObjectData(mapvert0,vx,infile);
										// close and open vertex data for the second time
										if (infile.is_open()) {  // i3
											infile.close();
										}  // 3i
										infile.clear();
										infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
										mapvert1.clear();
										readBinVertexObjectData(mapvert1,vx,infile);
										if (infile.is_open()) {  // i3
											infile.close();
										} // 3i
										infile.clear();
										// mmapVx0,1 - the vertex origin collection
										// mmapVxStop - the vertex - stop multi-map 
										// mmapStopVx - the stop - vertex multi-map 
										// mmapSver2 - the sorted walktime - vertex id multi-map 
										// mmapVerlng - vertex id with cost multi-map 
									} // 2i
									blnHistoric = false;
									if (tsidpmap.size()>0) 
									{	// i2
										xi =1;
										tsidpmap = stopspacing_analysis (mmapV1V2EOId, mmapV1V2EgId, mmaparced, mapvert, maped,mmaparStop,
										mmapSver2, mmapVxStop,mmapStopVx, mmapVx0, mmapVx1,mmapVx0R, mmapVx1R, blnHistoric, blnEuclid,
										strEdgeBaseName, strVertexBaseName, strParcelBaseName, strStopBaseName,	gc1, *phway, *luci,tsidpmap, q, logfile,
										InpFiles, kStop1,pTransOutDb, strADB,xi);
									}  // 2i
									tsidpmap.clear();
									cout <<endl<<"Run no. "<<(q-100000)<< " with candidate stops "<< endl<<tblDP <<" is completed."<<endl<<endl;
									q++;
								} // xi scene table has records
							} // table name is present
						} // loop over all rows
					
					} else 
					{  // read from a file with a list of stops in numeric order
						sceneifile.open(strsceneFileName.c_str());
						while(!sceneifile.eof())
						{ //w1
							cout <<endl<<"Starting Alternative Analysis."<<endl<<endl;

							getline(sceneifile,rec1);
							str1 = rec1;
							q++; //pattern serial id
							while (str1.length() >0)
							{ // w2
								// insert stops from the list into the alsternative set (tsidpmap) until all stops are included 
								j = str1.find_first_of("\t");
								if (j>0) { // i2
									string str2 = str1.substr(0,j);
									if (is_number(str2)) {
										i = fromString<long> (str2);
									}
									str1 = str1.substr(j+1);
								} else { // i2
									if (is_number(str1)) {
										i = fromString<long> (str1);
									}
									str1.clear();
								} // 2i
								tsidmit = tsidmap.find(i);
								if (tsidmit!=tsidmap.end()) { // i2
									tsidpmap.insert(tsidpair(i,tsidmit->second));
								} else { // i2
									outfile<<(q-100000)<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<i;
								} // 2i
							} // 2w
						// do the analysis using the stop set
							if (!blnEuclid) { // i2
								//fileName(tstopiname,"_Orig.bin",outfilename,-1);
								strStopFileName = strStopBaseName + "_Orig.bin";
								mmapSver2.clear();
								ifstream originin(strStopFileName.c_str(),  ios::in | ios::binary);
								readBindblong(mmapSver2,originin);
								if (originin.is_open()) {  // i3
									originin.close();
								} // 3i
								originin.clear();
								// get the vertex assignment data
								onoff = 0;
								mmapVx0.clear();	
								//fileName(vertinfname,"Vx.bin",outfilename,onoff);
								strVertexFileName = strVertexBaseName +"_" + to_string<long>(onoff) + "Vx.bin";
								originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
								readBindblong(mmapVx0,originin);
								if (originin.is_open()) { // i3
									originin.close();
								} // 3i
								originin.clear();  

								onoff = 1;
								mmapVx1.clear();	
								//fileName(vertinfname,"Vx.bin",outfilename,onoff);
								strVertexFileName = strVertexBaseName +"_" + to_string<long>(onoff) + "Vx.bin";
								originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
								readBindblong(mmapVx1,originin);
								if (originin.is_open()) { // i3
									originin.close();
								}  // 3i
								originin.clear();  

								// transfer the vertex network to the two on and off multi-maps
								// read the vertex object data as a binary file
								infile.clear();
								strVertexFileName = strVertexBaseName + ".bin";
								//fileName(vertinfname,".bin",outvertfilename);
								infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
								mapvert0.clear();
								readBinVertexObjectData(mapvert0,vx,infile);
								// close and open vertex data for the second time
								if (infile.is_open()) {  // i3
									infile.close();
								}  // 3i
								infile.clear();
								infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
								mapvert1.clear();
								readBinVertexObjectData(mapvert1,vx,infile);
								if (infile.is_open()) {  // i3
									infile.close();
								} // 3i
								infile.clear();
								// mmapVx0,1 - the vertex origin collection
								// mmapVxStop - the vertex - stop multi-map 
								// mmapStopVx - the stop - vertex multi-map 
								// mmapSver2 - the sorted walktime - vertex id multi-map 
								// mmapVerlng - vertex id with cost multi-map 
							} // 2i
							blnHistoric = false;
							if (tsidpmap.size()>0) 
							{	// i2
								xi = 1;
								tsidpmap = stopspacing_analysis (mmapV1V2EOId, mmapV1V2EgId, mmaparced, mapvert, maped,mmaparStop,
								mmapSver2, mmapVxStop,mmapStopVx, mmapVx0, mmapVx1,mmapVx0R, mmapVx1R, blnHistoric, blnEuclid,
								strEdgeBaseName, strVertexBaseName, strParcelBaseName, strStopBaseName,	gc1, *phway, *luci,tsidpmap, q, logfile,
								InpFiles, kStop1,pTransOutDb, strADB,xi);
							}  // 2i
							tsidpmap.clear();
							cout <<endl<<"Run no. "<<(q-100000)<< " with candidate stops "<< endl<<rec1 <<" is completed."<<endl<<endl;

						} // 1w scenefile has records

						if (sceneifile.is_open())
						{ // i2 close scene file 
							sceneifile.close();
						} // 2i
						sceneifile.clear();
						}
					sqlite3_close(pTransOutDb);

					str1 = datetimeStamp(logfile);
					cout<<str1<<endl;


					cout << " End of Bustop Spacing Analysis..."<<endl<<" ..."<<endl;
					//cout<<" Press any key and the Enter Key to exit !"<<endl;
					logfile << " End of Bustop Spacing Analysis..."<<endl<<" ..."<<endl;
					//cin>>ch;
					//exit(0);

				} else 
				{
					 cout <<endl<<"Starting Dynamic Programming Patterned Alternative Scenario Analysis."<<endl;

					// use the trip key collection to run all trips as alternatives
						q=0; //pattern serial id


						// perform the historic set run by using the largest stop set in the trips table 
						if (blnPd) {
							strSQLTable = "Select StOrdr, SchlName, RteName, DirName, TimePeriod, TripSTime, TripId," 
										"  TripNumber, SchlTime, AcTimeArr, AcTimeDep, SchlRunTm, Stop_Id, StrMainCrs, " 
										"  Ons, Offs, st_x(Geometry) x, st_y(Geometry) y, EdgeId, PosAlong, arrDelay, "
										" depDelay, NearDist,Hist,Geometry from  " + InpFiles.get_tblstop() + " "
										" where RteName like '" + InpFiles.get_rtename() + "' and "
										" dirName like '" + InpFiles.get_dirname() + "%' and " 
										" timePeriod = '" + kStop1.tripPeriod() + "' and "
										" tripNumber = " + to_string<long>(kStop1.tripNumber()) + " "
										" Order by StOrdr ;";
						} else {
							//strSQLTable = "Select StOrdr, SchlName, RteName, DirName, TimePeriod, TripSTime, TripKey," 
							//		"  TripNumber, SchlTime, AcTimeArr, AcTimeDep, SchlRunTm, Stop_Id, StrMainCrs, " 
							//		"  Ons, Offs, st_x(Geometry) x, st_y(Geometry) y, EdgeId, PosAlong, arrDelay, "
							//		" depDelay, NearDist,Hist, Geometry from  " + InpFiles.get_tblstop() + " "
							//		" where RteName like '" + InpFiles.get_rtename() + "' and "
							//		" dirName like '" + InpFiles.get_dirname() + "%' and " 
							//		" tripNumber = '" + to_string<long>(kStop1.tripNumber()) + "' "
							//		" Order by StOrdr ;";
							strSQLTable = "Select StOrdr, SchlName, RteName, DirName, TimePeriod, TripSTime, TripId," 
									"  TripNumber, SchlTime, cRdTm, undCRdTm, cumDist,  Stop_Id, Stopname, " 
									"  Ons, Offs, st_x(Geometry) x, st_y(Geometry) y, EdgeId, PosAlong, arrDelay, "
									" depDelay, SchlRunTm, Hist, dwlDelay, rideDelay, eliminate, Geometry from  " + InpFiles.get_tbltrip() + " "
									" where RteName like '" + InpFiles.get_rtename() + "' and "
									" dirName like '" + InpFiles.get_dirname() + "%' and  " 
									" tripNumber = " + to_string<long>(kStop1.tripNumber()) + " "
									" Order by StOrdr ;";

						}
						tsidmap.clear();
						// stop table name kStop1.route() +kStop1.schlName()+ to_string<long>(kStop1.tripId())+kStop1.tripPeriod()
						// choose only the historic set of stops first 
						
						if (blnPd) {
							tsidmap = stopTableRTD(pTransInpDb, strSQLTable,*pstop,tsidmap,kStop1,gcost,phw1,blnPd,intSRID,outfile);
						} else {
							tsidmap = stopTableRTD2(pTransInpDb, strSQLTable,*pstop,tsidmap,kStop1,gcost,phw1,blnPd,intSRID,outfile);
						}

						//						assert(tsidmap.size());
						if (tsidmap.size()<=1 )
						{
							fprintf(stderr, "\nSQLite error: %s\nSQL: %s \n", strSQLTable.c_str() ,  sqlite3_errmsg( pTransInpDb) );
							fprintf(txtStream, "\nSQLite error: %s\nSQL: %s \n", strSQLTable.c_str() ,  sqlite3_errmsg( pTransInpDb) );
							exit(1);
						}

						if (outfile.is_open()) {
							outfile.close();
						}
						outfile.clear();


						//tsidmap = stopTableRTD(pTransInpDb, strSQLTable,*pstop,tsidmap,txtStream);
						/* get the maximum ride time and calculate the ride time from boarding station 
						to the end of the line */
						sdptsidpmmap.clear();
						// re-read the stop table without the historic qualifier 						
						tsidmit = tsidmap.begin();
						tsrtmap.clear();
						while (tsidmit!=tsidmap.end()) {
								tsrtmap.insert(tsrtpair(tsidmit->second.get_CRdTm(),tsidmit->second));
								tsidmit++;
						}
						rideTime2End(tsrtmap);


					// save the stop data to a binary file after the ridetime to end is calculated
						//fileName(tstopiname,".bin",outfilename);
						strStopFileName = strStopBaseName + ".bin";
						outfile.open(strStopFileName.c_str(), ios::out | ios::binary|ios::trunc);
						writeBinObjectData(tsrtmap,stop0,xcost, outfile);

						if (outfile.is_open()) { // i1
							outfile.close();
						} // 1i
						outfile.clear();

					// transfer the objects with the calculated runtime to the end into the main collection (tsidmap)
						tsidmap.clear();
						remapObj2Id(tsrtmap,tsidmap,stop0,ip);

					// three dimensional pattern generator used only for reporting
						if (D==3) {
							PatternGenerator3d(tsidmap,M,strStopBaseName,strStopFileName);
						} else if (D==5) {
					//readtstopMapData(tsidmap,*pstop,ip,tstopfile,xs);
					// five dimensional pattern generator used only for reporting 
							PatternGenerator5d(tsidmap,M,strStopBaseName,strStopFileName);
						}
								 

						 //     gcost->set_dpdimension(5);
						 if ( D==3)
						 {
							str1 =	"( Id Integer PRIMARY KEY AUTOINCREMENT, TripId integer,Walk Integer, Ride Integer, " 
								" Q integer,I integer, J integer, K integer, L integer, M integer,StOrdr Integer, "
								" StopId Text,StopName Text, StopIdp Integer,StopIds Integer,"
								" EdgeId Integer,palong Double, lbl Integer, Hist Integer ,Inbound Integer, "
								" External Integer ,Included Integer ,Eliminate Integer , CumDist Double, "
								" CRdTm Double,CrdTmC Double,HistOns Double,HistOffs Double, "
								"Ons Double,Offs Double,DepVol Double,probStoph Double,probStop Double,"
								"depDelay Double, arrDelay Double,dwlDelay Double,rideDelay Double,PVal Double, "
								"AVal Double,CRdTmE Double, hWkTmOns Double,hWkTmOffs Double,WkTmOns Double,"
								"WkTmOffs Double,WalkCost Double,RideCost Double, OperCost Double, TCost Double,"
								"XC Double,YC Double ,ZC Double , Geometry Point ); ";
							if (blnHistoric) {
								tblTrip = InpFiles.get_tblstop() + "D" + to_string<short> (D) + "hi" + to_string<int> (q)+ "t" + to_string<long> (kStop1.tripId()) + "W" + to_string<int> ((int) (gcost->get_walkcost()*10))+ "R" + to_string<int> ((int) (gcost->get_ridecost()*10)); 
							} else 
							{
								tblTrip = InpFiles.get_tblstop() + "D" + to_string<short> (D) + "ai" + to_string<int> (q) + "t" + to_string<long> (kStop1.tripId())+ "W" + to_string<int> ((int) (gcost->get_walkcost()*10))+ "R" + to_string<int> ((int) (gcost->get_ridecost()*10)); 
							}
							replace(tblTrip.begin(),tblTrip.end(),' ','_');
							ReplaceAll2(tblTrip,"__","_");

							blnCreate = createSpaTbl(str1,tblTrip,pTransOutDb,intSRID,outfile);
							if (!blnCreate) {
							// exit since the DP Result table could not be created 
								logfile << "Could not create the DP result table" <<tblTrip<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
								cout << "Could not create the DP result table" <<tblTrip<<" . Exiting !"<<endl;
								cout << "Press enter to exit!"<<endl;
								cin>>strOutFileName;
								exit (0);
								
							}	
							strSQL = " SELECT recovergeometrycolumn('" + tblTrip + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Point',2);";
							blnCreate = recoverSpatialGeometry(strSQL,tblTrip,pTransOutDb,intSRID,outfile);
							stopout<<"Scid"<<","<<"M"<<","<<"j-i"<<","<<"k-j"<<","<<"o"<<",DPattern";

							 // 3d (i,j,k) then make DP pattern 
							for (i=z;i>=1;i--)
							{
								for (j=z+1;j<=i+M && j<=n;j++)
								{
									for (k=j+1;k<=j+M && k<=n;k++)
									{
												q++; //pattern serial id
												// first insert the i,j,k,l,m stops into the new set - tsidpmap 
												tsidmit = tsidmap.find(i);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(i,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<i;
												}
												tsidmit = tsidmap.find(j);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(j,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<j;
												}
												tsidmit = tsidmap.find(k);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(k,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<k;
												}
												// move forward until m = n according to the pattern by
												// repeating the i,j,k pattern itself  
												o=p=k;
												while (o<n)
												{
													o=o+j-i;
													if (o>=n) {
														tsidmit = tsidmap.find(n);
														tsidpmap.insert(tsidpair (n,tsidmit->second));
														break;
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o+k-j;
													if (o>=n) {
														tsidmit = tsidmap.find(n);
														tsidpmap.insert(tsidpair (n,tsidmit->second));
														break;
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
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
														tsidmit = tsidmap.find(p);
														tsidpmap.insert(tsidpair (p,tsidmit->second));
														break;}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o-j+i;
													if (o<=1) {
														p=1;
														tsidmit = tsidmap.find(p);
														tsidpmap.insert(tsidpair (p,tsidmit->second));
														break;
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
												}

												stopout<<q<<","<<M<<","<<j-i<<","<<k-j<<","<<o<<",";
												string txthdr = ""; 
												writeTextDPatternData(tsidpmap,stop0,stopout,txthdr);
												blnHistoric = false;
										
										if (!blnEuclid) {
												//fileName(tstopiname,"_Orig.bin",outfilename,-1);
												strStopFileName = strStopBaseName + "_Orig.bin";
												mmapSver2.clear();
												ifstream originin(strStopFileName.c_str(),  ios::in | ios::binary);
												readBindblong(mmapSver2,originin);
												if (originin.is_open()) {
													originin.close();
												}
												originin.clear();
												// get the vertex assignment data
												onoff = 0;
												mmapVx0.clear();	
												//fileName(vertinfname,"Vx.bin",outfilename,onoff);
												strVertexFileName = strVertexBaseName + "_" + to_string<int>(onoff) + "Vx.bin";
												originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
												readBindblong(mmapVx0,originin);
												if (originin.is_open()) {
													originin.close();
												}
												originin.clear();  

												onoff = 1;
												mmapVx1.clear();	
												//fileName(vertinfname,"Vx.bin",outfilename,onoff);
												strVertexFileName = strVertexBaseName + "_" + to_string<int>(onoff) + "Vx.bin";
												originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
												readBindblong(mmapVx1,originin);
												if (originin.is_open()) {
													originin.close();
												}
												originin.clear();  

												// transfer the vertex network to the two on and off multi-maps
												// read the vertex object data as a binary file
												infile.clear();
												//fileName(vertinfname,".bin",outvertfilename);
												strVertexFileName = strVertexBaseName + ".bin";
												infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
												mapvert0.clear();
												readBinVertexObjectData(mapvert0,vx,infile);
												// close and open vertex data for the second time
												if (infile.is_open()) {
													infile.close();
												}
												infile.clear();
												infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
												mapvert1.clear();
												readBinVertexObjectData(mapvert1,vx,infile);
												if (infile.is_open()) {
													infile.close();
												}
												infile.clear();
												// mmapVx0,1 - the vertex origin collection
												// mmapVxStop - the vertex - stop multi-map 
												// mmapStopVx - the stop - vertex multi-map 
												// mmapSver2 - the sorted walktime - vertex id multi-map 
												// mmapVerlng - vertex id with cost multi-map 
										}
										tsidpmap = stopspacing_analysis (mmapV1V2EOId, mmapV1V2EgId, mmaparced, mapvert, maped,mmaparStop,
											mmapSver2, mmapVxStop,mmapStopVx, mmapVx0, mmapVx1,mmapVx0R, mmapVx1R, blnHistoric,blnEuclid, 
											strEdgeBaseName, strVertexBaseName, strParcelBaseName, strStopBaseName,gc1, *phway, *luci,tsidpmap, q, logfile,
											InpFiles, kStop1,pTransOutDb, strADB,xi);
					
										// store the stop results for the dp routine
										sdptsidpmmap = insert3DPStops(tsidpmap,dptsidpmap,sdptsidpmmap,stop0,dpstop,q,i,str1);
										tsidpmap.clear();
										cout <<endl<<"Run no. "<<q<<" is completed."<<endl;

									} //k ranging from k up to j + M
								} //j ranging from z+1 up to i + M	 
							} //i ranging from z down to z-M 	 

						}  // if 3d state space
						 // if 5d (i,j,k,l,m ) make pattern accordingly for DP run
						if ( D==5)
						{
							// five state space output
							stopout<<"Scid"<<","<<"MaxSkip"<<","<<"j-i"<<","<<"k-j"<<","<<"l-k"<<","<<"m-l"<<","<<"o"<<",DPattern"<<endl;
							for (i=z;i>=1;i--)
							{
								for (j=z+1;j<=i+M && j<=n;j++)
								{
									for (k=j+1;k<=j+M && k<=n;k++)
									{
										for (l=k+1;l<=k+M && l<=n;l++)
										{
											for (m=l+1;m<=l+M && m<=n;m++)
											{
												q++; //pattern serial id
												// check against the list of the stops taht are fixed before generating the scenario

												// first insert the i,j,k,l,m stops into the new set - tsidpmap 
												tsidmit = tsidmap.find(i);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(i,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<i;
												}
												tsidmit = tsidmap.find(j);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(j,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<j;
												}
												tsidmit = tsidmap.find(k);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(k,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<k;
												}
												tsidmit = tsidmap.find(l);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(l,tsidmit->second));
												} else {
													outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<l;
												}
												tsidmit = tsidmap.find(m);
												if (tsidmit!=tsidmap.end()) {
													tsidpmap.insert(tsidpair(m,tsidmit->second));
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
														tsidmit = tsidmap.find(n);
														if (tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (n,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<n;
															break;
														}
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o+k-j;
													if (o>=n) {
														tsidmit = tsidmap.find(n);
														if(tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (n,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
															break;
														}

													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o+l-k;
													if (o>=n) {
														tsidmit = tsidmap.find(n);
														if (tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (n,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<n;
															break;
														}
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o+m-l;
													if (o>=n) {
														tsidmit = tsidmap.find(n);
														if (tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (n,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
															break;
														}
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
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
														tsidmit = tsidmap.find(p);
														if (tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (p,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
															break;
														}
													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o-l+k;
													if (o<=1) {
														p=1;
														tsidmit = tsidmap.find(p);
														tsidpmap.insert(tsidpair (p,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
														break;
													}
												}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o-k+j;
													if (o<=1) {
														p=1;
														tsidmit = tsidmap.find(p);
														if (tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (p,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<p;
															break;
														}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
													o=o-j+i;
													if (o<=1) {
														p=1;
														tsidmit = tsidmap.find(p);
														if (tsidmit!=tsidmap.end()) {
															tsidpmap.insert(tsidpair (p,tsidmit->second));
														} else {
															outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
															break;
														}

													}
													tsidmit = tsidmap.find(o);
													if (tsidmit!=tsidmap.end()) {
														tsidpmap.insert(tsidpair (o,tsidmit->second));
													} else {
														outfile<<q<<"\t"<<"Scenario lacks "<<"\t"<<"Stop Id"<<"\t"<<o;
													}
												}

												stopout<<q<<","<<M<<","<<j-i<<","<<k-j<<","<<l-k<<","<<m-l<<","<<o<<",";
												string txthdr = ""; 
												writeTextDPatternData(tsidpmap,stop0,stopout,txthdr);
												blnHistoric = false;
										
												if (!blnEuclid) {
												//fileName(tstopiname,"_Orig.bin",outfilename,-1);
												strStopFileName = strStopBaseName + "_Orig.bin";
												mmapSver2.clear();
												ifstream originin(strStopFileName.c_str(),  ios::in | ios::binary);
												readBindblong(mmapSver2,originin);
												if (originin.is_open()) {
													originin.close();
												}
												originin.clear();
												// get the vertex assignment data
												onoff = 0;
												mmapVx0.clear();	
												//fileName(vertinfname,"Vx.bin",outfilename,onoff);
												strVertexFileName = strVertexBaseName + "_" + to_string<long>(onoff) +"Vx.bin";
												originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
												readBindblong(mmapVx0,originin);
												if (originin.is_open()) {
													originin.close();
												}
												originin.clear();  

												onoff = 1;
												mmapVx1.clear();	
												//fileName(vertinfname,"Vx.bin",outfilename,onoff);
												strVertexFileName = strVertexBaseName + "_" + to_string<long>(onoff) +"Vx.bin";
												originin.open(strVertexFileName.c_str(),  ios::in | ios::binary);
												readBindblong(mmapVx1,originin);
												if (originin.is_open()) {
													originin.close();
												}
												originin.clear();  

												// transfer the vertex network to the two on and off multi-maps
												// read the vertex object data as a binary file
												infile.clear();
												//fileName(vertinfname,".bin",outvertfilename);
												strVertexFileName = strVertexBaseName +".bin";
												infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
												mapvert0.clear();
												readBinVertexObjectData(mapvert0,vx,infile);
												// close and open vertex data for the second time
												if (infile.is_open()) {
													infile.close();
												}
												infile.clear();
												infile.open(strVertexFileName.c_str(), ios::in | ios::binary);
												mapvert1.clear();
												readBinVertexObjectData(mapvert1,vx,infile);
												if (infile.is_open()) {
													infile.close();
												}
												infile.clear();
												// mmapVx0,1 - the vertex origin collection
												// mmapVxStop - the vertex - stop multi-map 
												// mmapStopVx - the stop - vertex multi-map 
												// mmapSver2 - the sorted walktime - vertex id multi-map 
												// mmapVerlng - vertex id with cost multi-map 
												}
												phw1 = (phsmap.find(toString<long>(kStop1.tripId()))->second); //kStop1.tripPeriod()))->second; 
												tsidpmap = stopspacing_analysis (mmapV1V2EOId, mmapV1V2EgId, mmaparced, mapvert, maped,mmaparStop,
													mmapSver2, mmapVxStop,mmapStopVx, mmapVx0, mmapVx1,mmapVx0R, mmapVx1R, blnHistoric,blnEuclid, 
													strEdgeBaseName, strVertexBaseName, strParcelBaseName, strStopBaseName,	gc1, *phway, *luci,tsidpmap, q, logfile,
													InpFiles,kStop1,pTransOutDb,strADB,xi);


												// store the stop results for the dp routine
												sdptsidpmmap = insert5DPStops(tsidpmap,dptsidpmap,sdptsidpmmap,stop0,dpstop,q,i,strDPSQLCreateTbl);
												// store the stop results in a sqlite table 

												//str1 =	"( Id Integer PRIMARY KEY, TripId integer,Walk Integer, Ride Integer, " 
												//	" I integer, J integer, K integer, L integer, M integer, StopId Integer,"
												//	" StopName Text, StopIdp Integer,StopIds Integer,EdgeId Integer,"
												//	" palong Double,StOrdr Integer, StopLbl Text,lbl Integer," 
												//	' Hist Integer ,Inbound Integer , External Integer , "
												//	" Included Integer ,Eliminate Integer , CumDist Double, "
												//	" CRdTm Double,CrdTmC Double,HistOns Double,HistOffs Double, "
												//	"Ons Double,Offs Double,DepVol Double,probStoph Double,probStop Double,"
												//	"depDelay Double, arrDelay Double,dwlDelay Double,rideDelay Double,PVal Double, "
												//	"AVal Double,CRdTmE Double, hWkTmOns Double,hWkTmOffs Double,WkTmOns Double,"
												//	"WkTmOffs Double,WalkCost Double,RideCost Double, OperCost Double, TCost Double,"
												//	"XC Double,YC Double ,ZC Double , Geometry Point ); ";
												//  intSRID = InpFiles.get_srid();
												// tblTrip = InpFiles.get_tblstop() + "D" + to_string<short> (D) + "t" + to_string<long> (kStop1.tripId())+ "W" + to_string<int> ((int) (gcost->get_walkcost()*10))+ "R" + to_string<int> ((int) (gcost->get_ridecost()*10)); 
												//replace(tblTrip.begin(),tblTrip.end(),' ','_');
												tsidpmap.clear();
												cout <<endl<<"Run no. "<<q<<" is completed."<<endl;
											} // m range l up to l + M
										} //l ranging from l up to k + M
									} //k ranging from k up to j + M
								} //j ranging from z+1 up to i + M	 
							} //i ranging from z down to z-M 	 
							// save the dp in a table of its own
							
						} // if D=5 do the quintuplet dimension run
						// save the dp in a table of its own

						strDPSQLCreateTbl =	"(  TripId integer,Walk Integer, Ride Integer, " 
							"Q integer, I integer, J integer, K integer, L integer, M integer, Id Integer,"
							" StopId Integer, StOrdr Integer, StopName Text, StopIdp Integer,StopIds Integer,"
							" EdgeId Integer, palong Double, lbl Integer," 
							" Hist Integer ,Inbound Integer , External Integer , "
							" Included Integer ,Eliminate Integer , CumDist Double, CRdTm Double,"
							" undCRdTm Double,CrdTmC Double,HistOns Double,HistOffs Double, "
							"Ons Double,Offs Double,DepVol Double,probStoph Double,probStop Double,"
							"depDelay Double, arrDelay Double,stopDelay Double, dwlDelay Double,rideDelay Double,PVal Double, "
							"AVal Double,CRdTmE Double, hWkTmOns Double,hWkTmOffs Double,WkTmOns Double,"
							"WkTmOffs Double,WalkCost Double,RideCost Double, OperCost Double, TCost Double,"
							"XC Double,YC Double , Geometry Point ); ";
						intSRID = InpFiles.get_srid();
						strDPTblName = InpFiles.get_tblstop() + "D" + to_string<short> (D) + "t" + to_string<long> (kStop1.tripId())+ kStop1.tripPeriod() + "W" + to_string<int> ((int) (gcost->get_walkcost()*10))+ "R" + to_string<int> ((int) (gcost->get_ridecost()*10)); 
						replace(strDPTblName.begin(),strDPTblName.end(),' ','_');
						ReplaceAll2(strDPTblName,"__","_");

						blnCreate = createSpaTbl(strDPSQLCreateTbl,strDPTblName,pTransOutDb,intSRID,outfile);
						if (!blnCreate) {
						// exit since the DP Result table could not be created 
							logfile << "Could not create the DP result table" <<strDPTblName<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
							cout << "Could not create the DP result table" <<strDPTblName<<" . Exiting !"<<endl;
							cout << "Press enter to exit!"<<endl;
							cin>>strOutFileName;
							exit (0);
							
						}	

						keyFld = "ID";
						
						strDPSQLInserTblDef =  " (\"TripId\","  
								"\"Walk\",\"Ride\",\"Q\",\"I\",\"J\",\"K\",\"L\",\"M\",\"ID\",\"StopID\"," 
								"\"StOrdr\",\"StopName\",\"StopIdp\",\"StopIds\",\"EdgeId\",\"palong\"," 
								"\"lbl\",\"Hist\",\"Inbound\",\"External\",\"Included\",\"Eliminate\"," 
								"\"CumDist\",\"CRdTm\",\"undCRdTm\",\"CrdTmC\",\"HistOns\",\"HistOffs\", "  
								"\"Ons\",\"Offs\",\"DepVol\",\"probStoph\",\"probStop\",\"depDelay\", "  
								"\"arrDelay\",\"dwlDelay\",\"rideDelay\",\"PVal\",\"AVal\",\"CRdTmE\", "  
								"\"hWkTmOns\",\"hWkTmOffs\",\"WkTmOns\",\"WkTmOffs\",\"WalkCost\",\"RideCost\", "  
								"\"OperCost\",\"TCost\",\"XC\",\"YC\",\"StopDelay\" ) "  
								" Values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);" ;
						
						strDPSQLInserTbl = "Insert into " + strDPTblName + "  " + strDPSQLInserTblDef ;
						
						// call insert table routine to populate the stop table
						//	blnInsert = inSpaTblStop(strSQLInserTbl,pTransOutDb,*pstop,stops,tblTrip, keyFld , intSRID,logfile);
						// make the trip dpmap
						ip = kStop1.tripId();

						tripDPsmap = dpTripWalkRideMap(sdptsidpmmap,tripDPsmap,ip,dpstop,tDPstp,gcost);
						blnS = inSpaTblTripDPStop(strDPSQLInserTbl,pTransOutDb,tDPstp,dpstop,stop0,sdptsidpmmap,gcost,ip,i,logfile);
						if (blnS) {
							strSQL = " Update " + strDPTblName + " set Geometry = MakePoint(XC,YC," + (to_string<int>(intSRID)) + " );";
							blnCreate = execSpatialQuery(strSQL,strDPTblName,pTransOutDb,intSRID,outfile);
					
							if (blnCreate) {

								strSQL = " SELECT recovergeometrycolumn('" + strDPTblName + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Point',2);";
								blnCreate = recoverSpatialGeometry(strSQL,strDPTblName,pTransOutDb,intSRID,outfile);
							}

						} else {
						// exit since the DP Result table could not be created 
							logfile << "Could not Insert into the Trip DP result table" <<strSQLInserTbl<<" . Exiting !"<<endl;
							cout << "Could not Insert into the Trip DP result table" <<strSQLInserTbl<<" . Exiting !"<<endl;
							cout << "Press enter to exit!"<<endl;
							cin>>strOutFileName;
							//exit (0);
						}	

			//			tskidmit++;
			//		} // loop over all the trips

				if (stopout.is_open()) {
					 stopout.close();
				 }
				 stopout.clear();

				 
				 // write the dp result container
				 rec1 = "D" + to_string<short>(D);
				 rec1.append("W" + to_string<double> ((int) (gcost->get_walkcost()*10)));
				 rec1.append("R" + to_string<double> ((int) (gcost->get_ridecost()*10)));
				 rec1.append("t" + to_string<long> (ip) + kStop1.tripPeriod());
				 replace(rec1.begin(),rec1.end(),' ','_');
    			 ReplaceAll2(rec1,"__","_");
				 //str1 =rec1;
				 //str1.append("dpts.txt");
				// fext[0]='\0';
				//strcpy(fext, str1.c_str());

				 //fileName(tstopiname,str1,outfilename);
				 strStopFileName = strStopBaseName  + rec1 +"dpts.txt";
				 stopout.open(strStopFileName.c_str(), ios::out |ios::trunc);
				sdptsidmmit = sdptsidpmmap.begin(); // sdptsidmit = sdptsidpmap.begin();
				/* 	string txthdr =  "Stopid\tStopidp\tStopids\tEdgeid\tpalong\tStOrdr\tStopLbl\tStopName\tlbl\tblnHist\tblnInbd\tblnExtr\tblnIncl\tblnElim\t";
					txthdr.append("CumDist\tCRdTm\tundCRdTm\tCRdTmC\tHistOns\tHistOffs\tHistDepVol\tOns\tOffs\tDepVol\tprobStoph\tprobStop\t");
					txthdr.append("depDelay\tarrDelay\tdwlDelay\trideDelay\tPVal\tAVal\tCRdTmE\thWkTmOns\thWkTmOffs\tWkTmOns\tWkTmOffs\t");
					txthdr.append("WalkCost\tRideCost\tOperCost\tTCost");
					stopout <<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"l"<<"\t"<<"m"<<"\t"<<"dpkey"<<"\t";
					stopout<<q<<"\t"<<txthdr<<endl;	*/
				dpstop.serializetexthdr(stopout);
				while (sdptsidmmit!=sdptsidpmmap.end())
				{// store the i-j-k-l-m data in the dp store
					(sdptsidmmit->second).serializetext(stopout);
					sdptsidmmit++;
				}

				 if (stopout.is_open()) {
					 stopout.close();
				 }
				 stopout.clear();

				 //replace(str1.begin(),str1.end(),' ','_');
				 //str1.append("dpts.bin");
				 //fext[0]='\0';
				 //strcpy(fext, str1.c_str());

				 //fileName(tstopiname,str1,outfilename);
				 str1 = strStopBaseName  + rec1;
				 strStopFileName = str1 +"dpts.bin";

				 stopout.open(strStopFileName.c_str(), ios::out | ios::binary |ios::trunc);

				 dptStop dpts; 
				 string strkey,strij,strjk;
				 writeBinObjectCol(sdptsidpmmap,dpts,strkey,stopout);
				 if (stopout.is_open()) {
					stopout.close();
				 }
				 stopout.clear();
				 sdptsidmmit = sdptsidpmmap.begin();

				 //at stop 1 , f*(0) = 0
					double icost=0, fcost=0,fcostar=0;
				// run dp routine and trace the optimal path from the above result
					char ijkSep[] = "|";

					size_t j1 = sdptsidpmmap.size();
					if (j1>0) {
						j1=tsidmap.size();
						str1 = datetimeStamp(logfile);
						cout<<str1<<endl;
						if (D==3) {
							cout << " Begin One-Dimensional (i,j,k) Optimzation..."<<endl<<" ..."<<endl;
							logfile << " Begin One Dimensional (i,j,k) Optimzation..."<<endl<<" ..."<<endl;
							sdptsidprmap.clear();
							sdptsidprmap = dpbarun3d(sdptsidprmap,sdptsidpmmap,dpstop,j1,M,strkey,ijkSep,str1);
							cout << " End One-Dimensional (i,j,k) Optimzation..."<<endl;
							logfile << " End Multi-Dimensional (i,j,k) Optimzation..."<<endl;
						} else if ( D==5)
						{
							cout << " Begin Multi-Dimensional (i,j,k,l,m) Optimzation..."<<endl<<" ...  "<<endl;
							logfile << " Begin Multi-Dimensional (i,j,k,l,m) Optimzation..."<<endl<<" ...  "<<endl;
							strStopFileName = InpFiles.get_filestop();
							sdptsidprmap.clear();
							sdptsidpmap.clear();
							//sdptsidpmap = dpbarun5d(sdptsidpmap,sdptsidpmmap,dpstop,j1,M,strkey,ijkSep,str1,InpFiles,gcost,kStop1,pTransOutDb);
							sdptsidpmap = dpOptimalMultiPd5dbxe(sdptsidpmap,sdptsidpmmap,dpstop,j1,M,strDPTblName,strStopBaseName,ijkSep,InpFiles,gcost,phmap,pTransOutDb,dpstop,stop1,sElim,logfile);
							//sdptsidpmap = dpOptimalMultiPd5db(sdptsidpmap,sdptsidpmmap,dpstop,j1,M,strDPTblName,strStopBaseName,ijkSep,InpFiles,gcost,phmap,pTransOutDb,dpstop,stop1,sElim,logfile);
							// write the DP Run Result Map into a new table 
							strDPTblName.append("ResultTrace"); 
							blnCreate = createSpaTbl(strDPSQLCreateTbl,strDPTblName,pTransOutDb,intSRID,outfile);
							if (!blnCreate) {
							// exit since the DP Result table could not be created 
								logfile << "Could not create the DP result table" <<strDPTblName<<" . Exiting !"<<endl<<"Query " <<str1<<endl;
								cout << "Could not create the DP result table" <<strDPTblName<<" . Exiting !"<<endl;
								cout << "Press enter to exit!"<<endl;
								cin>>strOutFileName;
								exit (0);
								
							}
							strDPSQLInserTbl = "Insert into " + strDPTblName + "  " + strDPSQLInserTblDef ;
							// call insert table routine to populate the DP result table
							ip = kStop1.tripId();
							blnS = inSpaTblTripDPStop(strDPSQLInserTbl,pTransOutDb,tDPstp,dpstop,stop0,sdptsidpmap,gcost,ip,i,logfile);
							if (blnS) {
								strSQL = " Update " + strDPTblName + " set Geometry = MakePoint(XC,YC," + (to_string<int>(intSRID)) + " );";
								blnCreate = execSpatialQuery(strSQL,strDPTblName,pTransOutDb,intSRID,outfile);
						
								if (blnCreate) {
	
									strSQL = " SELECT recovergeometrycolumn('" + strDPTblName + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Point',2);";
									blnCreate = recoverSpatialGeometry(strSQL,strDPTblName,pTransOutDb,intSRID,outfile);
									// create the DP summary view for the dp result
									strSQL = " SELECT count( t1.K) NumStops, t1.Walk, t1.Ride,max( t1.CRdTm) CRDTM, "
										" max( t1.undCRdTm) UnCRDTM, max( t1.CrdTmC) CRDTMC,  sum( t1.Ons) Ons, sum( t1.Offs) offs, "
										" min( t1.DepVol) MinDepVol, max( t1.DepVol) MaxDepVol, "
										" sum(( t1.depDelay + t1.arrDelay)* t1.probStop) StopDelay, " 
										" sum( t1.dwlDelay) DwlDelay, sum( t1.rideDelay) RDelay, sum( t1.WkTmOns) WkTmOns, " 
										" Sum( t1.WkTmOffs) WkTmOffs, sum( t1.WalkCost) WkCost, sum( t1.RideCost) RdCost, "
										" sum( t1.OperCost) OpCost, sum( t1.TCost) TCost, St_Collect( t1.Geometry)  Geometry "  
										" FROM " + strDPTblName + " t1 " ;
										string vwName =  "VW_" + strDPTblName ; 
										blnCreate = createView(strSQL, vwName, pTransOutDb, outfile);
								}
							} else {
							// exit since the DP Result table could not be created 
								logfile << "Could not Insert into the Trip DP result table" <<strSQLInserTbl<<" . Exiting !"<<endl;
								cout << "Could not Insert into the Trip DP result table" <<strSQLInserTbl<<" . Exiting !"<<endl;
								cout << "Press enter to exit!"<<endl;
								cin>>strOutFileName;
								exit (0);
								
							}
							// update the table with geometry data

							dptsidpmap.clear();
							sdptsidpmmap.clear();
							sdptsidprmap.clear();
							sdptsidpmap.clear();
							cout << " End Multi-Dimensional (i,j,k,l,m) Optimzation for "<<kStop1.tripPeriod()<< " , Id "<<kStop1.tripId() <<endl;
							logfile << " End Multi-Dimensional (i,j,k,l,m) Optimization for "<<kStop1.tripPeriod()<< " , Id "<<kStop1.tripId() <<endl;
						}
					}

				} // multi dimensional generation


				 if (stopout.is_open()) {
					 stopout.close();
				 }
				 stopout.clear();
				 if (!blnScene) {
					sqlite3_close(pTransOutDb);
				 }
				fprintf( stderr,"\n\tTrip/Period %s and trip Number %s   table %s ; Route %s ;  and Dirn %s DP run is completed \n", kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str());
				fprintf( txtStream,"\n\tTrip/Period %s and trip Number %s   table %s ; Route %s ;  and Dirn %s DP run is completed \n", kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str());
				   // increment at the begining tskidmit++;
			} // check the period headway table if this trip is included in the set of trips to be run
			else { // trip not included
				fprintf( stderr,"\tTrip/Period %s and trip Number %s is not included (skipped) for trip table %s ; Route %s ;  and Dirn %s \n", kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str());
				fprintf( txtStream,"\tTrip/Period %s and trip Number %s is not included (skipped) for trip table %s ; Route %s ;  and Dirn %s \n", kStop1.tripPeriod().c_str(), to_string<long>(kStop1.tripNumber()).c_str(),InpFiles.get_tbltrip().c_str(),InpFiles.get_rtename().c_str(),InpFiles.get_dirname().c_str());

			}
		 } // if the period is found in the period data collection 
	} // loop over all the trips/periods
	str1 = datetimeStamp(logfile);
	cout<<str1<<endl;


	cout << " End of Bustop Spacing Analysis..."<<endl<<" ..."<<endl;
	//cout<<" Press any key and the Enter Key to exit !"<<endl;
	logfile << " End of Bustop Spacing Analysis..."<<endl<<" ..."<<endl;
	//cin>>ch;
				logfileName = (strStopBaseName + "_opt.log"); 
			txtStreamName = (strStopBaseName + "_opt.out");
			txtStdoutName = (strStopBaseName + "_std.out");
			txtStderrName = (strStopBaseName + "_err.out");
	if( txtStream)
	{
	  if ( fclose( txtStream ) )
	  {
		 printf( "The file %s was not closed\n",txtStreamName.c_str() );
	  }
	}

   // Close ll other files 
	int numclosed;

	numclosed = _fcloseall( );
	printf( "Number of files closed : %u\n", numclosed );
	cout << " End of Bustop Spacing Analysis..."<<endl<<" ..."<<endl;
	//cout<<" Press any key and the Enter Key to exit !"<<endl;

	exit(0);
}// end main





 void strAppend( char string[], char suffix[], size_t n )
{
   strncat( string, suffix, __min( n, MaxStrLen-strlen(string)) );
}

 edgev& readedge1(string& rec1, edgev& evx, map<int , string>::iterator maphdrit, globcost& gc, char *seps ) //*seps = "," 
{
	int i=0; int ib=10; int j=0;
    static long sid=0; // serial id
	//	         class vertexp *vertx[2];
	unsigned long eid;  // edge id corresponding to the street object id
	         long esid; // edge serial id
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
        	string frfc; //Vertex From FC pointer
        	string tofc; //Vertex To FC pointer
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
   
  mapflds1 = recoread(rec1,seps);
  mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;
	f1 = maphdrit->second;
    std::transform(f1.begin(), f1.end(), f1.begin(), to_upper());

//						
//ORIGINFC		ORIGINNAME				
//	
	if (f1 == "EDGEOID" || f1 == "EDGEID" || f1 == "OBJECTID" )  // edge id
	{  
		eid = fromString<long>(fldval);   //strtol(strx,&stop1,ib);
	  if (eid>0) {
		  evx.set_id(eid);
		  evx.set_esid(++sid);
	  }
	  else
	  {
		  return evx;
	  }
	}
	else if (f1 == "EDGEFC") // edge feature class
	{
		evx.set_efc(efc);
	}
	else if (f1 == "NDSTARTOID" || f1 == "VX1") // vertex from id 
	{
		frid = fromString<long>(fldval); //strtol(strx,&stop1,ib);
		evx.set_frid(frid);
	} 
	else if (f1 == "NDENDOID" || f1 == "VX2" ) // vertex to id
       {
			toid = fromString<long>(fldval); //strtol(strx,&stop1,ib);
			evx.set_toid(toid);
	  }
	else if (f1 == "EDGEWT" || f1 == "COST") // edge weight(cost)
    {
		ecost = fromString<double>(fldval);   //strtod(strx,&stop1);
	    //evx.set_cost(ecost);
	   evx.set_cost(gc.get_walkcost()/gc.get_ridecost()*ecost);
	}
	else if (f1 == "POSALONG" || f1 == "PALONG") // position along
	{
		palong = fromString<double>(fldval);   //strtod(strx,&stop1);
		evx.set_palong(palong);
	}
	else if (f1 == "ORIGINOID") // origin id
	{
		orig=fromString<long>(fldval);  //strtol(strx,&stop1,ib);
	    evx.set_orig(orig);
	}
	else if (f1 == "EDGESTARTC" || f1 == "EDGESTARTCOST") // edge start cost
	{
	   scost = fromString<double>(fldval); //strtod(strx,&stop1);
	   evx.set_scost(scost);
	}
	else if (f1 == "EDGEENDCOS" || f1 == "EDGEENDCOST") // edge end cost
	{
		tcost = fromString<double>(fldval); //strtod(strx,&stop1);
		evx.set_tcost(tcost);
	}
	else if (f1 == "LABELED") // label
	{
		evx.set_lbl(lbl=fromString<short>(fldval)); //strtol(strx,&stop1,ib));
	}

	else if (f1 == "EDGENOTE") //Edge Note
	{
		enote = fldval;
	   evx.set_enote(enote);
	}
	else if (f1 == "SHAPE_LENG" || f1 == "SHAPE_LENGTH" || f1 == "LENGTH" )	// shape length
	{
		 slen=fromString<double>(fldval); //strtod(strx,&stop1);
		 evx.set_slen(slen);
	}            
   // get the next token
	mapfldsit++;
	maphdrit++;
}
   if (evx.get_orig()>0)
   {
	cout<<"Stop Edge id "<<eid<<", orig : "<<orig<<", Cost "<<ecost<<", i "<<frid<<", j "<<toid<<" Ride Time= "<<evx.get_scost()<<endl;
   }

   return evx;
}
/*
parcel& readparc1(string& rec1, parcel& parc, map<int , string>::iterator maphdrit, char *seps ) // *seps = "\t" 
{
	  static int rno = 0; //record number

	int i=0; int ib=10; int j=0;
//	         class edgev *edg[2];  //edge pointer for ons & offs
	unsigned long pid;  // parcel id long
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
         	 double palong; // total path cost
//	         long evp; // edge pointer id
//			 long frid;  // vertex from id
//			 long toid;  // vertex to id
//			 long origon;  // origin for ons of parc id
//			 long origoff;  // origin for offs of parc id
//			 long ons;  // ons strength
//			 long offs;  // offs strength
        	string ptype; // parcel land use type 
        	string pfc; // parcel FC 
        	string frfc; //Vertex From FC 
        	string tofc; //Vertex To FC 
        	string orfc; //Origin FC 
        	string ornm; //Origin Name 
        	string pnote; //parcel Note 
//        	parcel* par1;

// PARCELOID,PIDLONG,PTYPE,LANDSF,GROSSAREA,LIVINGAREA,TOTALVAL,LANDVAL,
// BLDGVAL,COMPFACT,PVAL,AVAL,INONS,INOFFS,OUTONS,OUTOFFS,EDGEID,POSALONG

   const int bufsz = 1000; // Buffer size;

   string f1, fldval,fldhdr;
   typedef map<int, string, less<int>> mapflds;
// vertex map iterator objects
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;

   rno++;

   mapflds1 = recoread(rec1,seps);
   mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;
	f1 = maphdrit->second;
    std::transform(f1.begin(), f1.end(), f1.begin(), to_upper());
//	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "PARCELOID" || f1 == "OBJECTID" ) // parcel object id
	{
		pid = fromString<long>(fldval); //strtol(val1,&stop1,ib);
	  if (pid>0) {
		  parc.set_id(pid);
	  }
	  else
	  {
		  return parc;
	  }
	}
	else if (f1 == "PID_LONG" || f1 == "PIDLONG" || f1 == "PID" || f1 == "PARCELID") // parcel String Id
	{
		parc.set_pacid(fldval);
	  }
	else if (f1 == "PTYPE") // parcel land use type - PTYPE
	{
		parc.set_ptype(fldval);
	  }
	else if (f1 == "LU") // parcel land use type - PTYPE
	{
		parc.set_lucd(fldval);
	  }
	else if (f1 == "LANDSF" || f1 == "LAND_SF") // LAND_SF
	{ 
		dblVal = fromString<double>(fldval); //dblVal = strtod(val1,&stop1);
		parc.set_lndsf(dblVal);
	} 
	else if (f1 == "GROSSAREA" || f1 == "GROSS_AREA") // Gross Area - 
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	   parc.set_grarea(dblVal);
	}
	else if (f1 == "LIVINGAREA" || f1 == "LIVING_ARE") // LIVING_ARE
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      parc.set_lvarea(dblVal);
	  }
	else if (f1 == "TOTALVAL" || f1 == "FY2003_TOT") // Total valuation
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		parc.set_valtot(dblVal);
	}
	else if (f1 == "LANDVAL" || f1 == "FY2003_LAN") // Land valuation
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      parc.set_valand(dblVal);
	  }
	else if (f1 == "BLDGVAL" || f1 == "FY2003_BLD") // Bldg valuation
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		parc.set_valbld(dblVal);
	  }
	else if (f1 == "COMPFACT") // competition factor
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      parc.set_cfact(dblVal);
		}
	else if (f1 == "PVAL") // Production strength
	{	  dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      parc.set_pval(dblVal);
	  }
	else if (f1 == "AVAL") // parcel attraction strength
	{	  dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
	      parc.set_aval(dblVal);
	  }
	else if (f1 == "HINONS" || f1 == "HONS" || f1 == "HOUTONS") // in Boardings - inOns
	{		
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		parc.set_ons(0.0);
		}
	else if (f1 == "HINOFFS" || f1 == "HOFFS" || f1 == "HOUTOFFS")  // Aligtings - inOffs 
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
			parc.set_offs(0.0);
		}
	else if (f1 == "INONS" || f1 == "ONS" || f1 == "OUTONS") // in Boardings - inOns
	{		
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
		parc.set_ons(0.0);
		}
	else if (f1 == "INOFFS" || f1 == "OFFS" || f1 == "OUTOFFS")  // Aligtings - inOffs 
	{
		dblVal = fromString<double>(fldval); //strtod(val1,&stop1);
			parc.set_offs(0.0);
		}
	else if (f1 == "EDGEID") // Edge Object Id  
	{	
		eoid=fromString<long>(fldval); //strtol(val1,&stop1,ib);
	    parc.set_eoid(eoid);
	  }
	else if (f1 == "PALONG" || f1 == "POSALONG") // position along edge
	{	
		palong = fromString<double>(fldval); //strtod(val1,&stop1);
		parc.set_palong(palong);
	}
   // get the next token
	mapfldsit++;
	maphdrit++;
}
   if (parc.get_origon()>0)
   {
	// cout<<"Stop parcel id "<<pid<<", orig : "<<orig<<", Cost "<<cost<<", i "<<frid<<", j "<<toid<<endl;
   }

   return parc;
}
*/
map<int , string> readatahdr(string& rec1, char *seps) // *seps = "," 
{
	int i=0; int ib=10; int j=0;


   string f1, fldval;
   typedef map< int, string > mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int ,string  >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
  bool done = false;

  while(!done) 
  {
      size_t cPos = rec1.find(seps);
	  if(cPos == string::npos) {
     // last piece of record
		  i++;
		  fldval = rec1.substr(0,cPos);
         std::transform(fldval.begin(), fldval.end(), fldval.begin(), to_upper());
		 mapflds1.insert(fld_pair(i,fldval));
		 done=true;
	  }
	  else
	  {
		  fldval = rec1.substr(0,cPos);
         std::transform(fldval.begin(), fldval.end(), fldval.begin(), to_upper());
		  i++;
		  mapflds1.insert(fld_pair(i,fldval));
          rec1.erase(0,cPos+1);
	  }
  }
  return mapflds1;
 }

globcost& readgcost1(string& rec1, globcost& gc1, map<int , string>::iterator maphdrit, char *seps ) // *seps = "," 
{
   static int rno = 0; //record number

   int i=0; int ib=10;

//COSTWALK,COSTRIDE,UNITONTM,UNITOFFTM,MAXWLKDIST,PROPENSITY,FILESTEM,NOPERIODS,
//WALKSPD,FILEPATH
			 double walkcost;  // vertex id
	         short nopds; // no of periods
         	 float unitontm;  // unit boarding time
             float unitofftm;  // unit alighting time
			 float maxwalkdist; // total cost
			 double ridecost;  // ride cost 
			 float propensity;  // propensity ratio
	         float walkspd; // walk speed 
	         string  filestem; // stem for file naming 
	         string  filepath; // path name for output files
			 int dpdimension;

   string f1, fldval,fldhdr;
   typedef map<int, string, less<int>> mapflds;
   mapflds :: iterator fld1map_Iter, fld2map_Iter;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
   mapflds::iterator mapfldsit;
 
   rno++;

  mapflds1 = recoread(rec1,seps);
  mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;

	fldhdr = maphdrit->second;


     std::transform(fldhdr.begin(), fldhdr.end(), fldhdr.begin(), to_upper());
	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "COSTWALK")   // Walk Cost
	{
		walkcost = fromString<double>(fldval);//val1,&stop1);
	  if (walkcost>0) {
		  gc1.set_walkcost(walkcost);
	  }
	  else
	  {
          gc1.set_walkcost(0);	 
          gc1.set_nopds(0);	 
		  return gc1;
	  }
	}
	else if (f1 == "COSTRIDE") // Ride Cost
	{
      ridecost = fromString<double>(fldval); //strtod(val1,&stop1);
	  if (walkcost>0) {
		  gc1.set_ridecost(ridecost);
	  }
	}
	else if (f1 == "UNITONTM") // Unit On Time
	{
		unitontm = fromString<float>(fldval);
		  gc1.set_unitontm(unitontm);
	}
	else if (f1 == "UNITOFFTM") // Unit Off Time
	{
		unitofftm = fromString<float>(fldval); //strtod(val1,&stop1);
		  gc1.set_unitofftm(unitofftm);
	}
	else if (f1 == "MAXWLKDIST") // Max Walk Distance
	{
			maxwalkdist = fromString<float>(fldval); //strtod(val1,&stop1);
			gc1.set_maxwalkdist(maxwalkdist);
	}
	else if (f1 == "PROPENSITY") // propensity
	{
		propensity = fromString<float>(fldval); //strtod(val1,&stop1);
		gc1.set_propensity(propensity); 
	}
	else if (f1 == "FILESTEM") // FILE STEM used for naming output files
	{
		    gc1.set_filestem(filestem);
	}
	else if (f1 == "NOPERIODS" || f1 == "PERIOD") // Analysis Period used to extract the period parameters from the Period/headway table
	{
		nopds = fromString<short>(fldval);  //strtol(val1,&stop1,ib);
			gc1.set_nopds(nopds); 
	}
	else if (f1 == "WALKSPD") // WALKing SPeeD
	{
			walkspd = fromString<float>(fldval); //strtod(val1,&stop1);
		    gc1.set_walkspd(walkspd);
	}
	else if (f1 == "FILEPATH") // FILE path for saving output files
	{
		    gc1.set_filepath(filepath);
	}
	else if (f1 == "DPDIMENSION") // FILE path for saving output files
	{
			dpdimension = fromString<float>(fldval); //strtod(val1,&stop1);
		    gc1.set_dpdimension(dpdimension);
	}
   // get the next fld
	mapfldsit++;
	maphdrit++;
}
   if (gc1.get_walkcost()>0)
   {
	   gc1.show_globcosthdr(cout);
	   gc1.show_globcost(cout);
   }

   return gc1;
}

 pdhdway& readpdhdway1(string& rec1, pdhdway& pdh1, map<int , string>::iterator maphdrit, char *seps ) // *seps = "," 
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

  mapflds1 = recoread(rec1,seps);
 
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
		  pdh1.set_pdId(pdId);
	  }
	  else
	  {
		  pdh1.set_pdId(0);
		  return pdh1;
	  }
	}
	if (f1 == "FLDBEGT") // Period Begining Time
	{
		begtm = fromString<float>(fldval); //strtod(val1,&stop1);
		  pdh1.set_begtm(begtm);
	  }
	if (f1 == "FLDHDWAY") // headway
	{
		hdway = fromString<float>(fldval); //strtod(val1,&stop1);
		  pdh1.set_hdway(hdway);
	  }
	if (f1 == "FLDPDLEN") // period length
	{
		pdlen = fromString<float>(fldval); //strtod(val1,&stop1);
		  pdh1.set_pdlen(pdlen);
	  }
	if (f1 == "COSTOPER") // operating cost
	{
		opercost = fromString<float>(fldval); //strtod(val1,&stop1);
		pdh1.set_opercost(opercost);
	} 
   // get the next fld
	mapfldsit++;
	maphdrit++;
}

   return pdh1;
}


 lucodes& readlucodes1x(string& rec1, lucodes& luc1, map<int , string>::iterator maphdrit, char *seps ) // *seps = "," 
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
  mapflds1 = recoread(rec1,seps);

  mapfldsit = mapflds1.begin();
while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;
	fldval = mapfldsit->second;

	fldhdr = maphdrit->second;
     std::transform(fldhdr.begin(), fldhdr.end(), fldhdr.begin(), to_upper());
	 f1 = (fldhdr);

	if (f1 == "PTYPE") // Property Type
	{  
	  ci = fromString<long>(fldval); // strtol(val1,&stop1,ib);
	  if (ci>0) {
		  luc1.set_pType(fldval);
	  }
	  else
	  {
		  luc1.set_pType("0");
		  return luc1;
	  }
	}
	else if (f1 == "DESC") // Description
	{  Desc = fldval;
		  luc1.set_Desc(Desc);
	}
	else if (f1 == "LU") // Land Use Code
	{
		LUC = fldval;
		luc1.set_LUC(LUC);
	  }
	else if (f1 == "KEYFLD" )   // Property Type Text Code
	{
		KeyProp = fldval;
		luc1.set_KeyProp(KeyProp);
	}
	else if (f1 == "ONCOEF") // Unit On rate per key field
	{  
		OnCoeff = fromString<float>(fldval); //strtod(val1,&stop1);
		luc1.set_OnCoeff(OnCoeff);
	}
	else if (f1 == "OFFCOEF") // Unit Off rate per key field
	{
		OffCoeff = fromString<float>(fldval); //strtod(val1,&stop1);
		luc1.set_OffCoeff(OffCoeff);
	} 
   // get the next fld
	mapfldsit++;
	maphdrit++;
}
   if (fromString<long>(luc1.get_pType())>0)
   {
	   //luc1.show_lucodeshdr(cout);
	   //luc1.show_lucodes(cout);
   }

   return luc1;
}

map<int , string> recoread(string& rec1, char *seps ) // *seps = "," 
{
	int i=0; 
    string fldval;
   typedef map<int, string, less<int>> mapflds;
   typedef pair < int, string >  fld_pair;
   mapflds mapflds1;
 
  bool done = false;

  while(!done) 
  {
      size_t cPos = rec1.find(seps);
	  if(cPos == string::npos) {
     // last piece of record
		  i++;
		  fldval = rec1.substr(0,cPos);
		 mapflds1.insert(fld_pair(i,fldval));
		  done=true;
	  }
	  else
	  {
		  fldval = rec1.substr(0,cPos);
		  i++;
		  mapflds1.insert(fld_pair(i,fldval));
          rec1.erase(0,cPos+1);
	  }
  }
  return mapflds1;
}

maplngvx& vxv (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,mmaplng& mmapved,
				  mmaplng& mmapVxStop,short onoff) 
{
// vertex voronoi analysis
	// variables
	double ecost,vcost,tcost,xcost=0;
    long k=0, vid=0,eid=0,eoid=0,ip=0,pid=0,origin=0,v1=0,v2=0,vcnt=0,o1=0; // vid - vertex id, eid - edge id, pid - parcel id
    vertexp* pVx1,*pVx2;
    edgev* pEg;
	bool bnInf=true,forwd=true,rev=true,bnVxUp=false,bnEdge=false,bnSameOrg=true;
	size_t j1=0;
/* start the search from the sorted vertex origin collection  
 and process the vertex voronoi partition */
mmapSverit = mmapSver.begin();
ip = 0;
while(mmapSverit != mmapSver.end()) 
{  
	  // get the minimum of the origin list
    vcost =(*mmapSverit).first;
    v1 =(*mmapSverit).second;
	mapverit = mapvert.find(v1);
	pVx1 = &(mapverit->second);
	pVx1->set_lbl(-1);
	pVx1->set_cost(vcost);
	if (pVx1->get_idp()==-1) 
	{
		pVx1->set_index(ip);
		pVx1->set_lowlink(ip);
	} else {
		if (pVx1->get_idp()<-1) 
		{
			pVx1->set_index(ip);
			pVx1->set_lowlink(ip);
		}
	}

    // erase/remove the fixed vertex from the sorted vertex list (the heap) 
	mmapSver.erase(mmapSverit);
    outfile<<"Vertex Fixed  "<<v1<<" orig "<<o1<<" cost "<<vcost<<endl;

    // search through the edges coming out of this vertex - vid
	pair<mumved_AIter, mumved_AIter> vedrange = mmapved.equal_range(v1);
	//size_t j = distance(vedrange.first,vedrange.second);
	for (mumavedit = vedrange.first; mumavedit!=vedrange.second;mumavedit++)
	{
		// get the edge  
		eoid = (*mumavedit).second;
		mapedit = maped.find(eoid);
		pEg = &(mapedit->second);
		// find vertex 2
		v2 = pEg->get_toid();
    	forwd=(pEg->get_frid()==v1);  // if the vertex being scanned is the edge's 1st vertex (forward direction)
		ecost = pEg->get_cost();
		mapverit = mapvert.find(v2);
		// calculate vertex cost, dirn and set origin 
		pVx2 = &(mapverit->second);
        outfile << " i |"<<v1<<"| VxCost |"<<pVx1->get_cost()<<" |edge Id |"<<eoid<< " |EgCost |"<<ecost<< " |j |"<<v2<<" |Vxcost |"<<pVx2->get_cost()<<endl;
		tcost = vcost + ecost;
		bnInf=(pVx2->get_cost()>=inf); // if v2 has not been visited
		bnVxUp = tcost < pVx2->get_cost();
		bnEdge= !bnInf && !bnVxUp; 
		
		if (bnInf || bnVxUp) {	// end vertex is not visited i.e.(cost=INF) or it can improve
			pVx2->set_idp(v1);
			pVx2->set_index(pVx1->get_index()+1);
			pVx2->set_lowlink(pVx1->get_lowlink());
			pVx2->set_orig(pVx1->get_orig());
			pVx2->set_cost(tcost);
			pEg->set_orig(pVx1->get_orig());
			upEdgeCostDirn(pEg, pVx1,pVx2,forwd);
			pVx2->show_vertof(outfile);
			if (bnInf) { // end vertex is not visited i.e. if (cost=INF)
				mmapSver.insert(dblng_Pair(tcost,v2)); // insert it into the heap
			}
		} else if (bnEdge) { // boundary edge found
			bnSameOrg = (pVx1->get_orig()==pVx2->get_orig()); // they two v's belong to the same stop
			if (bnSameOrg) {
				pEg->set_orig(pVx1->get_orig());
			} else {
				pEg->set_orig(-2);
			}
			upEdgeCostDirn(pEg, pVx1,pVx2,forwd);
			pEg = calcEdgeBndPos(pEg,pVx1,pVx2,xcost);
		} 

                //pEg->show_edge(outedgefile);
    			outfile<<" eid "<<pEg->get_id()<<" i "<<pEg->get_frid()<<" j "<<
					pEg->get_toid()<<" cost "<<pEg->get_cost()<<" Scost "<<pEg->get_scost()<<endl;
						  pVx2->set_fid(pVx1->get_id());
	}  // end vertex vertex (for current edge) is being scanned/fixed , i.e. vid = pEg->get_toid 
        mmapSverit = mmapSver.begin();  // pick the next lowest cost vertex 
} // while map heap list is present loop over all vertices in the container 

cout << "End Vertex Voronoi for onoff = "<<onoff<<endl;
	mmapvert.clear();
	mmapvert = remapVertId2Orig(mapvert,mmapvert,*pVx1,ip);	

	VxVorAccum VxVAcc;
	deque<VxVorAccum> VxVorSumAcc;
	deque<VxVorAccum>::iterator it;
	j1=mmapVxStop.size();
	k=-10;
	VxVorSumAcc = vertaggbyMaxKey (j1,k,mmapvert,*pVx1,VxVAcc);

// write out the vertex summary results
	cout << "Summary of Vertex attributes by id " <<endl
    << "Id" << "\t"<< "Stop" << "\t" << " id " << "\t" << " cost "<< "\t" << "TotCost "
	 << "\t"<< "ParCnt"<< "\t"<<"index "<< "\t"<<"LowLink "<< endl;
			for (it=VxVorSumAcc.begin();it!=VxVorSumAcc.end();++it)
			{
			    cout << *it;
			}


cout << "No. of unassigned Vertices = "<<o1 << " Assigned = " <<ip<<" for travel direction = "<<onoff<<endl;

return mapvert;

}; //End vxv routine


maplngvx& vxvtarjan (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,mmaplng& mmapved,
				  mmaplng& mmapVxStop,short onoff) 
{
// vertex voronoi analysis
	// variables
	double ecost,vcost,tcost,xcost=0;
    long k=0, vid=0,eid=0,eoid=0,ip=0,pid=0,origin=0,v1=0,v2=0,vcnt=0,o1=0; // vid - vertex id, eid - edge id, pid - parcel id
    vertexp* pVx1,*pVx2;
    edgev* pEg;
	bool bnInf=true,forwd=true,rev=true,bnVxUp=false,bnEdge=false,bnSameOrg=true;
/* start the search from the sorted vertex origin collection  
 and process the vertex voronoi partition */
mmapSverit = mmapSver.begin();
ip = 0;
while(mmapSverit != mmapSver.end()) 
{  
	  // get the minimum of the origin list
    vcost =(*mmapSverit).first;
    v1 =(*mmapSverit).second;
	mapverit = mapvert.find(v1);
	pVx1 = &(mapverit->second);
	pVx1->set_lbl(-1);
	pVx1->set_cost(vcost);
	if (pVx1->get_idp()==-1) 
	{
		pVx1->set_index(ip);
		pVx1->set_lowlink(ip);
	} else {
		if (pVx1->get_idp()<-1) 
		{
			pVx1->set_index(ip);
			pVx1->set_lowlink(ip);
		}
	}

    // erase/remove the fixed vertex from the sorted vertex list (the heap) 
	mmapSver.erase(mmapSverit);
    outfile<<"Vertex Fixed  "<<v1<<" orig "<<o1<<" cost "<<vcost<<endl;

    // search through the edges coming out of this vertex - vid
	pair<mumved_AIter, mumved_AIter> verng = mmapV1V2EgId.equal_range(v1);
	size_t j1 = distance(verng.first,verng.second);
	for (mumavedit = verng.first; mumavedit!=verng.second;mumavedit++)
	{
		// get the edge  
		eoid = (*mumavedit).second;
		mapedi = maped.find(eoid);
		if (mapedi!=maped.end()) {
			pEg = &(mapedi->second);
			// find vertex 2
			v2 = pEg->get_toid();
    		forwd=(pEg->get_frid()==v1);  // if the vertex being scanned is the edge's 1st vertex (forward direction)
			ecost = pEg->get_cost();
			mapverit = mapvert.find(v2);
		// calculate vertex cost, dirn and set origin
			if (mapverit != mapvert.end()) {
				pVx2 = &(mapverit->second);
				outfile << " i |"<<v1<<"| VxCost |"<<pVx1->get_cost()<<" |edge Id |"<<eoid<< " |EgCost |"<<ecost<< " |j |"<<v2<<" |Vxcost |"<<pVx2->get_cost()<<endl;
				tcost = vcost + ecost;
				bnInf=(pVx2->get_cost()>=inf); // if v2 has not been visited
				bnVxUp = tcost < pVx2->get_cost();
				bnEdge= !bnInf && !bnVxUp; 
		
				if (bnInf || bnVxUp) {	// end vertex is not visited i.e.(cost=INF) or it can improve
					pVx2->set_idp(v1);
					pVx2->set_index(pVx1->get_index()+1);
					pVx2->set_lowlink(pVx1->get_lowlink());
					pVx2->set_orig(pVx1->get_orig());
					pVx2->set_cost(tcost);
					pEg->set_orig(pVx1->get_orig());
					upEdgeCostDirn(pEg, pVx1,pVx2,forwd);
					pVx2->show_vertof(outfile);
					if (bnInf) { // end vertex is not visited i.e. if (cost=INF)
						mmapSver.insert(dblng_Pair(tcost,v2)); // insert it into the heap
					}
				} else if (bnEdge) { // boundary edge found
					bnSameOrg = (pVx1->get_orig()==pVx2->get_orig()); // they two v's belong to the same stop
					if (bnSameOrg) {
					pEg->set_orig(pVx1->get_orig());
					} else {
						pEg->set_orig(-2);
					}
					upEdgeCostDirn(pEg, pVx1,pVx2,forwd);
					pEg = calcEdgeBndPos(pEg,pVx1,pVx2,xcost);
				} 

                //pEg->show_edge(outedgefile);
    			outfile<<" eid "<<pEg->get_id()<<" i "<<pEg->get_frid()<<" j "<<
				pEg->get_toid()<<" cost "<<pEg->get_cost()<<" Scost "<<pEg->get_scost()<<endl;
				pVx2->set_fid(pVx1->get_id());
			} // 2nd vertex is not found
			else {
				cout << " j "<<v2<<" is not found for vertex "<<pVx1->get_id()<<" |edge Id |"<<eoid<<endl;
			}
		} // edge is found
		else { // edge is not found
				cout << " Edge "<< eoid<< " is not found for vertex "<<v1<<endl;
		}

	}  // end vertex1 - vertex2 (for current edge) is being scanned/fixed , i.e. vid = pEg->get_toid 
        mmapSverit = mmapSver.begin();  // pick the next lowest cost vertex 
} // while map heap list is present loop over all vertices in the container 

cout << "End Vertex Voronoi for onoff = "<<onoff<<endl;
	mmapvert.clear();
	mmapvert = remapVertId2Orig(mapvert,mmapvert,*pVx1,ip);	

	VxVorAccum VxVAcc;
	deque<VxVorAccum> VxVorSumAcc;
	deque<VxVorAccum>::iterator it;
	size_t j1=mmapVxStop.size();
	k=-10;
	VxVorSumAcc = vertaggbyMaxKey (j1,k,mmapvert,*pVx1,VxVAcc);

// write out the vertex summary results
	cout << "Summary of Vertex attributes by id " <<endl
    << "Id" << "\t"<< "Stop" << "\t" << " id " << "\t" << " cost "<< "\t" << "TotCost "
	 << "\t"<< "ParCnt"<< "\t"<<"index "<< "\t"<<"LowLink "<< endl;
			for (it=VxVorSumAcc.begin();it!=VxVorSumAcc.end();++it)
			{
			    cout << *it;
			}


cout << "No. of unassigned Vertices = "<<o1 << " Assigned = " <<ip<<" for travel direction = "<<onoff<<endl;

return mapvert;

}; //End vxv routine





void edv (maplngvx& mapvert, mmaplnged& mmaped,mmaplng& mmapved,mmaplng& mmapEdgeStop,short onoff)
{
    vertexp* pVx1,*pVx2,*pVx;
    edgev* pEg;
    long o1=0,o2=0, eid=0; //  eid - edge id
    int j=0,ne=0,nv=0,ne1=0,ne2=0,ne3=0,ne4=0; // ne - no of edges, nv - no of vertices,# fixed, #loop, # bound, # free
    string str1;
	double vcost=0,xcost=0,ecost=0,vcost1=0;
	bool fixed,bound,free,free1,free2,forwd,bacwd,loop;
// process the edge voronoi routine now
j=0;

mmapedit = mmaped.begin();
while(mmapedit != mmaped.end()) 
{ 
     eid =(*mmapedit).first;
	 // get the edge
	pEg =  &(mmapedit->second);
	if (pEg->get_lbl()==0) {j++;} // unlabeled edge count
	if (pEg->get_lbl()==-1) {ne++;} // labeled edge count
	if (pEg->get_dirn()==-2) {o1++;} // boundary 	
	if (pEg->get_dirn()==-1) {o2++;} // reverse 	
	mmapedit++;
}
	cout<<" Vertex Voronoi # of unlableled Edges  = "<< j<<endl;
	cout<<" Vertex Voronoi # of lableled Edges  = "<< ne<<endl;
	cout<<" # of boundary Edges  = "<< o1<<endl;
	cout<<" # of reverse Edges  = "<< o2<<endl;
j=0;ne=0;o1=0;o2=0;
mmapedit = mmaped.begin();
while(mmapedit != mmaped.end()) 
{ 
     eid =(*mmapedit).first;

//		pVx2 = *(&pEg->vertices[1]);
	 
	 // get the edge
	pEg =  &(mmapedit->second);
	eid = pEg->get_id();
	mumESit = mmapEdgeStop.find(eid);
	mapverit = mapvert.find(pEg->get_frid());
	pVx1 = &(mapverit->second);
	mapverit = mapvert.find(pEg->get_toid());
	pVx2 = &(mapverit->second);
		o1 = pVx1->get_orig();
		o2 = pVx2->get_orig();
		fixed = o1>0 && o2>0 && o1==o2;
		bound = o1>0 && o2>0 && o1!=o2;
		free1 = o1<=0 && o2>0; 
		free2 = o2<=0 && o1>0;
		free = o1<=0 && o2<=0;
		if (fixed || bound) 
		{ // both ends are fixed including boundary edges
			forwd = pVx1->get_id() == pVx2->get_idp();
			bacwd = pVx1->get_idp() == pVx2->get_id();
			loop = pVx1->get_id() == pVx2->get_id();
				//ceil( ( num * pow( 10,x ) ) - 0.5 ) / pow( 10,x );
				vcost1=	ceil((abs(pVx1->get_cost()-pVx2->get_cost())*pow(10.0,3.0))-0.5)/pow(10.0,3.0);
				ecost = ceil((pEg->get_cost()*pow(10.0,3.0))-0.5)/pow(10.0,3.0);
					ne++;
					if (forwd)  
					{ // if it is a forward labelled arc 
						pEg->set_orig(pVx1->get_orig());
						pEg->set_scost(pVx1->get_cost());
						pEg->set_tcost(pVx2->get_cost());
						pEg->set_enote( "fixed" );
					}
					else if (bacwd)
					{ // if it is a reverse labelled arc 
						pEg->set_orig(pVx2->get_orig());
						pEg->set_scost(pVx2->get_cost());
						pEg->set_tcost(pVx1->get_cost());
						pEg->set_enote( "fixed" );
					}
					else if ( loop) // loop edge where begining & end are same
					{
						pEg->set_orig(pVx1->get_orig());
						pEg->set_scost(pVx1->get_cost());
						pEg->set_tcost(pVx2->get_cost());
						pEg->set_palong(0.5); // travel is made in either directions to half the distance
						str1.erase();
						str1= "loop node " + to_string(pVx1->get_id()) + " on edge " + to_string(pEg->get_id());
						pEg->set_enote( str1 );
						ne1++;
					}
					else
					{ // Vertices are fixed but edge is not 
						if (pVx1->get_cost()<pVx2->get_cost()) {
							pEg->set_orig(pVx1->get_orig());
							pEg->set_scost(pVx1->get_cost());
							pEg->set_tcost(pVx2->get_cost());
						}
					}
			if (bound) // boundary edge - beg. & end vertices are assigned to different stops
			{ // find the indifference point between the vertices 
				pEg->set_orig(-2);
				xcost = (abs(pVx1->get_cost()-pVx2->get_cost())+ pEg->get_cost())/(2*pEg->get_cost());
				pEg->set_palong(xcost);
				str1.erase();
				str1= "Bdry " + to_string(pVx1->get_orig()) + " - " + to_string(pVx2->get_orig());
				pEg->set_enote( str1 );
				pEg->show_edge(outedgebdyfile);
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
		//		pEg->show_edge(cout);
				nv++;
		}

				pEg->show_edge(outfile);
    if ((pEg->get_orig()>0)||(pEg->get_orig()==-2))
	{
		pEg->show_edge(outedgefile);
       // pEg->show_edge(cout);
	}
	mmapedit++;
} // while map origin is present
		cout<<" End of Edge Voronoi Process "<<endl;
		cout<<" # of Fixed Edges  = "<< ne<<endl;
		cout<<" # of Loop Edges =  "<< ne1<<endl;
		cout<<" # of Boundary Edges = "<< ne2<<endl;
		cout<<" # of Free Vertex 1 Edges = "<< ne3<<endl;
		cout<<" # of Free Vertex 2 Edges = "<< ne4<<endl;
		cout<<" # of Free Edges = "<< nv<<endl;
		cout<<" Total # Edges read = "<< (nv+ne+ne1+ne2+ne3+ne4)<<endl;
		cout<<" # of Missed Edges = "<< (distance(mmaped.begin(), mmaped.end())-(nv+ne+ne1+ne2+ne3+ne4))<<endl;
o1=0;
o2=0;
j=0;
ne2=0;
ne=0;
ne3=0;
mmapedit = mmaped.begin();
while(mmapedit != mmaped.end()) 
{ 
     eid =(*mmapedit).first;
	 // get the edge
	pEg =  &(mmapedit->second);
	if (pEg->get_lbl()==0) {
		j++;
		if (pEg->get_orig()>0) {
			pEg->set_lbl(-1);
			ne3++;
		}
		else {
			pVx1 = &(mapvert.find(pEg->get_frid()))->second;
			pVx2 = &(mapvert.find(pEg->get_toid()))->second;
			free2=(pVx1->get_orig()>0 && pVx2->get_orig()<=0);
			free1=(pVx1->get_orig()<=0 && pVx2->get_orig()>0);
			fixed=(pVx1->get_orig()>0 && pVx2->get_orig()>0);
			loop=(pVx1->get_id()==pVx2->get_id());

			if (fixed&&!loop) { // relabel this edge
				if (roundupx((pEg->get_cost()),3)==roundupx(abs((pEg->get_tcost()-pEg->get_scost())),3)) {
					pEg->set_lbl(-1);
				}
				else {
					pEg->set_lbl(-2);
				}
				pEg->set_orig(pVx1->get_orig());
				if (pVx1->get_cost()<=pVx2->get_cost()) {
					pEg->set_dirn(1);
				} else
				{
					pEg->set_dirn(-1);
				}
			}
			else if(free2) {// fixed at vertex 1 only
				pEg->set_lbl(-1);
				pEg->set_orig(pVx1->get_orig());
				pEg->set_scost(pVx1->get_cost());
				pEg->set_tcost(pVx1->get_cost()+pEg->get_cost());
				pEg->set_dirn(1);
			}
			else if(free1) {// fixed at vertex 2 only
				pEg->set_lbl(-1);
				pEg->set_orig(pVx2->get_orig());
				pEg->set_scost(pVx2->get_cost());
				pEg->set_tcost(pVx2->get_cost()+pEg->get_cost());
				pEg->set_dirn(-1);
			}
			else if ( loop) // loop edge where begining & end are same
			{
						if (pVx1->get_idp() > 0) 
						{ pVx = &((mapvert.find(pVx1->get_idp()))->second);
							pEg->set_orig(pVx->get_orig());
							pEg->set_dirn(2);
							pEg->set_scost(pVx1->get_cost());
							pEg->set_tcost(pVx2->get_cost());
							pEg->set_palong(0.5); // travel is made in either directions to half the distance
							str1.erase();
							str1= "loop node " + to_string(pVx1->get_id()) + " on edge " + to_string(pEg->get_id());
							pEg->set_enote( str1 );
							ne1++;
						}
						else  
						{
							cout <<" Unfixed loop Edge id ="<<pEg->get_id()<<" Vx1 ="<<pVx1->get_id()<<" orig 1="<<
							pVx1->get_orig()<<" Vx2 ="<<pVx2->get_id()<<" orig 2="<<pVx2->get_orig()<<endl;
						}
			}
			else {
//				cout <<" Unfixed Edge id ="<<pEg->get_id()<<" Vx1 ="<<pVx1->get_id()<<" orig 1="<<
//				pVx1->get_orig()<<" Vx2 ="<<pVx2->get_id()<<" orig 2="<<pVx2->get_orig()<<endl;
			}
		} // unlabeled edge count
//	else {nv++;}
	}
		if (pEg->get_dirn()==-2) {ne2++;} // boundary 	
		else if (pEg->get_dirn()==-1) {o2++;} // reverse 	
		else if (pEg->get_dirn()==1) {o1++;} // forward 	
		else if (pEg->get_dirn()==2) {o1++;} // loop 	
		else  {
		ne++;
		} // unlabelled
	mmapedit++;
}
	cout<<" Edge Voroni # of unlableled Edges  = "<< j<<endl;
	cout<<" # of already labeled Edges  = "<< ne3<<endl;
	cout<<" # of boundary Edges  = "<< ne2<<endl;
	cout<<" # of forward Edges  = "<< o1<<endl;
	cout<<" # of reverse Edges  = "<< o2<<endl;
	cout<<" # of relabeled edges  = "<< (j-ne)<<endl;

} // end edge voronoi

// calculate the trip time from the the current stop to the end of the line

// parcel list processing
template <typename s>
s& rideTimetoEnd(s& tsrtmap)
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
// parcel list processing
template <typename o>
void parcstopwalkx(mmaplng& mmapEOIdEgId, mmaplngpar& mmaparced,maplngvx& mapvert, 
				   mmaplnged& mmaped,mmaplngpar& mmaparStop,short onoff,bool blnHistoric, o& logfile )
{
	// mmapV1V2EgId - Vertex 1-2 -> Edge Id multi-map 
	// mmapEOIdEgId - Edge Object Id -> Edge Id multi-map 
	// mmaparced - Parcel ID -> Edge Id Multi-Map 
	// mapvert - Vertex id -> Vertex Object Map 
	// mmaped - edge id -> Edge Object multi-map
	// mmaparStop - Stop Assigned -> Parcel Object multi-map
	// onoff - Offs (away from stops)=0, Ons (toward stops)=1
	// blnHistoric = true if this is historic run, false for alternative runs
    parcel* parp;
    edgev* evp;
    vertexp* vx1,*vx2;
	int j=0;
	long k=0,eid=0,id=0,p1=0,p2=0,p3=0,p4=0; // vid - vertex id, eid - edge id,id - edge seq. id
	double c1=0, wkCost=0;
// Parcels are assigned to edges and edges are assigned to stops prior to this routine being invoked, 
// this routine assignes parcels to stops via the relationshp of edges to stops 
// It also computes walking cost based on edge cost to stops. Cost from parcel to edge is neglected.

mmAPEi = mmaparced.begin();
// loop over all parcels
while(mmAPEi != mmaparced.end()) 
 {  
	  // get the next parcel object and its edge 
      eid =mmAPEi->first;
	  parp = &(mmAPEi->second);
      parp->set_on(onoff); 
	  p1++;
	    if (eid>0) 
	    {
				p2++;
	    	if ((parp->get_origon()<=0 && onoff==1) || (parp->get_origoff()<=0 && onoff==0) ) 
			{  // parcel has not been assigned to a stop yet
			  pair<edfmap_Iter, edfmap_Iter> edrange = mmaped.equal_range((eid));
		      for (mmapedi = edrange.first; mmapedi!=edrange.second;mmapedi++)
			  {
				evp = &(mmapedi->second);
				if (evp->get_orig()>0) { //current edge has valid stop assignment
					mapverit = mapvert.find(evp->get_frid()); // get the end vertex
					vx1 = &(mapverit->second);
					mapverit = mapvert.find(evp->get_toid()); // get the end vertex
					vx2 = &(mapverit->second);
					p3++;
					if(parp->get_palong()<0) {
						double d = (rand() / ((double)RAND_MAX));
						parp->set_palong(d);
					}
					k=evp->get_orig();
					if(evp->get_dirn()==1) {
						if (parp->get_on()==0) 
						{
							parp->set_wkoffcost((evp->get_tcost()-evp->get_scost())*evp->get_palong()+evp->get_scost());
							parp->set_origoff(k);
						} else
						{
							parp->set_wkoncost((evp->get_tcost()-evp->get_scost())*evp->get_palong()+evp->get_scost());
							parp->set_origon(k);
						}
					}
					else if(evp->get_dirn()==-1) {
						if (parp->get_on()==0) 
						{
							parp->set_wkoffcost((evp->get_tcost()-evp->get_scost())*(1-evp->get_palong())+evp->get_scost());
						parp->set_origoff(k);
						} else
						{
							parp->set_wkoncost((evp->get_tcost()-evp->get_scost())*(1-evp->get_palong())+evp->get_scost());
							parp->set_origon(k);
						}
					}
					else if(evp->get_dirn()==-2) // this is a boundary edge
					{ // edge is not assigned to a particular stop and but has valid cost
					 // get the vertices of the edge to see if they have valid cost
						mapverit = mapvert.find(evp->get_frid());
						vx1 = &(mapverit->second);

						mapverit = mapvert.find(evp->get_toid());
						vx2 = &(mapverit->second);

						 // find the indifference point and see which side the parcel is located
							c1 = (abs(vx1->get_cost()-vx2->get_cost())+evp->get_cost())/(2*evp->get_cost());
							if (c1 < parp->get_palong()) // this is closer to vx1
							{
								if (parp->get_on()==0) 
								{
									parp->set_origoff(vx1->get_orig());
									parp->set_wkoffcost(vx1->get_cost()+evp->get_cost()*c1);
								}
								else
								{
									parp->set_origon(vx1->get_orig());
									parp->set_wkoncost(vx1->get_cost()+evp->get_cost()*c1);
								}
							}
							else  // closer to vx2
							{
								if (parp->get_on()==0) 
								{
									parp->set_origoff(vx2->get_orig());
									parp->set_wkoffcost(vx2->get_cost()+evp->get_cost()*c1);
								}
								else
								{
									parp->set_origon(vx2->get_orig());
									parp->set_wkoncost(vx2->get_cost()+evp->get_cost()*c1);
								}
							}
							parp->show_parcel(outparcfile);
							mmaparStop.insert(PE_Pair(k,*parp));
							p4++;
					}
					else 
					{
						// no cost info could be found for parcel
							cout<<"No cost info could be found for parcel "<<parp->get_pacid()<<endl<<
						"Parcel Obj. Id "<<parp->get_id()<<" edge id "<<parp->get_eoid()<<endl;
					} // if origin is fixed for edge
//						parp->set_wkoffcost((abs(vx2->get_cost()-vx1->get_cost()))*parp->get_palong()+evp->get_scost());
					wkCost = (abs(evp->get_tcost()-evp->get_scost()))*parp->get_palong()+evp->get_scost();
					if (parp->get_on()==0) {
						parp->set_origoff(vx2->get_orig());
						parp->set_wkoffcost(wkCost);
					}
					else {
						parp->set_origon(vx1->get_orig());
						parp->set_wkoncost(wkCost);
					}
					parp->show_parcel(outparcfile);
					mmaparStop.insert(PE_Pair(k,*parp));
				} // current edge is assigned to a stop
			  } // loop for multiple edges with same id (edge divided by individual links)   
			} // skip if parcel is assigned to an origin already 
		} // skip if parcel is not mapped to an edge 
   mmAPEi++; // increment parcel iterator
 }; // while parcel is present
 //cout << "  " <<p1<< " Parcels read"<<endl;
 //cout << "  " <<p2<< " Parcels with assigned edges."<<endl;
 //cout << "  " <<p3<< " Parcels with assigned stops."<<endl;
 //cout << "  " <<p4<< " Parcels have assigned edges with out stops."<<endl;
 //cout << " Parcels with stops = " <<(p3+p4)<<" Parcels without stops = "<<(p1-p3-p4)<<endl;
 logfile << "  " <<p1<< " Parcels read"<<endl;
 logfile << "  " <<p2<< " Parcels with assigned edges."<<endl;
 logfile << "  " <<p3<< " Parcels with assigned stops."<<endl;
 logfile << "  " <<p4<< " Parcels assigned to boundary edges."<<endl;
 logfile << " Parcels with stops = " <<(p3+p4)<<" Parcels without stops = "<<(p1-p3-p4)<<endl;

} // end parcel cost calculator
template <typename l, typename p, typename v, typename e, typename o>
p& parcstopwalk(l& mmapEOIdEgId, p& mmaparced,v& mapvert, 
				  e& mmaped,p& mmaparStop,short onoff,bool blnHist, o& logfile)
{
// Parcels are assigned to edges and edges are assigned to stops prior to this routine being invoked, 
// this routine assignes parcels to stops via the relationshp of edges to stops 
// It also computes walking cost based on edge cost to stops. Cost from parcel to edge is neglected.
	// mmapV1V2EgId - Vertex 1-2 -> Edge Id multi-map 
	// mmapEOIdEgId - Edge Object Id -> Edge Id multi-map 
	// mmaparced - Parcel ID -> Edge Id Multi-Map 
	// mapvert - Vertex id -> Vertex Object Map 
	// mmaped - edge id -> Edge Object multi-map
	// mmaparStop - Stop Assigned -> Parcel Object multi-map
	// onoff - Offs (away from stops)=0, Ons (toward stops)=1
	// blnHist = true if this is historic run, false for alternative runs
    parcel* parp;
    edgev* evp;
    vertexp* vx1,*vx2;
	int j=0;
	long k=0,eid=0,id=0,p1=0,p2=0,p3=0,p4=0; // vid - vertex id, eid - edge id,id - edge seq. id
	double c1=0, wkCost=0;

 mmAPEi = mmaparced.begin();
// loop over all parcels
 while(mmAPEi != mmaparced.end()) 
 {  
	  // get the next parcel object and its edge 
    eid =mmAPEi->first;
	parp = &(mmAPEi->second);
	parp->set_on(onoff); 
	p1++;
	    if (eid>0) 
	    {
				p2++;
			// do the same search and match with edge object based on the position along the edge
			pair<edfmap_Iter, edfmap_Iter> edrange = mmaped.equal_range(eid);
		  	size_t j1 = distance(edrange.first,edrange.second);
			if (j1>1)
			{
				evp = &findEdgePart(parp,edrange.first,edrange.second,mmaped,*evp,mapvert,vx1);
			} else if (j1==1) {
				mmapedi = edrange.first;
				evp = &(mmapedi->second);
			}
		      for (mmapedi = edrange.first; mmapedi!=edrange.second;mmapedi++)
			  {
				evp = &(mmapedi->second);
				if (evp->get_orig()>0 ) { //current edge has valid stop assignment
					p3++;
					if(parp->get_palong()<=0) {
						double d = (rand() / ((double)RAND_MAX));
						parp->set_palong(d);
					}
					k=evp->get_orig();
					if(evp->get_dirn()==1 ||evp->get_dirn()==0) {
						if (parp->get_on()==0) 
						{
							if (blnHist) {
								parp->set_hwkoffcost((evp->get_tcost()-evp->get_scost())*evp->get_palong()+evp->get_scost());
							} else
							{
								parp->set_wkoffcost((evp->get_tcost()-evp->get_scost())*evp->get_palong()+evp->get_scost());
							}
							parp->set_origoff(k);
						} else
						{
							if (blnHist) {
								parp->set_hwkoncost((evp->get_tcost()-evp->get_scost())*evp->get_palong()+evp->get_scost());
							} else 
							{
								parp->set_wkoncost((evp->get_tcost()-evp->get_scost())*evp->get_palong()+evp->get_scost());
							}

							parp->set_origon(k);
						}
					}
					else if(evp->get_dirn()==-1) {
						if (parp->get_on()==0) {	
							if (blnHist) {
	  							parp->set_hwkoffcost((evp->get_tcost()-evp->get_scost())*(1-evp->get_palong())+evp->get_scost());
							} else
							{
	  							parp->set_wkoffcost((evp->get_tcost()-evp->get_scost())*(1-evp->get_palong())+evp->get_scost());
							}
							parp->set_origoff(k);
						} else
						{
							if (blnHist) {
		  						parp->set_hwkoncost((evp->get_tcost()-evp->get_scost())*(evp->get_palong())+evp->get_scost());
							} else
							{
		  						parp->set_wkoncost((evp->get_tcost()-evp->get_scost())*(evp->get_palong())+evp->get_scost());
							}
							parp->set_origon(k);
						}

					}
					else if (evp->get_dirn()==-2)   
					{
						mapverit = mapvert.find(evp->get_frid());
						vx1 = &(mapverit->second);

						mapverit = mapvert.find(evp->get_toid());
						vx2 = &(mapverit->second);

						assignParcStop( vx1, vx2, evp, parp,blnHist);

						parp->show_parcel(outparcfile);
					}
					else {
						// no cost info could be found for parcel
							cout<<"No cost info could be found for parcel "<<parp->get_pacid()<<endl<<
						"Parcel Obj. Id "<<parp->get_id()<<" edge id "<<parp->get_eoid()<<endl;
					} // if origin is fixed for edge
				} // current edge is assigned to a stop
				else if(evp->get_orig()==-2) // this is a boundary edge
					{ // edge is not assigned to a particular stop and but has valid cost
					 // get the vertices of the edge to see if they have valid cost
						mapverit = mapvert.find(evp->get_frid());
						vx1 = &(mapverit->second);

						mapverit = mapvert.find(evp->get_toid());
						vx2 = &(mapverit->second);
						assignParcStop( vx1, vx2, evp, parp,blnHist);

							parp->show_parcel(outparcfile);
							p4++;
					}
					if (parp->get_on()==0) {
						if (parp->get_origoff()>0 && parp->get_frid()!=-999) {
							mmaparStop.insert(PE_Pair(parp->get_origoff(),*parp));
							parp->set_frid(-999);
						}
					} else if (parp->get_on()==1)
					{
						if (parp->get_origon()>0 && parp->get_toid()!=-999) {
							mmaparStop.insert(PE_Pair(parp->get_origon(),*parp));
							parp->set_toid(-999);
						}
					}

			  } // loop for multiple edges with same id (edge divided by individual links)   
			//} // skip if parcel is assigned to an origin already 
		} // skip if parcel is not mapped to an edge 
   mmAPEi++; // increment parcel iterator
 }; // while parcel is present
/* cout << "  " <<p1<< " Parcels read"<<endl;
 cout << "  " <<p2<< " Parcels with assigned edges."<<endl;
 cout << "  " <<p3<< " Parcels with assigned stops."<<endl;
 cout << "  " <<p4<< " Parcels assigned to boundary edges."<<endl;
 cout << " Parcels with stops = " <<(p3+p4)<<" Parcels without stops = "<<(p1-p3-p4)<<endl;
*/
 logfile <<endl<< "  " <<p1<< " Parcels read"<<endl;
 logfile << "  " <<p2<< " Parcels with assigned edges."<<endl;
 logfile << "  " <<p3<< " Parcels with assigned stops."<<endl;
 logfile << "  " <<p4<< " Parcels assigned to boundary edges."<<endl;
 logfile << " Parcels with stops = " <<(p3+p4)<<" Parcels without stops = "<<(p1-p3-p4)<<endl<<endl;

} // end parcel cost calculator

template <typename l, typename p, typename v, typename e, typename o>
p& parcstopwalknew(l& mmapEOIdEgId, p& mmaparced, v& mapvert, e& mmaped,
				  p& mmaparStop,short onoff,bool blnHist, o& logfile )
{
// Parcels are assigned to edges and edges are assigned to stops prior to this routine being invoked, 
// this routine assignes parcels to stops via the relationshp of edges to stops 
// It also computes walking cost based on edge cost to stops. Cost from parcel to edge is neglected.
	// mmapV1V2EgId - Vertex 1-2 -> Edge Id multi-map 
	// mmapEOIdEgId - Edge Object Id -> Edge Id multi-map 
	// mmaparced - Parcel ID -> Edge Id Multi-Map 
	// mapvert - Vertex id -> Vertex Object Map 
	// mmaped - edge id -> Edge Object multi-map
	// mmaparStop - Stop Assigned -> Parcel Object multi-map
	// onoff - Offs (away from stops)=0, Ons (toward stops)=1
	// blnHist = true if this is historic run, false for alternative runs
    parcel* parp;
    edgev* pEg;
    vertexp* pVx1,*pVx2;
	int j=0;
	long k=0,eid=0,id=0,p1=0,p2=0,p3=0,p4=0; // vid - vertex id, eid - edge id,id - edge seq. id
	double c1=0, wkCost=0,posx=0,d=0;

 mmAPEi = mmaparced.begin();
// loop over all parcels
 while(mmAPEi != mmaparced.end()) 
 {  
	  // get the next parcel object and its edge 
    eid =mmAPEi->first;
	parp = &(mmAPEi->second);
	parp->set_on(onoff); 
	p1++;
	    if (eid>0) 
	    {
				p2++;
			// do the same search and match with edge object based on the position along the edge
			pair<edfmap_Iter, edfmap_Iter> edrange = mmaped.equal_range(eid);
		  	size_t j1 = distance(edrange.first,edrange.second);
			if (j1==0) 
			{// skip edge
			} 
			else { // parcel is attached to an edge in the list 	
				if (j1>1)
				{
					pEg = &findEdgePart(parp,edrange.first,edrange.second,mmaped,*pEg,mapvert,pVx1);
				} else if (j1==1) {
					mmapedi = edrange.first;
					pEg = &(mmapedi->second);
				}
				pEg = &(mmapedi->second);
				mapverit = mapvert.find(pEg->get_frid());
				pVx1 = &(mapverit->second);
				mapverit = mapvert.find(pEg->get_toid());
				pVx2 = &(mapverit->second);
				if (pEg->get_orig()>0 ) { //current edge has valid stop assignment
					p3++;
					k=pEg->get_orig();
				} else {
					k=-1;
				}
				if(parp->get_palong()<0) {
					d = (rand() / ((double)RAND_MAX));
					parp->set_palong(d);
				} else {
					d = parp->get_palong();
				}
				if (k>0) {
						wkCost=min((pEg->get_cost())*(1-parp->get_palong())+pEg->get_tcost(),
							(pEg->get_cost())*(parp->get_palong())+pEg->get_scost());
					if(pEg->get_dirn()==1 || pEg->get_dirn()==0) {// forward 
						if (parp->get_on()==0) 
						{
							if (blnHist) {
								parp->set_hwkoffcost(wkCost);
							} else
							{
								parp->set_wkoffcost(wkCost);
							}
							parp->set_origoff(k);
						} else
						{
							if (blnHist) {
								parp->set_hwkoncost(wkCost);
							} else 
							{
								parp->set_wkoncost(wkCost);
							}

							parp->set_origon(k);
						}
					}
					else if(pEg->get_dirn()==-1) {
						if (parp->get_on()==0) {	
							if (blnHist) {
	  							parp->set_hwkoffcost(wkCost);
							} else {
	  							parp->set_wkoffcost(wkCost);
							}
							parp->set_origoff(k);
						} else
						{
							if (blnHist) {
		  						parp->set_hwkoncost(wkCost);
							} else {
		  						parp->set_wkoncost(wkCost);
							}
							parp->set_origon(k);
						}

					}
					else if (pEg->get_dirn()==-2)   
					{
						mapverit = mapvert.find(pEg->get_frid());
						pVx1 = &(mapverit->second);

						mapverit = mapvert.find(pEg->get_toid());
						pVx2 = &(mapverit->second);

						assignParcStop( pVx1, pVx2, pEg, parp,blnHist);

						parp->show_parcel(outparcfile);
					}
					else {
						// no cost info could be found for parcel
							cout<<"No edge direction info could be found for parcel "<<parp->get_pacid()<<endl<<
						"Parcel Obj. Id "<<parp->get_id()<<" edge id "<<parp->get_eoid()<<endl;
							logfile<<"No edge direction info could be found for parcel "<<parp->get_pacid()<<endl<<
						"Parcel Obj. Id "<<parp->get_id()<<" edge id "<<parp->get_eoid()<<endl;
					} // if origin is fixed for edge
				} // current edge is assigned to a stop
				else if(pEg->get_orig()==-2) // this is a boundary edge
				{ // edge is not assigned to a particular stop and but has valid cost
					 // get the vertices of the edge to see if they have valid cost
						mapverit = mapvert.find(pEg->get_frid());
						pVx1 = &(mapverit->second);

						mapverit = mapvert.find(pEg->get_toid());
						pVx2 = &(mapverit->second);
						assignParcStop( pVx1, pVx2, pEg, parp,blnHist);

							parp->show_parcel(outparcfile);
							p4++;
				}
					if (parp->get_on()==0) {
						if (parp->get_origoff()>0 && parp->get_frid()!=-999) {
							mmaparStop.insert(PE_Pair(parp->get_origoff(),*parp));
							parp->set_frid(-999);
						}
					} else if (parp->get_on()==1)
					{
						if (parp->get_origon()>0 && parp->get_toid()!=-999) {
							if (parp->get_origon()!=parp->get_origoff()) {
								d = parp->get_wkoffcost() ;
								posx = parp->get_wkoncost() ;
							}
							mmaparStop.insert(PE_Pair(parp->get_origon(),*parp));
							parp->set_toid(-999);
						}
					}
			  } // if edge is found in the edge set j>0   
		} // skip if parcel does not have an edge id in the object data (not mapped to an edge) 
   mmAPEi++; // increment parcel iterator
 }; // while parcel is present
 cout <<endl<< " End of Parcel Assignment for " <<mapMsg.find(onoff)->second<< endl<<endl;
 /* cout << "  " <<p1<< " Parcels read"<<endl;
 cout << "  " <<p2<< " Parcels with assigned edges."<<endl;
 cout << "  " <<p3<< " Parcels with assigned stops."<<endl;
 cout << "  " <<p4<< " Parcels assigned to boundary edges."<<endl;
 cout << " Parcels with stops = " <<(p3+p4)<<" Parcels without stops = "<<(p1-p3-p4)<<endl<<endl;*/
 logfile <<endl<< " End of Parcel Assignment for " <<mapMsg.find(onoff)->second<<endl<< "  " <<p1<< " Parcels read!"<<endl;
 logfile << "  " <<p2<< " Parcels with assigned edges."<<endl<< "  " <<p3<< " Parcels with assigned stops."<<endl;
 logfile << "  " <<p4<< " Parcels assigned to boundary edges."<<endl<< " Parcels with stops = " 
	 <<(p3+p4)<<endl<<" Parcels without stops = "<<(p1-p3-p4)<<endl<<endl;

return mmaparStop;
} // end parcel cost calculator


deque<ParcelOffAccum> pagwalkstopoff (maplngstop& tsidmap, mmaplngpar& mmaparStop,lucodes* luci,
				  globcost* gc1,short onoff,ParcelOffAccum& pa,short blnHist, ofstream& logfile)
{
    long ip;
	tstop* pstop;
    parcel* parp;
	double propon=0,propoff=0;
deque<parcel> vpar;
deque<ParcelOffAccum> vpa;
deque<ParcelOffAccum>::iterator vpait;
 tsidmit = tsidmap.begin();

// aggregate parcel walk time impact by transit stop
 while (tsidmit !=tsidmap.end()) {
	 pstop = &(tsidmit->second);
	 ip = pstop->get_id();
// get the range of parcel objects of the current stop id
 		pair<mmapP_AIter, mmapP_AIter> pastrange = mmaparStop.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mmAPEi=pastrange.first; mmAPEi!=pastrange.second;mmAPEi++)
	    {
//  		  par1 = mmAPEi->second;
    	   parp = &(mmAPEi->second);

   // compute the on/off strengths per parcel
			if (onoff == 0 && parp->get_hwkoffcost()>=0) {
					propoff = exp(-gc1->get_propensity()/300*(parp->get_hwkoffcost() * 60 * gc1->get_walkspd()));
			} else	if (onoff == 1 && parp->get_hwkoncost()>=0) {
					propon = exp(-gc1->get_propensity()/300*(parp->get_hwkoncost())* 60 * gc1->get_walkspd());
//					propon = -gc1->get_propensity()*parp->get_hwkoncost()*gc1->get_walkspd()*60/300;
			} 
/*			else if (onoff == 0 && parp->get_wkoffcost()>=0) {
					propoff = log(gc1->get_propensity())/300*(parp->get_wkoffcost());//* 60 * gc1->get_walkcost());
//					propoff = -gc1->get_propensity()*parp->get_wkoffcost()*gc1->get_walkspd()*60/300;
			} else	if (onoff == 1 && parp->get_wkoncost()>=0) {
					propon = log(gc1->get_propensity())/300*(parp->get_wkoncost());//* 60 * gc1->get_walkcost());
//					propon = -gc1->get_propensity()*parp->get_wkoncost()*gc1->get_walkspd()*60/300;
			}
*/
			string strPtype = parp->get_ptype();
		       lucmit = lumap.find(strPtype);
			   if (lucmit!=lumap.end()) {
			     luci = &(lucmit->second);
				 parp->set_lucd(luci->get_LUC());
	 		     if (luci->get_KeyProp()=="Gross_Area" || luci->get_KeyProp()=="GROSS_AREA"
					 || luci->get_KeyProp()=="GROSSAREA" || luci->get_KeyProp()=="AREA" || luci->get_KeyProp()=="Area") {
						if (onoff) {
							propon=(luci->get_OnCoeff()*parp->get_grarea()*parp->get_cfact()* propon);
							parp->set_pval(propon);
						}
						else 
						{
							propoff=(luci->get_OffCoeff()*parp->get_grarea()*parp->get_cfact()* propoff);
							parp->set_aval(propoff);
						}
				 }
				 else if (luci->get_KeyProp()=="Living_are" || luci->get_KeyProp()=="LIVINGAREA" 
					 || luci->get_KeyProp()=="LIVING_ARE") {
						if (onoff) {
							propon = (luci->get_OnCoeff() * parp->get_lvarea()*parp->get_cfact()* propon);
							parp->set_pval(propon);
						}
						else {
							propoff = (luci->get_OffCoeff() * parp->get_lvarea()*parp->get_cfact()* propoff);
							parp->set_aval(propoff);
						}
				 }
				vpar.push_back((*parp));
			   }
			   else
			   { // skip since land use code data is not found
				cout<<"Land use code data not found for Type "<<strPtype<< " for parcel "<<parp->get_pacid();
				logfile<<"Land use code data not found for Type "<<strPtype<< " for parcel "<<parp->get_pacid();
			   }
		}
		
  vpa.push_back(for_each(vpar.begin(),vpar.end(), ParcelOffAccum()));
  vpar.clear();
  tsidmit++;
 }
return vpa;
}


deque<ParcelOnAccum> pagwalkstopon (maplngstop& tsidmap, mmaplngpar& mmaparStop,lucodes* luci,
				  globcost* gc1,short onoff,ParcelOnAccum& pa,short blnHist, ofstream& logfile)
{
    long ip;
	tstop* pstop;
    parcel* parp;
	double propon=0,propoff=0;
deque<parcel> vpar;
deque<ParcelOnAccum> vpa;
deque<ParcelOnAccum>::iterator vpait;
 tsidmit = tsidmap.begin();

// aggregate parcel walk time impact by transit stop

///*
 while (tsidmit !=tsidmap.end()) {
	 pstop = &(tsidmit->second);
	 ip = pstop->get_id();
// get the range of parcel objects of the current stop id
 		pair<mmapP_AIter, mmapP_AIter> pastrange = mmaparStop.equal_range(ip);
	   	size_t j = distance(pastrange.first,pastrange.second);
		for (mmAPEi=pastrange.first; mmAPEi!=pastrange.second;mmAPEi++)
	    {
//  		  par1 = mmAPEi->second;
    	   parp = &(mmAPEi->second);

   // compute the on/off strengths per parcel
			if (onoff == 0 && parp->get_hwkoffcost()>=0) {
//					propoff = -gc1->get_propensity()*parp->get_hwkoffcost()*gc1->get_walkspd()*60/300;
					propoff = exp(-gc1->get_propensity()*parp->get_hwkoffcost()*gc1->get_walkspd()*60/300);
			} else	if (onoff == 1 && parp->get_hwkoncost()>=0) {
					propon = exp(-gc1->get_propensity()*parp->get_hwkoncost()*gc1->get_walkspd()*60/300);
			}
               string strPtype = parp->get_ptype();
		       lucmit = lumap.find(strPtype);
			   if (lucmit!=lumap.end()) {
			     luci = &(lucmit->second);
				 parp->set_lucd(luci->get_LUC());
	 		     if (luci->get_KeyProp()=="Gross_Area" || luci->get_KeyProp()=="GROSS_AREA"
					 || luci->get_KeyProp()=="GROSSAREA") {
						if (onoff) {
							propon=(luci->get_OnCoeff()*parp->get_grarea()*parp->get_cfact()* propon);
							parp->set_pval(propon);
						}
						else 
						{
							propoff=(luci->get_OffCoeff()*parp->get_grarea()*parp->get_cfact()*propoff);
							parp->set_aval(propoff);
						}
				 }
				 else if (luci->get_KeyProp()=="Living_are" || luci->get_KeyProp()=="LIVINGAREA" 
					|| luci->get_KeyProp()=="LIVING_ARE") {
						if (onoff) {
							propon = (luci->get_OnCoeff() * parp->get_lvarea()*parp->get_cfact()* propon);
							parp->set_pval(propon);
						}
						else {
							propoff = (luci->get_OffCoeff() * parp->get_lvarea()*parp->get_cfact()* propoff);
							parp->set_aval(propoff);
						}
				 }
				vpar.push_back((*parp));
			   }
			   else
			   { // skip since land use code data is not found
				cout<<"Land use code data not found for Type "<<strPtype<< " for parcel "<<parp->get_pacid();
				logfile<<"Land use code data not found for Type "<<strPtype<< " for parcel "<<parp->get_pacid();
			   }
		}
		
  vpa.push_back(for_each(vpar.begin(),vpar.end(), ParcelOnAccum()));
  vpar.clear();
  tsidmit++;
 }
return vpa;
}


void fileName (char  fname[],char  ext[],char  outfname[], int onoff,int alt)
{
	char* dot = strrchr(fname,'.');
	char ch1[15];
    int	ix = dot-fname;
	if (ix>0) {
		strncpy(outfname,fname, ix );
        outfname[ix]='\0';
	}
	else
    { ix = sizeof(fname);
		outfname[0]='\0';
	strcat(outfname,fname);}

	if (onoff>=0) {
		strcat(outfname,"_");
		ix = sizeof(fname);
		dot = ch1;
		strcat(outfname,_itoa(onoff,dot,10));
	}
	if (alt>0) {
		strcat(outfname,"_");
		ix = sizeof(fname);
		dot = ch1;
		strcat(outfname,_itoa(alt,dot,10));
	}
	{strcat(outfname,ext);}
	ext='\0';
}

void fileName (char  fname[],string  ext,char  outfname[], int onoff,int alt)
{
	char* dot = strrchr(fname,'.');
	char ch1[15];
	char ch2[30];
    int	ix = dot-fname;
	if (ix>0) {
		strncpy(outfname,fname, ix );
        outfname[ix]='\0';
	}
	else
    { ix = sizeof(fname);
		outfname[0]='\0';
	strcat(outfname,fname);}

	if (onoff>=0) {
		strcat(outfname,"_");
		ix = sizeof(fname);
		dot = ch1;
		strcat(outfname,_itoa(onoff,dot,10));
	}
	if (alt>0) {
		strcat(outfname,"_");
		ix = sizeof(fname);
		dot = ch1;
		strcat(outfname,_itoa(alt,dot,10));
	}
	strcpy(ch2, ext.c_str());
	{strcat(outfname,ch2);}
	ext='\0';
}




maplngstop& stopspacing_analysis (mmaplng& mmapV1V2EOId, mmaplng& mmapV1V2EgId,
			mmaplngpar& mmaparced,maplngvx& mapvert, maplnged& maped,mmaplngpar& mmaparStop,
			mmapdblng& mmapSver,mmaplng& mmapVxStop,mmaplng& mmapStopVx,mmapdblng& mmapVx0,
			mmapdblng& mmapVx1,mmaplngdbl& mmapVx0R,mmaplngdbl& mmapVx1R, bool blnHist, bool blnEuclid,
			string strEdgeBaseName, string strVertexBaseName, string strParcelBaseName,string strStopBaseName,
			globcost& gc1, pdhdway phway, lucodes luci,maplngstop& stops, int q, ofstream& logfile,
			inputfilelist& inpList, stopKey& kStop,sqlite3* db, string sADB, int xi) 
{ 
	ofstream outfile; // output file name
	ifstream infile; // input file name
	//char outfilename[ MaxStrLen +1]=""; // data output file 
	//char infilename[ MaxStrLen +1]=""; // data input file 

// object declarations
  
	vertexp vx; //vertex object
	vertexp* pVx1=&vx;

	edgev ev; // edge obejct
	edgev* pEg=&ev;
	maplngvx mapvert0;
	maplngvx mapvert1;
	maplngvx mapvert2;
	mmaplngvx mmapvert0;
	mmaplngvx mmapvert1;
	maplnged maped0;
	maplnged maped1;
	mmaplnged mmaped0;
	mmaplnged mmaped1;
//initialize the parcel object
	parcel par1;
	parcel* parp=&par1;
	long i1=0,k1=0, ip;
	int  alt=0, r=1; // r - indicates if the run is historic or alternative  
	double x1;
	short D=gc1.get_dpdimension();
	//char ext2[20]; // extension for creating unique output filenames between the various dimensions
	short onoff = 0;
	string str1 = "",strSQL="", tblDef, keyFld;
	bool qNet=true, blnCreate=false, blnInsert;
//initialize the stop variables needed for the program
	tstop stop0;
	tstop* pstop=&stop0;	
	sqlite3_stmt *stmti = NULL;
	maplngstop::iterator sit;
	int intSRID=0;
	EdgeAccum edgeAcc;
	deque<EdgeAccum> edgeSumAcc;
	string tblEdge = "" , tblVertex="", tblParcel="", tblTrip="" , viewName="" , strInFname="",strOutFname="";
	if (q>0) {
		r=1; // if this is not a historic run a dynamic programming run 
	}  else {r=0;}

	cout <<endl<< " Starting stop spacing analysis... run "<<q<<endl;
	logfile <<endl<< " Starting stop spacing analysis... run "<<q<<endl;
	if (!blnEuclid) {
// read the original edge data from the edge file
		// fileName(edgeinfname,".bin",infilename);
		strInFname = strEdgeBaseName + ".bin";

		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		maped0.clear();
		readBinEdgeObjectData(maped0,ev,infile);
//		readBinIdMap(maped0,ev,ip,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

// transfer the vertex network to the two on and off multi-maps

		//fileName(vertinfname,".bin",infilename);
		strInFname = strVertexBaseName + ".bin";
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		mapvert0.clear();
		readBinVertexObjectData(mapvert0,vx,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

// read the origin vertices data from the origin file
		//fileName(tstopiname,"_Orig.bin",infilename);
		strInFname = strStopBaseName + "_Orig.bin";
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		mmapSver2.clear(); 
		readBindblong(mmapSver2,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

		onoff=0;
		//fileName(vertinfname,"Vx.bin",infilename,onoff);
		strInFname = strVertexBaseName + "_" + to_string<long>(onoff) +"Vx.bin";

		mmapVx0.clear();
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		readBindblong(mmapVx0,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

		//fileName(vertinfname,"VxR.bin",infilename,onoff);
		strInFname = strVertexBaseName + "_" + to_string<long>(onoff) +"VxR.bin";

		mmapVx0R.clear();
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		mmapVx0R = readBindualvar(mmapVx0R,i1,x1,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear(); 

		//fileName(tstopiname,"_VxStop.bin",infilename);
		strInFname = strStopBaseName + "_VxStop.bin";
		mmapVxStop.clear();
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		mmapVxStop = readBindualvar(mmapVxStop,i1,k1,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear(); 

		//fileName(tstopiname,"_StopVx.bin",infilename);
		strInFname = strStopBaseName + "_StopVx.bin";

		mmapStopVx.clear();
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		mmapStopVx = readBindualvar(mmapStopVx,i1,k1,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear(); 
	
// search and remove the start vertices that are not in the current stop set
  		if (qNet) {
			 maped0 =  fixVxEdgeCoStop(maped0,mapEidOid,ev,stops,stop0,mapvert0,*pVx1,
				mmapVxStop,mmapStopVx,mmapStopEdge,mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,
				mmapVx0R,mmapVx1R,mapvertor0,mapvertor1,outedgefile);
		} else {
			maped0 =  funCoStopVxEdge(maped0,mapEidOid,ev,stops,stop0,mapvert0,*pVx1,
				mmapVxStop,mmapStopVx,mmapStopEdge,mmapEdgeStop,mmapSver,mmapVx0,mmapVx1,
				mmapVx0R,mmapVx1R,mapvertor0,mapvertor1,outedgefile);
		}

//		mmapVxStop.clear();
//		mmapVxStop = reverseMapKeys(mmapStopVx,mmapVxStop,i1,k1);

//update/set the start vertex costs of the set of stops before starting the voronoi routine
//		mmapVx0 = dpSynchDblVertexMap(mmapStopVx,stops,mmapVx0R,mmapVx0,blnHist);
// vertex, vertex-origin,stop vertex
		updateVxOrigins0(mmapStopVx,mmapVx0R,mmapVx0,mapvertor0,mapvert0,*pVx1,i1);
//		updateVxOrigins2(stops,stop0,mmapStopVx,mapvertor0,mapvert0,*pVx1,i1,onoff,maped0,ev);

// call the vertex voronoi routine
 
// mapvert0 = vxv (mmapVx0,mapvert0,mmaped0,mmapV1V2EgId,mmapVxStop,onoff);

// first remap the edge map to the edge id as key 
		mmaped0 = remapObj2Id (maped0,mmaped0,ev,i1);
		mapvert2.clear();
		mapvert2 = vxvtarx (mmapV1V2EgId,mmapVx0,mmapVxStop,mapvert0,mapvert2,mmaped0,onoff,gc1);

// aggregate the vertex results
// first remap the key of the vertices by stop
		mmapvert0.clear();
		mmapvert0 = remapVertId2Orig(mapvert0,mmapvert0,*pVx1,i1);	

// Subtract out the Ride time from the begining for offs & the Ride time to the end for ons from the vertex cost

		mmapvert0 = reduceVxToTime2WalkTime(mmapvert0,*pVx1,stops,stop0,onoff,logfile);	

		mapvert0.clear();
		mapvert0 = remapObj2Id(mmapvert0,mapvert0,*pVx1,i1);	

	// write out the vertex results
		//fileName(vertinfname,".vxv",outfilename,onoff,r);
		strOutFname = strVertexBaseName + "_" + to_string<short>(onoff)+ "_" + to_string<long>(r) + ".vxv";
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

		if (q>0) {
			outfile.open(strOutFname.c_str(),ios::out |ios::app);
		} else {
			outfile.open(strOutFname.c_str(),ios::out |ios::trunc);
		}
		VertAccum vertAcc;
		deque<VertAccum> vertSumAcc;
		vertSumAcc = aggbystop (stops,mmapvert0,*pstop,*pVx1,vertAcc);
//// create table of vertex voronoi result
//		str1 =	"( VERTID INTEGER  PRIMARY KEY, EDGEID  INTEGER default -1, IDP INTEGER default -1, "
//			" COST DOUBLE default 0.0, IdHigh INTEGER default -1, IdLow Integer default -1, "
//			" Stop INTEGER default -1, Geometry POINT ); ";
//		intSRID = InpFiles.get_srid();
//		tblVertex = InpFiles.get_tblvertex() + "D" + to_string<short> (D) + "i" + to_string<int> (q); 
//		blnCreate = createSpaTbl(str1,tblVertex,db,intSRID,outfile);
//		if (blnCreate) {
//		// write the vedrtex output to the table
//			keyFld = "VERTID";
//			string strSQLInserTbl = "Insert into " + tblVertex + " (\"VertId\",\"IDP\"," + 
//				"\"COST\",\"IdHigh\",\"IdLow\",\"Stop\",\"Geometry\" ) " + 
//				" Values ( ?, ?, ?, ?, ?, ?, ?, ? );" ;
//		intSRID = InpFiles.get_srid();
//		// call insert table routine to populate the vertex table
//			blnInsert = inSpaTbl(strSQLInserTbl,db,*pVx1,mmapvert0,tblVertex, keyFld, intSRID ,logfile);
//		}

		if ((q%xi)==0 || blnHist || q==1) {
			datetimeStamp(outfile);
			str1 = "Summary of Vertex attributes by id iteration " + to_string(q) + 
			"\nId\tStop\tcost\tTotCost\tParCnt\tindex\tLowLink";
			outfile << str1<<endl;
			writeObjData(vertSumAcc,vertAcc,str1,outfile);
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();

			str1 = "D" + to_string<short>(D);
			str1.append("vorvx.txt");
	//		ext2[0]='\0';
	//		strcpy(ext2, str1.c_str());
			//fileName(vertinfname,ext2,outfilename,onoff,q);
			strOutFname = strVertexBaseName + "_" + to_string<short>(onoff)+ "_" + to_string<long>(q) + str1;
			outfile.open(strOutFname.c_str(), ios::out | ios::trunc);
			vx.serializetexthdr(outfile);
			writeTextObjectData(mmapvert0,vx,i1,outfile,"");

		}
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

  // save the vertex voronoi result 
		//fileName(vertinfname,"vorvx.bin",outfilename,onoff,r);
		strOutFname = strVertexBaseName + "_" + to_string<short>(onoff)+ "_" + to_string<long>(r) + "vorvx.bin";

		outfile.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
		writeBinObjectData(mmapvert0,vx,ip,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();  


// write the edge object data
		//fileName(edgeinfname,"eg.bin",outfilename,onoff,r);
		strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff)+ "_" + to_string<long>(r) + "eg.bin";
		outfile.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
		writeBinObjectData(maped0,ev,ip,outfile);

		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

// do the ons vertex voronoi data i.e towards the stop using the ride time from the stop to the end	
		onoff = 1;
// read the original edge data from the edge file
		//fileName(edgeinfname,".bin",infilename);
		strInFname = strEdgeBaseName + ".bin";
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		maped1.clear();
		readBinEdgeObjectData(maped1,ev,infile);
		//		readBinIdMap(maped1,ev,ip,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

// read vertex data for the second time

		//fileName(vertinfname,".bin",infilename);
		strInFname = strVertexBaseName + ".bin";
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		mapvert1.clear();
		readBinVertexObjectData(mapvert1,vx,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

// search and remove the start vertices that are not in the current stop set
// first read origin vertx data 
		//fileName(vertinfname,"Vx.bin",infilename,onoff);
		strInFname = strVertexBaseName + "_" + to_string<short>(onoff) + "Vx.bin";
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();
		mmapVx1.clear();
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		readBindblong(mmapVx1,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear(); 

		//fileName(vertinfname,"VxR.bin",infilename,onoff);
		strInFname = strVertexBaseName + "_" + to_string<short>(onoff) + "VxR.bin";
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

		mmapVx1R.clear();
		infile.open(strInFname.c_str(),  ios::in | ios::binary);
		mmapVx1R = readBindualvar(mmapVx1R,i1,x1,infile);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear(); 

//  mmaped1 =  upCoStopVxEdge(mmaped1,mapEidOid,ev,stops,stop0,mapvert0,pVx1,mapvertor1,outedgefile);

 //update set the start vertex costs of the set of stops before starting the voronoi routine
/*  maped1 = funCoStopVxEdge(maped1,mapEidOid,ev,stops,stop0,mapvert1,*pVx1,
	  mmapVxStop,mmapStopVx,mmapStopEdge,mmapEdgeStop,mmapSver2,mmapVx0,mmapVx1,
	  mmapVx0R,mmapVx1R,mapvertor0,mapvertor1,outedgefile);
	  */
//		mmapVx1 = dpSynchDblVertexMap(mmapStopVx,stops,mmapVx1R,mmapVx1,blnHist);
//update set the start vertex costs of the set of stops before starting the voronoi routine
		updateVxOrigins0(mmapStopVx,mmapVx1R,mmapVx1,mapvertor1,mapvert1,*pVx1,i1);
		//updateVxOrigins(mmapVx1,mmapVxStop,mapvertor1,mapvert1,*pVx1,i1);
//		updateVxOrigins2(stops,stop0,mmapStopVx,mapvertor1,mapvert1,*pVx1,i1,onoff,maped1,ev);
  
  // call the voronoi vertext routine
		// remap edge id as the key 
		mmaped1 = remapObj2Id (maped1,mmaped1,ev,i1);
 
		mapvert2.clear();
		mapvert2 = vxvtarx (mmapV1V2EgId,mmapVx1,mmapVxStop,mapvert1,mapvert2,mmaped1,onoff,gc1);
 
// aggregate the vertex results
// first remap the key of the vertices by stop
		mmapvert1.clear();
		mmapvert1 = remapVertId2Orig(mapvert1,mmapvert1,*pVx1,i1);	

		mmapvert1 = reduceVxToTime2WalkTime(mmapvert1,*pVx1,stops,stop0,onoff,logfile);	
		mapvert1.clear();
		mapvert1 = remapObj2Id(mmapvert1,mapvert1,*pVx1,i1);	

		//fileName(vertinfname,"vorvx.bin",outfilename,onoff,r);
		strOutFname = strVertexBaseName + "_" + to_string<long>(onoff) + "_" + to_string<long>(r) +"vorvx.bin";
		outfile.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
		writeBinObjectData(mapvert1,vx,ip,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

		if ((q%xi)==0 || blnHist  || q==1) {

			// write out the vertex results
			str1 = to_string<short>(onoff) + "_" + to_string<short>(r) + "D" + to_string<short>(D) + ".vxv";
			//ext2[0]='\0';
			//strcpy(ext2, str1.c_str());
			//fileName(vertinfname,ext2,outfilename,onoff,r);
			strOutFname = strVertexBaseName + "_" +  str1;
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();
			if (q>0) {
				outfile.open(strOutFname.c_str(),ios::out |ios::app);
			} else {
				outfile.open(strOutFname.c_str(),ios::out |ios::trunc);
			}
			vertSumAcc.clear();
			vertSumAcc = aggbystop (stops,mmapvert1,*pstop,*pVx1,vertAcc);
			str1 = "Summary of Vertex attributes by id " + to_string(q) + 
			" \nId\tStop\tcost\tTotCost\tParCnt\tindex\tLowLink";
			datetimeStamp(outfile);
			outfile << str1<<endl;
			writeObjData(vertSumAcc,vertAcc,str1,outfile);
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();

			str1 = "D" + to_string<short>(D);
			str1.append("vorvx.txt");
			//ext2[0]='\0';
			//strcpy(ext2, str1.c_str());
			//fileName(vertinfname,ext2,outfilename,onoff,q);
			strOutFname = strVertexBaseName + "_" + to_string<long>(onoff) + "_" + to_string<long>(q) +str1;
			outfile.open(strOutFname.c_str(), ios::out | ios::trunc);
			vx.serializetexthdr(outfile);
			writeTextObjectData(mapvert1,vx,i1,outfile,"");
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();

// create table of vertex voronoi result for on and offs 
			str1 =	"( VERTID INTEGER PRIMARY KEY , IDP INTEGER default -1, StopOn INTEGER default -1,"
				" CostOn DOUBLE default 0.0, IdHighOn INTEGER default -1, IdLowOn Integer default -1, "
				" StopOff INTEGER default -1,CostOff DOUBLE default 0.0, IdHighOff INTEGER default -1," 
				"  IdLowOff Integer default -1); ";
			intSRID = inpList.get_srid();
			if (blnHist) {
				tblVertex = inpList.get_tblvertex() + "D" + to_string<short> (D) + "hi" + to_string<int> (q)+ "t" + to_string<long> (kStop.tripId()) + kStop.tripPeriod();  
			} else {
				tblVertex = inpList.get_tblvertex() + "D" + to_string<short> (D) + "ai" + to_string<int> (q)+ "t" + to_string<long> (kStop.tripId()) + kStop.tripPeriod();  
			}
			replace(tblVertex.begin(),tblVertex.end(),' ','_');
			ReplaceAll2(tblVertex,"__","_");
			blnCreate = createSpaTbl(str1,tblVertex,db,intSRID,outfile);
			if (blnCreate) {
			// write the vedrtex output to the table
				keyFld = "VERTID";
				string strSQLInserTbl = "Insert into " + tblVertex + " (\"VertId\",\"IDP\","  
					"\"StopOn\",\"CostOn\",\"IdHighOn\",\"IdLowOn\"," 
					"\"StopOff\",\"CostOff\",\"IdHighOff\",\"IdLowOff\" ) "  
					" Values ( ?, ?, ?, ?, ?, ?, ?, ?,?,? );" ;

			// call insert table routine to populate the vertex table
				blnInsert = inSpaTblVx(strSQLInserTbl,db,*pVx1,mapvert1,mapvert0,tblVertex, keyFld , intSRID,logfile);
				//if (blnInsert) {
				//	qrySQL = "Update " + tblNm + " set geometry = MakePoint(" 
				//		" " + to_string<double>(o1.get_x()) + " , " + to_string<double>(o1.get_y()) + " , " 
				//		" " + to_string<s> (srid ) + " ); " ; // where " + keyFld + " = " + to_string<int> (id ) + " ;";
				//	if ( sqlite3_exec( db, qrySQL.c_str(), NULL, NULL, NULL ) != SQLITE_OK )
				//	{
				//		fprintf(stdout,"\n Geometry field MakePoint failed for SQL %s  \n Error :  %s \n", qrySQL.c_str(), sqlite3_errmsg(db)) ;
				//		logFile<<"Geometry field MakePoint update failed!" <<endl<<" Sql command : "<<qrySQL<<endl<<"Error : "<<to_string<const char *>(sqlite3_errmsg(db))<<endl;
				//	}

				//// if insert executed without error	
				//} else { // if insert casued error
				//	logFile<<"\n\tInsert SQL "<< qrySQL<< endl<< " error: " <<to_string<const char *>(sqlite3_errmsg(db));
				//	sqlite3_free(zErrMsg);
				//	return false;
				//}

				//strSQL = " SELECT recovergeometrycolumn('" + tblVertex + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Point',2);";
				//blnCreate = recoverSpatialGeometry(strSQL,tblVertex,db,intSRID,outfile);

				// create summary view for blocks / vertices 
				str1 = " SELECT  VertID, StopOn, sum(CostOn) OnCoSum, Count(StopOn) OnCnt, "
					"  StopOff, sum(CostOff) OffCoSum, Count(StopOff) OffCnt "
					" FROM " + tblVertex + " "  // RTD_BLOCKS_CDP_Rte12D5ai100004t8695498 
					" group by stopon, stopoff "
					" Order by stopon, stopoff;";
					viewName =  "VW_" + tblVertex ;
					blnCreate = createView(str1,viewName,db,outfile);
			}

		}
		
	// export the processed vertices
	//mapvert0 = VertexReport(mapvert0,pVx1,cout);
	// write the Edge data
		onoff = 1;
		//fileName(edgeinfname,"eg.bin",outfilename,onoff,r);
		strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(r) +"eg.bin";
		outfile.open(strEdgeBaseName.c_str(), ios::out | ios::binary|ios::trunc);
		writeBinObjectData(maped1,ev,ip,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

// read the edge data
		onoff = 0;

// write the result statistics into an output file
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

		str1 = "D" + to_string<short>(D) + "eg.out";
		//ext2[0]='\0';
		//strcpy(ext2, str1.c_str());
		//fileName(edgeinfname,ext2,outfilename,onoff);
		strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + str1;
		outfile.open(strOutFname.c_str(), ios::out |ios::app);

// edge voronoi process
		onoff=0; // off - from stop to parcel

//edv(mapvert0, mmaped0,mmapV1V2EgId,mmapEdgeStop,onoff);

		maped0 = edgeVoronoi (mapvert0, maped0,mmapV1V2EgId,mmapEdgeStop,vx,ev,onoff,q,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

		////fileName(edgeinfname,"voreg.bin",outfilename,onoff,r);
		//strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(r) + "voreg.bin";
		//outfile.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
		//writeBinObjectData(maped0,ev,ip,outfile);

		//if (outfile.is_open()) {
		//	outfile.close();
		//}
		//outfile.clear();

		//str1 = "D" + to_string<short>(D) + "voreg.txt";

		////ext2[0]='\0';
		////strcpy(ext2, str1.c_str());
		////fileName(edgeinfname,ext2,outfilename,onoff,q);
		//strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) +str1;
		//outfile.open(strOutFname.c_str(), ios::out |ios::trunc);
		//ev.serializetexthdr(outfile);
		//writeTextObjectData(maped0,ev,i1,outfile,"");

		//if (outfile.is_open()) {
		//	outfile.close();
		//}
		//outfile.clear();

// aggregate the edge results
// first remap the edge key by stop
		mapegor0.clear();
		mapegor0 = remapVertId2Orig(maped0,mapegor0,*pEg,i1);	

		if ((q%xi)==0 || blnHist  || q==1) {
			// write out the edge results
			str1 = "D" + to_string<short>(D) + ".egv";
			//fileName(edgeinfname,ext2,outfilename,onoff,q);
			strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) + str1;
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();

			outfile.open(strOutFname.c_str(),ios::out |ios::app);
			edgeSumAcc.clear();
	//	edgeSumAcc = aggbystop (stops,mapegor0,*pstop,*pEg,edgeAcc);
			str1 = "Summary of Edge attributes by stop id iteration " + to_string(q) + 
			" \nId\tStop\tecost\tsCost\ttCost\tParCnt\tlbl\tdirn";
			k1=-10;
			sit = stops.end();
			sit--;
			i1 =  sit->first;
			datetimeStamp(outfile);
			outfile << str1<<endl;

			edgeSumAcc = vertaggbyMaxKey (i1,k1,mapegor0,*pEg,edgeAcc);
			writeObjData(edgeSumAcc,edgeAcc,str1,outfile);
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();

	// write the result statistics into an output file
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();
			//fileName(edgeinfname,"eg.out",outfilename);
			strOutFname = strEdgeBaseName + "eg.out";
			outfile.open(strOutFname.c_str(), ios::out |ios::app);
		}
		// process the ons edge voronoi routine
		onoff=1;
		//edv (mapvert1, mmaped1,mmapV1V2EgId,mmapEdgeStop,onoff);
		maped1 = edgeVoronoi (mapvert1, maped1,mmapV1V2EgId,mmapEdgeStop,vx,ev,onoff,q,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

// aggregate the edge results
// first remap the edge key by stop
		mapegor1.clear();
		mapegor1 = remapVertId2Orig(maped1,mapegor1,*pEg,i1);	

		if ((q%xi)==0 || blnHist || q==1) {
			// write out the edge results
			str1 = "D" + to_string<short>(D) + ".egv";
			//ext2[0]='\0';
			//strcpy(ext2, str1.c_str());
			//fileName(edgeinfname,ext2,outfilename,onoff,q);
			strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) + str1;

			outfile.open(strOutFname.c_str(),ios::out |ios::app);
			edgeSumAcc.clear();
	//	edgeSumAcc = aggbystop (stops,mapegor0,*pstop,*pEg,edgeAcc);
			sit = stops.end();
			sit--;
			i1 =  sit->first;
			k1=-10;
			edgeSumAcc = vertaggbyMaxKey (i1,k1,mapegor1,*pEg,edgeAcc);
			str1 = "Summary of Edge attributes by stop id iteration " + to_string(q) + 
			" \nId\tStop\tecost\tsCost\ttCost\tParCnt\tlbl\tdirn";
			datetimeStamp(outfile);
			outfile << str1<<endl;
			writeObjData(edgeSumAcc,edgeAcc,str1,outfile);
			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();


			//fileName(edgeinfname,"voreg.bin",outfilename,onoff,r);
			strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(r) + "voreg.bin";
			outfile.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
			writeBinObjectData(maped1,ev,ip,outfile);

			if (outfile.is_open()) {
				outfile.close();
			}
			outfile.clear();

			//str1 = "D" + to_string<short>(D) + "voreg.txt";
			////ext2[0]='\0';
			////strcpy(ext2, str1.c_str());
			////fileName(edgeinfname,ext2,outfilename,onoff,q);
			//strOutFname = strEdgeBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) + str1;
			//outfile.open(strOutFname.c_str(), ios::out | ios::trunc);
			//ev.serializetexthdr(outfile);
			//writeTextObjectData(maped1,ev,i1,outfile,"");

			//if (outfile.is_open()) {
			//	outfile.close();
			//}
			//outfile.clear();
// create table of edge voronoi result for on and offs 

			str1 =	"( EdgeID INTEGER PRIMARY KEY , ESID INTEGER default -1, StopOn INTEGER default -1,"
				" CostOn DOUBLE default 0.0, lblOn INTEGER default 0, dirnOn Integer default -1, "
				" StopOff INTEGER default -1,CostOff DOUBLE default 0.0, lblOff INTEGER default 0," 
				"  dirnOff Integer default -1 ); ";
			intSRID = inpList.get_srid();
			if (blnHist) {
				tblEdge = inpList.get_tbledge() + "D" + to_string<short> (D) + "hi" + to_string<int> (q)+ "t" + to_string<long> (kStop.tripId()) + kStop.tripPeriod(); 
			} else 
			{
				tblEdge = inpList.get_tbledge() + "D" + to_string<short> (D) + "ai" + to_string<int> (q)+ "t" + to_string<long> (kStop.tripId()) + kStop.tripPeriod();  
			}
			replace(tblEdge.begin(),tblEdge.end(),' ','_');
			ReplaceAll2(tblEdge,"__","_");

			blnCreate = createSpaTbl(str1,tblEdge,db,intSRID,outfile);
			if (blnCreate) {
	//			strSQL = " SELECT recovergeometrycolumn('" + tblEdge + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'LineString',2);";
	//			blnCreate = recoverSpatialGeometry(strSQL,tblEdge,db,intSRID,outfile);
			// write the vedrtex output to the table
				keyFld = "EdgeID";
				string strSQLInserTbl = "Insert into " + tblEdge + " (\"EdgeId\",\"ESID\","  
					"\"StopOn\",\"CostOn\",\"lblOn\",\"dirnOn\"," 
					"\"StopOff\",\"CostOff\",\"lblOff\",\"dirnOff\" ) "  
					" Values ( ?, ?, ?, ?, ?, ?, ?, ?,?,? );" ;

			// call insert table routine to populate the vertex table
				// skip remaping the containers from ESID to EdgeID (id) and use ESID as the key between offs and ons and create table

				blnInsert = inSpaTblEdge(strSQLInserTbl,db,*pEg,maped1,maped0,tblEdge, keyFld , intSRID,logfile);
				// create summary view for blocks / vertices 
				str1 = " SELECT  EdgeID, StopOn, sum(CostOn) OnCoSum, Count(StopOn) OnCnt, "
					" StopOff, sum(CostOff) OffCoSum, Count(StopOff) OffCnt  "
					" FROM " + tblEdge + " "  // RTD_BLOCKS_CDP_Rte12D5ai100004t8695498 
					" group by StopOn, StopOff "
					" Order by StopOn, StopOff;";
					viewName =  "VW_" + tblEdge ;
					blnCreate = createView(str1,viewName,db,outfile);

			}
		}
		cout <<endl<< " End of street partitioning..."<<endl<<endl<<" Start Parcel to Stop Assignment ..."<<endl<<endl;
	} // if this is not euclidean metrics
 //open parcel data output file - empty file if it already exists before writing new data

	if (blnHist) {
// read the original parcel data from the parcel file with the edge id as the key for the multi-map
		strInFname = strParcelBaseName + ".bin" ;
		//fileName(pariname,".bin",infilename);
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		mmaparced0.clear();
		readBinParcObjectData(mmaparced0,par1,infile,blnHist);
	} else {
// Read the historic parcel data from the binary file format into the multi-map memory (mmaparced0)
		onoff = 1;
		//fileName(pariname,"out.bin",infilename,onoff,alt);
		strInFname = strParcelBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(alt) + "out.bin";
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		mmaparced0.clear();
		readBinParcObjectData(mmaparced0,par1,infile,blnHist);
	}

	if (infile.is_open()) {
		infile.close();
	}
	infile.clear();


//Assign Parcel to stops and update the walk time with the vertex/edge voronoi result for the onoff=0 (offs)
		onoff=0;

	if (!blnEuclid) {
		mmaped0.clear();
		mmaped0 = remapObj2Id(maped0,mmaped0,*pEg,i1);	

		mmaparStop.clear();

		mmaparStop = parcstopwalknew(mmapeidvid,mmaparced0,mapvert0, mmaped0,mmaparStop,onoff,blnHist,logfile);
	} else {
	// perform the euclid voronoi calculation
		mmapStopVx.clear();
		mmapVx1.clear();
		mapvert0.clear();
		mapvert1.clear();
		// for Euclidean matrix solve the network using stop and parcel data
		vxveuclid(mmaparced0,stops,mmapStopVx,mapvert0,mapvert1,mmapVx1,gc1,blnHist);
		mmaparStop.clear();
		mmaparStop = SortOnOffIdParcObject(mmaparced0,mmaparStop,onoff);	
	}

//calculate the summation for the onoff=0 (offs)
		ParcelOffAccum paoff;
		deque<ParcelOffAccum> vpaoff;
		deque<ParcelOffAccum>::iterator vpait;
 // calculate summary of parcel impact by tranist stop - (aggregate parcel walk time impact by transit stop)
//vpaoff = pagwalkstopoff (stops,mmaparStop,luci,gcost,onoff,paoff,blnHist);

		if (blnHist) {
			mmaparStop = pacalc (stops,mmaparStop,&luci,&gc1,onoff,blnHist);
		}

		vpaoff = pagstop (stops,mmaparStop,paoff);
		str1 = "D" + to_string<short>(D) + ".pac";
		//ext2[0]='\0';
		//strcpy(ext2, str1.c_str());
		//fileName(pariname,ext2,outfilename,onoff);
		strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + str1;
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();
		str1 = " Scenario id  " + to_string(q);
		if (q>0) {
			outfile.open(strOutFname.c_str(),ios::out |ios::app);
		} else {
			outfile.open(strOutFname.c_str(),ios::out |ios::trunc);
		}
		datetimeStamp(outfile);
		paoff.AccumHeader(outfile);
		writeObjData(vpaoff,paoff,str1,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

// for historic run update the offs, hist. walk time offs 
		updStopDmndWlkTm(vpaoff, stops, onoff,blnHist);

// calculate the parcel offs using the aggregate 
		if (blnHist) {
			calcParcelDemand0 ( stops, mmaparStop, onoff, blnHist);	
		}

// aggregate the parcel cost by transit stops
		vpaoff = pagstop (stops,mmaparStop,paoff); //pagwalkstop<maplngstop,mmaplngpar,ParcelOffAccum> (stops,mmaparStop,paoff);

		str1 = "D" + to_string<short>(D) + ".pac";
		//ext2[0]='\0';
		//strcpy(ext2, str1.c_str());
		//fileName(pariname,ext2,outfilename,onoff);
		strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + str1;
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();
		str1 = " Scenario id  " + to_string(q); //"Stop\tPVal \tAval \th.ons \th.In offs \tIn ons \tIn offs \tParCnt\n";
		datetimeStamp(outfile);
		if (q>0) { 
			outfile.open(strOutFname.c_str(),ios::out |ios::app);
		} else {
			outfile.open(strOutFname.c_str(),ios::out |ios::trunc);
		}
		paoff.AccumHeader(outfile);	
		writeObjData(vpaoff,paoff,str1,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

		// write the parcel data into binary file format
		//fileName(pariname,"out.bin",outfilename,onoff,r);
		strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + "_" + to_string<short>(r) + "out.bin";
		outfile.open(strOutFname.c_str(), ios::out | ios::binary | ios::trunc);
		writeBinObjectData(mmaparStop,par1,ip,outfile);
		if (outfile.is_open()) {
			outfile.close();
		}
		outfile.clear();

		//fileName(pariname,"out.bin",infilename,onoff,r);
		strInFname = strParcelBaseName + "_" + to_string<short>(onoff) + "_" + to_string<short>(r) + "out.bin";
		infile.open(strInFname.c_str(), ios::in | ios::binary);
		mmaparced1.clear();
		readBinParcObjectData(mmaparced1,par1,infile,blnHist);
		if (infile.is_open()) {
			infile.close();
		}
		infile.clear();

		onoff = 1;

// Read the parcel data from the binary file format into the multi-map memory (mmaparced1)
	if (!blnEuclid) {

//Assign Parcel to stops and update the walk time with the vertex/edge voronoi result for the onoff=1 (ons)

		mmaped1.clear();
		mmaped1 = remapObj2Id(maped1,mmaped1,*pEg,i1);	

		mmaparStop.clear();
		mmaparStop = parcstopwalknew(mmapeidvid,mmaparced1,mapvert1, mmaped1,mmaparStop,onoff,blnHist,logfile);
	} else {
		// copy map container
		mmaparStop.clear();
		mmaparStop = SortOnOffIdParcObject(mmaparced1,mmaparStop,onoff);
	}
	// calculate summary of parcel impact by transit stop - (aggregate parcel walk time impact by transit stop)
	ParcelOnAccum paon;
	deque<ParcelOnAccum> vpaon;
	deque<ParcelOnAccum>::iterator vpaonit;

	if (blnHist) {
		mmaparStop = pacalc (stops,mmaparStop,&luci,&gc1,onoff,blnHist);
	}

	vpaon = pagstop (stops,mmaparStop,paon);

	updStopDmndWlkTm(vpaon, stops, onoff,blnHist);

// calculate the parcel ons using the aggregate 
	if (blnHist) {
		mmaparStop = calcParcelDemand0 ( stops, mmaparStop, onoff, blnHist);	
	}

// reaggregate the parcel data as a check against the historic demand per stop 
//vpaon = pagwalkstop<maplngstop,mmaplngpar,ParcelOnAccum> (stops,mmaparStop,paon);
	vpaon = pagstop (stops,mmaparStop,paon);

	str1 = "D" + to_string<short>(D) + ".pac";
	//ext2[0]='\0';
	//strcpy(ext2, str1.c_str());
	//fileName(pariname,ext2,outfilename,onoff);
	strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + str1;
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	str1 = ""; //"Stop\tPVal \tAval \th.ons \th.In offs \tIn ons \tIn offs \tParCnt\n";
	if (q>0) {
		outfile.open(strOutFname.c_str(),ios::out |ios::app);
	} else {
		outfile.open(strOutFname.c_str(),ios::out |ios::trunc);
	}
	datetimeStamp(outfile);
	paon.AccumHeader(outfile);	
	writeObjData(vpaon,paon,str1,outfile);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

//open parcel data output file - empty file if it already exists before writing new data
	//fileName(pariname,"out.bin",outfilename,onoff,r);
	strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(r) + "out.bin";
	outfile.clear();
	outfile.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
// write the parcel data to file
	writeBinObjectData(mmaparStop,par1,ip,outfile);
// close parcel output file
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

//// write the parcel data to text file format file
//	str1 = "D" + to_string<short>(D) + "pass.txt";
//	//ext2[0]='\0';
//	//strcpy(ext2, str1.c_str());
//	//fileName(pariname,ext2,outfilename,onoff,q);
//	strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) + str1;
//	outfile.open(strOutFname.c_str(), ios::out | ios::trunc);
//	par1.serializetexthdr(outfile);
//	writeTextObjectData(mmaparStop,par1,i1,outfile,"");
//	if (outfile.is_open()) {
//		outfile.close();
//	}
//	outfile.clear();
	maparcid0.clear();
	// create assignment result parcel table for on and offs
		maparcid0 = remapObj2Id(mmaparStop,maparcid0,par1,i1);	
		if ((q%xi)==0 || blnHist || q==1) {
			//  the primary key requirement since SQLite seems to not allow update using statement
			str1 =	"( ID INTEGER PRIMARY KEY , Parcid text, StopOn INTEGER default -1, "
				" Ons Double default 0,CostOn DOUBLE default 0.0, lblOn INTEGER default 0,  "
				" StopOff INTEGER default -1,Offs Double default 0, CostOff DOUBLE default 0.0, " 
				" lblOff INTEGER default 0 ); ";
			intSRID = inpList.get_srid();
			if (blnHist) {
				tblParcel = inpList.get_tblparcel() + "D" + to_string<short> (D) + "hi" + to_string<int> (q)+ "t" + to_string<long> (kStop.tripId()) + kStop.tripPeriod();  
			} else {
				tblParcel = inpList.get_tblparcel() + "D" + to_string<short> (D) + "ai" + to_string<int> (q)+ "t" + to_string<long> (kStop.tripId()) + kStop.tripPeriod();  
			}
			replace(tblParcel.begin(),tblParcel.end(),' ','_');

			blnCreate = createSpaTbl(str1,tblParcel,db,intSRID,outfile);
			if (blnCreate) {
			// write the vedrtex output to the table
				keyFld = "ID";
				string strSQLInserTbl = "Insert into " + tblParcel + " (\"Id\",\"ParcId\","  
					"\"StopOn\",\"Ons\",\"CostOn\",\"lblOn\"," 
					"\"StopOff\",\"Offs\",\"CostOff\",\"lblOff\" ) "  
					" Values ( ?, ?, ?, ?, ?, ?, ?, ?,?,? );" ;

			// call insert table routine to populate the vertex table
				blnInsert = inSpaTblParcel(strSQLInserTbl,db,par1,maparcid0,tblParcel, keyFld , intSRID,logfile,blnHist);
				//strSQL = " SELECT recovergeometrycolumn('" + tblParcel + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Polygon',2);";
				//blnCreate = recoverSpatialGeometry(strSQL,tblParcel,db,intSRID,outfile);

				// create summary view for blocks  
				str1 = " SELECT  ID, Parcid, StopOn, sum(Ons) OnSum, sum(CostOn) OnCoSum,lblOn,count(StopOn) OnCnt, "
					"  StopOff, sum(Offs) OffSum, sum(CostOff) OffCoSum, lblOff,count(StopOff) OffCnt "
					" FROM " + tblParcel + " "  // RTD_BLOCKS_CDP_Rte12D5ai100004t8695498 
					" group by stopon, stopoff "
					" Order by stopon, stopoff;";
					viewName =  "VW_" + tblParcel ;
					blnCreate = createView(str1,viewName,db,outfile);

			}

		}

// write the stop assignment and parcel data to a binary format file
	//fileName(pariname,"stopid.bin",outfilename,onoff,r);
	strOutFname = strParcelBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(r) + "stopid.bin";
	outfile.open(strOutFname.c_str(), ios::out | ios::binary | ios::trunc);
	writeBinIdplusObjectData(mmaparStop,par1,outfile);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

	if (blnHist) {
//vpaon.clear();
 // calculate the undelayed run time
		tsidmit = stops.begin();
		tsidmit2 = stops.end();
//    undelTm<maplngstop::iterator>(tsidmit, tsidmit2,gc1.get_unitontm(),
//			gc1.get_unitofftm(), phway.get_hdway());
		stops = undelTripMap<maplngstop,maplngstop::iterator>(stops, gc1.get_unitontm(),
			gc1.get_unitofftm(), phway.get_hdway(),logfile);
	}

// loop over the stops and calculate the riding and operating details


	// for this run update the offs, walk time offs 
		onoff = 0;
		mmaparStop0.clear();
		SortOnOffIdParcObject(mmaparStop,mmaparStop0,onoff);
		vpaoff = pagstop (stops,mmaparStop0,paoff);
		stops = updStopDmndWlkTm(vpaoff, stops, onoff,blnHist);

	// for this run update the ons, walk time ons 
		onoff = 1;
		mmaparStop1.clear();
		mmaparStop1 = SortOnOffIdParcObject(mmaparStop,mmaparStop1,onoff);
		vpaon = pagstop (stops,mmaparStop1,paon);
		stops = updStopDmndWlkTm(vpaon, stops, onoff,blnHist);

//	ImpactCalc<maplngstop::iterator>(tsidmit,tsidmit2,gc1,phway);
		stops = ImpactCalcMap<maplngstop,globcost,pdhdway,int,maplngstop::iterator>(stops,gc1,phway,q);

//open stop data output file - empty file if it already exists before writing new data
		onoff=2;
		str1 = "D" + to_string<short>(D) + "out.txt";
		//ext2[0]='\0';
		//strcpy(ext2, str1.c_str());
		//fileName(tstopiname,str1,outfilename,onoff,q);
		strOutFname = strStopBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) + str1;

		ofstream stopout(strOutFname.c_str(), ios::out |ios::trunc);
// write the global header data
		gc1.show_globcosthdr(stopout);
		gc1.show_globcost(stopout);

// write the global header data
		phway.show_pdhdwayhdr(stopout);
		phway.show_pdhdway(stopout);
// write the stop data to file
		stop0.serializetexthdr(stopout);
		string txthdr = "";
		writeTextObjectData(stops,stop0,i1,stopout,txthdr);
// close parcel output file
		if (stopout.is_open()) {
			 stopout.close();
		}
		stopout.clear();

//open binary output stop set file - empty file if it already exists, write new data,
		onoff=2;
		//fileName(tstopiname,"alt.bin",outfilename,onoff,q);
		strOutFname = strStopBaseName + "_" + to_string<short>(onoff) + "_" + to_string<int>(q) + "alt.bin";
		stopout.open(strOutFname.c_str(), ios::out | ios::binary|ios::trunc);
		writeBinObjectData(stops,stop0,ip,stopout);
		if (stopout.is_open()) {
			stopout.close();
		}
		stopout.clear();
 
// create result table of Stops 
	    intSRID = inpList.get_srid();

		str1 =	"( Id Integer PRIMARY KEY, TripId Integer, StopId Integer, StOrdr Integer, "
				"StopName Text,StopIdp Integer,StopIds Integer,EdgeId Integer, palong Double, "
				"lbl Integer,Hist Integer ,Inbound Integer , External Integer ,Included Integer ,"
				" Eliminate Integer, CumDist Double, CRdTm Double, undCRdTm Double,CrdTmC Double, "
				" HistOns Double,HistOffs Double, Ons Double,Offs Double,DepVol Double, "
				" probStoph Double,probStop Double, depDelay Double, arrDelay Double, "
				" StopDelay Double, dwlDelay Double,rideDelay Double," 
				"PVal Double,AVal Double,CRdTmE Double,hWkTmOns Double,hWkTmOffs Double,WkTmOns Double,"
				"WkTmOffs Double,WalkCost Double,RideCost Double, OperCost Double, TCost Double,"
				"XC Double,YC Double ,ZC Double , Geometry Point ); ";
		if (blnHist) {
			tblTrip = inpList.get_tbltrip() + "D" + to_string<short> (D) + "hi" + to_string<int> (q) + "t" + to_string<long> (kStop.tripId())+ kStop.tripPeriod(); 
		} else 
		{
			tblTrip = inpList.get_tbltrip() + "D" + to_string<short> (D) + "ai" + to_string<int> (q) + "t" + to_string<long> (kStop.tripId())+ kStop.tripPeriod(); 
		}
		replace(tblTrip.begin(),tblTrip.end(),' ','_');
		ReplaceAll2(tblTrip,"__","_");

		blnCreate = createSpaTbl(str1,tblTrip,db,intSRID,outfile);
		if (blnCreate) {
			// Create Summary View for Result Table
			str1 = " SELECT count(id) NumStops, max(CRdTm) CRDTM, max(undCRdTm) UnCRDTM, max(CrdTmC) CRDTMC, "
				" sum(Ons) Ons, sum(Offs) offs, DepVol, sum((depDelay+ arrDelay)*probStop) StopDelay, " 
				" sum(dwlDelay) DwlDelay, sum(rideDelay) RDelay, sum(WkTmOns) WkTmOns, Sum(WkTmOffs) WkTmOffs, " 
				" sum(WalkCost) WkCost, sum(RideCost) RdCost, sum(OperCost) OpCost, sum(TCost) TCost, "
				" St_Collect(makePoint(XC,YC," + to_string<int>(intSRID) + "))  Geometry  "
				" FROM " + tblTrip + "  ;" ;
				viewName =  "VW_" + tblTrip ;
				blnCreate = createView(str1,viewName,db,outfile);
		// write the stop output to the table
			keyFld = "ID";
			string strSQLInserTbl = "Insert into " + tblTrip + " (\"ID\",\"TripId\","  
				"\"StopId\",\"StOrdr\",\"StopName\",\"StopIdp\",\"StopIds\",\"EdgeId\",\"palong\"," 
				"\"lbl\",\"Hist\",\"Inbound\",\"External\",\"Included\",\"Eliminate\"," 
				"\"CumDist\",\"CRdTm\",\"undCRdTm\",\"CrdTmC\",\"HistOns\",\"HistOffs\", "  
				"\"Ons\",\"Offs\",\"DepVol\",\"probStoph\",\"probStop\",\"depDelay\", "  
				"\"arrDelay\",\"StopDelay\",\"dwlDelay\",\"rideDelay\",\"PVal\",\"AVal\",\"CRdTmE\", "  
				"\"hWkTmOns\",\"hWkTmOffs\",\"WkTmOns\",\"WkTmOffs\",\"WalkCost\",\"RideCost\", "  
				"\"OperCost\",\"TCost\",\"XC\",\"YC\",\"ZC\",\"Geometry\") "  
				" Values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);" ;

		// call insert table routine to populate the stop table
			blnInsert = inSpaTblStop(strSQLInserTbl,db,*pstop,stops,tblTrip, keyFld , intSRID,logfile);
		
			if (blnInsert) {
				strSQL = " Update " + tblTrip + " set Geometry = MakePoint(XC,YC," + (to_string<int>(intSRID)) + " );";
				blnCreate = execSpatialQuery(strSQL,tblTrip,db,intSRID,outfile);
				
				if (blnInsert) {
					strSQL = " SELECT recovergeometrycolumn('" + tblTrip + "', 'Geometry'," + (to_string<int>(intSRID)) + " ,'Point',2);";
					blnCreate = recoverSpatialGeometry(strSQL,tblTrip,db,intSRID,outfile);
				}
			}
		}

		return stops;
}




vertexp& readvertex1(string& rec1, vertexp& vx, map<  int, string >::iterator maphdrit, char *seps) // *seps= "," 
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

  mapflds1 = recoread(rec1,seps);

  mapfldsit = mapflds1.begin();

while ( mapfldsit != mapflds1.end())
{
	i = mapfldsit->first;

	fldval = mapfldsit->second;
//    strx = fldval.c_str(); //strtok(stredge,seps);
	fldhdr = maphdrit->second;

   f1 = (fldhdr); 
   std::transform(f1.begin(), f1.end(), f1.begin(), to_upper());
   if (f1 ==  "VERTOID" || f1 == "VERTID" || f1 == "OBJECTID")
   {
      vid = fromString<long>(fldval); //strtol(strx,&stop1,ib);
	  if (vid>0) 
	  {
		  vx.set_id(vid);
	  }
	  else
	  {
		  return vx;
	  }
	}
	else if (f1 == "PVERTOID" || f1 == "IDP" )  //predecssor node id
	{
		vidp = fromString<long>(fldval); //strtol(strx,&stop1,ib);
	    vx.set_idp(vidp);
	}
	else if (f1 == "ORIGINOID")  //voronoi attractor/source/facility id
	{
		vp = fromString<long>(fldval); //strtol(strx,&stop1,ib);
			vx.set_idp((int) vp);
			vx.set_orig(vp);
	}
	else if (f1 == "TCOST") // total cost at node (the least cost path from the source to the node)
	{
		tcost = fromString<double>(fldval); //strtod(strx,&stop1);
		vx.set_tcost(tcost);
	}
	else if (f1 == "X" || f1 == "XC" || f1 == "XCOORD")  // X-Coordinate
	{
		  vcost = fromString<double>(fldval); //strtod(strx,&stop1);
//		  vx.pt.set_x(vcost);
	      vx.set_x(vcost);
	} 
	else if (f1 == "Y" || f1 == "YC" || f1 == "YCOORD")  // X-Coordinate
	{
		  vcost = fromString<double>(fldval); //strtod(strx,&stop1);
//		  vx.pt.set_y(vcost);
	      vx.set_y(vcost);
	} 
	else if (f1 == "Y" || f1 == "YC" || f1 == "YCOORD")  // X-Coordinate
	{
		  vcost = fromString<double>(fldval); //strtod(strx,&stop1);
//		  vx.pt.set_y(vcost);
	      vx.set_z(vcost);
	} 
	else if (f1 == "LABELED") // label 
	{
		if (vid>0) {vlbl = fromString<short>(fldval); //strtol(strx,&stop1,ib); 
		vx.set_lbl((short) vlbl);}
	}
             
   // get the next fld data
	mapfldsit++;
	maphdrit++;
} // close while
   if (vx.get_orig()>0)
   {
	cout<<"Stop point id "<<vid<<", origin : "<<vx.get_orig()<<", Cost "<<vcost<<endl;
   }
   return vx;
}

	list<string>  loadFields(sqlite3* netdb ,string mTableName,bool isQuery)
	{
	  int ret;
	  int i;
	  sqlite3_stmt *stmt = NULL;
	  char **results;
	  int rows;
	  int columns;
	  char *errMsg = NULL;
	  string pkName;
	  int pkCount = 0;
	  int fldNo = 0;
	  string sql;
	  list<string> attributeFields;
	  if (isQuery)
	  {
		mTableName = "\"" + mTableName + "\"";
		sql = "PRAGMA table_info( " +  mTableName + "; ";

		ret = sqlite3_get_table( netdb, sql.c_str(), &results, &rows, &columns, &errMsg );
		if ( ret != SQLITE_OK )
		  goto error;
		if ( rows < 1 )
		  ;
		else
		{
		  for ( i = 1; i <= rows; i++ )
		  {
			string name = to_string( results[( i * columns ) + 1] );
			const char *type = results[( i * columns ) + 2];
			const char *pk = results[( i * columns ) + 5];
			if ( from_string<char*> (pk) != 0 )
			{
			  // found a Primary Key column
			  pkCount++;
			  pkName = name;
			}

			 // for sure any SQLite value can be represented as SQLITE_TEXT
			  string fieldType ;

			  // making some assumptions in order to guess a more realistic type
			  if ( strcasecmp( type, "int" ) == 0 ||
				   strcasecmp( type, "integer" ) == 0 ||
				   strcasecmp( type, "bigint" ) == 0 ||
				   strcasecmp( type, "smallint" ) == 0 ||
				   strcasecmp( type, "tinyint" ) == 0 ||
				   strcasecmp( type, "boolean" ) == 0 )
			  {
				fieldType = "int";
			  }
			  else if ( strcasecmp( type, "real" ) == 0 ||
						strcasecmp( type, "double" ) == 0 ||
						strcasecmp( type, "double precision" ) == 0 ||
						strcasecmp( type, "float" ) == 0 )
			  {
				fieldType = "double";
			  }

			  attributeFields.push_back(name);
		  }
		}
		sqlite3_free_table( results );
	  }
	  else
	  {
		sql =  "select * from " + mTableName + " limit 1";

		if ( sqlite3_prepare_v2( netdb, sql.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
		{
		  // some error occurred
		  fprintf(stderr, "SQLite error: %s\n\nSQL: %s", sql ,  sqlite3_errmsg( netdb) );
		  return attributeFields;
		}

		ret = sqlite3_step( stmt );
		if ( ret == SQLITE_DONE )
		{
		  // there are no rows to fetch
		  sqlite3_finalize( stmt );
		  return attributeFields;
		}

		if ( ret == SQLITE_ROW )
		{
		  // one valid row has been fetched from the result set
		  columns = sqlite3_column_count( stmt );
		  for ( i = 0; i < columns; i++ )
		  {
			string name = to_string<const char*> ( sqlite3_column_name( stmt, i ) );
			const char *type = sqlite3_column_decltype( stmt, i );
			if ( type == NULL )
			  type = "TEXT";

			  // for sure any SQLite value can be represented as SQLITE_TEXT
			  string fieldType;

			  // making some assumptions in order to guess a more realistic type
			  if ( strcasecmp( type, "int" ) == 0 ||
				   strcasecmp( type, "integer" ) == 0 ||
				   strcasecmp( type, "bigint" ) == 0 ||
				   strcasecmp( type, "smallint" ) == 0 ||
				   strcasecmp( type, "tinyint" ) == 0 ||
				   strcasecmp( type, "boolean" ) == 0 )
			  {
				fieldType = "int";
			  }
			  else if ( strcasecmp( type, "real" ) == 0 ||
						strcasecmp( type, "double" ) == 0 ||
						strcasecmp( type, "double precision" ) == 0 ||
						strcasecmp( type, "float" ) == 0 )
			  {
				fieldType = "double";
			  }

			  attributeFields.push_back(name );
		  }
		}
		sqlite3_finalize( stmt );
	  }


	  return attributeFields;

	error:
	  // unexpected error
	  if ( errMsg != NULL )
	  {
		fprintf(stderr, "SQL error: %1 " , errMsg  );
		sqlite3_free( errMsg );
	  }
	  return attributeFields;

	}



