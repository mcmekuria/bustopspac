using namespace std;

class kTrip {
	public:
	// trip seq, RouteId, TripId, DirnId 
	

		kTrip(long id1,  string route1, string trip1, string dirn1 ) ;

		kTrip( int id1, string route1, string trip1, string dirnId1);

		kTrip() : mid(0),  mRoute(""), mTrip(""), mDirn("")
			 {};
		kTrip(const kTrip& ts )  
		{
			kTrip::mid = ts.mid;
			kTrip::mRoute = ts.mRoute;
			kTrip::mTrip= ts.mTrip; 
			kTrip::mDirn= ts.mDirn; 
			//  cout << "Ts[" << id << "]" << endl;
			//    ++copycons;
	  }
	  kTrip& operator=(const kTrip& ts) 
	  {
		//cout << "(" << id << ")=[" << ts.mid << "]" << endl;
		kTrip::mid = ts.mid;
		kTrip::mRoute = ts.mRoute;
		kTrip::mTrip= ts.mTrip; 
		kTrip::mDirn= ts.mDirn; 
		return *this;
	  }

		bool operator==(kTrip ts);

	  // Create a stop from a stop pointer:
	  kTrip(const kTrip* ts, const int& ix) 
		: mid(ts->mid+ix) {
		cout << "Copied stop " << *this << " from "
			 << *ts << endl;
	  }

	 
	   void setId (long id1) {mid=id1;}

	   long id () {return mid;}

		void setRoute (string Route1) {mRoute=Route1;}

		string route () {return mRoute;}

		void setTrip (string trip1) {mTrip=trip1;}

		string trip () {return mTrip;}

		void setDirn (string dirn1) {mDirn=dirn1;}

		string dirn () {return mDirn;}

		friend std::ostream& operator<<(
			std::ostream& os, const kTrip& ts) {
				return os << ts.mid << "\t" << ts.mRoute << "\t"<< ts.mTrip <<"\t"<<ts.mDirn<<endl ;
	  }

		void ts_hdr(ostream& out)
		{ 
		 out << "id"<< "\t" << "Route"<<"\t" << "Trip Id" <<"\t"<< "Dirn"<<endl ;
	  }
	  
		void show_tsf(ostream& out)
		{ 
			out << id() << "\t" << route()<< "\t"  <<trip() << "\t" << dirn() <<"\t"<< endl;
		};


		void show_tsdh(void)
		{ 
			cout << "Stop Id: " << id() << endl;
			cout << "Direction: " << dirn() << endl;
		};

		void show_tsd(void)
		{ 
			cout << "Stop Id: " << id() << endl;
			cout << "Route : " << route() << endl;
			cout << "Trip Id: " << trip() << endl;
			cout << "Dirn : " << dirn() << endl;
		};

			void kTrip::serialize(ofstream& pts);
			void kTrip::serializetext(ofstream& pts);
			void kTrip::serializetexthdr(ofstream& pts);
			void kTrip::deserialize(ifstream& pts);
			void kTrip::deserializetext(ifstream& pts);
//			string kTrip::makey( kTrip& kt);

		~kTrip(){
				//cout << "Deleting stop: " << *this << endl;
		}

	protected:
				long mid; // Stop id j
				string mRoute;  // Route
				string mTrip; // Trip Id
				string mDirn; // direction 
	};

	void kTrip::serialize(ofstream& pts)
	{
	 pts.write(reinterpret_cast<char *>(&mid), sizeof(mid));
	 streamsize sizet=mRoute.size();// store Route ID's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mRoute.c_str(), sizet+1); // write Route Id including '\0' too
	 sizet=mTrip.size();// store Trip ID's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mTrip.c_str(), sizet+1); // write Trip Id including '\0' too
	 sizet=mDirn.size();// store Trip ID's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mDirn.c_str(), sizet+1); // write Trip Id including '\0' too

	}


	void kTrip::deserialize(ifstream& pts)
	{
	 int len=0;
	 char *p=0;
	 pts.read(reinterpret_cast<char *>(&mid), sizeof(mid));

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mRoute=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mTrip=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mDirn=p; // copy temp to data member
	 delete[] p;
	}


	kTrip::kTrip(long id1, string route1="",
				 string trip1="", 
				 string dirn1="" ) 
	{
		kTrip::mid = id1;
		kTrip::mRoute = route1;
		kTrip::mTrip = trip1;
		kTrip::mDirn=dirn1; 
	}


	void kTrip::serializetext(ofstream& pts)
	{
	 pts <<mid<<"\t"<<mRoute<<"\t"<<mTrip<<"\t"<<mDirn<<endl;
	}
	void kTrip::serializetexthdr(ofstream& pts)
	{
	 pts <<"mid"<<"\t"<<"mRoute"<<"\t"<<"mTrip"<<"\t"<<"mDirn"<<"\t"<<endl;
	}

	void kTrip::deserializetext(ifstream& pts)
	{
	 pts >>mid>>mRoute>>mTrip>>mDirn;
	}

class stopKey {
	public:
	// RteName,timePeriod,TripSTime,TripKey,tripNumber, count(stop_Id) 
	

		stopKey( string route1,string schlName1,string dir1, string timePd1,string tripSTime1,double tripTime1, long tripKey1 , long tripNum1, long stopCnt1, double timeNexTrip1,double hOnSum1,double hOffSum1 ) ;


		//stopKey() :  mRoute(""), mtripPeriod(""), mtripSTime(""), mtripId(-1), mtripNumber(-1), mstopCount(-1)
		//	 {};
		stopKey(const stopKey& ts )  
		{
			stopKey::mRoute = ts.mRoute;
			stopKey::mschlName= ts.mschlName; 
			stopKey::mtripPeriod= ts.mtripPeriod; 
			stopKey::mtripSTime= ts.mtripSTime; 
			stopKey::mtripTime= ts.mtripTime; 
			stopKey::mtripId= ts.mtripId; 
			stopKey::mtripNumber= ts.mtripNumber; 
			stopKey::mstopCount= ts.mstopCount; 
			stopKey::mtimeNexTrip= ts.mtimeNexTrip; 
			stopKey::mhOnSum= ts.mhOnSum; 
			stopKey::mhOffSum= ts.mhOffSum; 
			//  cout << "Ts[" << id << "]" << endl;
			//    ++copycons;
	  }
	  stopKey& operator=(const stopKey& ts) 
	  {
		//cout << "(" << id << ")=[" << ts.mid << "]" << endl;
			stopKey::mRoute = ts.mRoute;
			stopKey::mschlName= ts.mschlName; 
			stopKey::mtripPeriod= ts.mtripPeriod; 
			stopKey::mtripSTime= ts.mtripSTime; 
			stopKey::mtripTime= ts.mtripTime; 
			stopKey::mtripId= ts.mtripId; 
			stopKey::mtripNumber= ts.mtripNumber; 
			stopKey::mstopCount= ts.mstopCount; 
			stopKey::mtimeNexTrip= ts.mtimeNexTrip; 
			stopKey::mhOnSum= ts.mhOnSum; 
			stopKey::mhOffSum= ts.mhOffSum; 
		return *this;
	  }

		bool operator==(stopKey ts);

	  // Create a stopkey object from a stopkey pointer:
	  stopKey(const stopKey* ts, const int& ix) 
		{
			cout << "Copied stop " << *this << " from "
			 << *ts << endl;
		}

	 
	   void settripId (long tripId1) {mtripId=tripId1;}

	   long tripId () {return mtripId;}

		void setRoute (string Route1) {mRoute=Route1;}

		string route () {return mRoute;}

		void setschlName (string schlName1) {mschlName=schlName1;}

		string schlName () {return mschlName;}

		void setdir (string dir1) {mdir=dir1;}

		string dir () {return mdir;}

		void settripPeriod (string tripPeriod1) {mtripPeriod=tripPeriod1;}

		string tripPeriod () {return mtripPeriod;}

		void settripSTime (string tripSTime1) {mtripSTime=tripSTime1;}

		string tripSTime () {return mtripSTime;}

		void settripTime (double tripTime1) {mtripTime=tripTime1;}

		double tripTime () {return mtripTime;}

		void settripNumber (long tripNumber1) {mtripNumber=tripNumber1;}

		long tripNumber () {return mtripNumber;}

		void setstopCount (long stopCount1) {mstopCount=stopCount1;}

		long stopCount () {return mstopCount;}

		void settimeNexTrip (long timeNexTrip1) {mtimeNexTrip=timeNexTrip1;}

		long timeNexTrip () {return mtimeNexTrip;}

		void sethOnSum (double hOnSum1) {mhOnSum=hOnSum1;}

		double hOnSum () {return mhOnSum;}

		void sethOffSum (double hOffSum1) {mhOffSum=hOffSum1;}

		double hOffSum () {return mhOffSum;}

		friend std::ostream& operator<<(
			std::ostream& os, const stopKey& ts) {
				return os << ts.mRoute << "\t"<<ts.mdir << "\t" <<ts.mtripPeriod <<"\t"<<ts.mtripSTime<< "\t"
					<<ts.mtripSTime<< "\t"<< ts.mtripId << "\t" << ts.mtripNumber << "\t"
					<< ts.mstopCount<<endl ;
	  }

		void ts_hdr(ostream& out)
		{ 
		 out << "Route"<<"\t"<<"Direction" << "\t" "Trip Pd" <<"\t" << "Trip Start Time" <<"\t"<< " Trip Id"
			 <<"\t"<< " Trip Number"<<"\t"<< " Stop Count"<<endl ;
	  }
	  
		void show_tsf(ostream& out)
		{ 
			out << route() << "\t"<<dir() << "\t"<< tripPeriod() <<"\t"<<tripSTime()<< "\t"<< tripId() 
					<< "\t" << tripNumber() << "\t"<< stopCount()<<endl ;
		};


		void show_tsdh(void)
		{ 
			cout << "Route : " << route()<< "\t Dirn "<<dir() ;
			cout << " Num. Stops : " << stopCount() << endl;
		};

		void show_tsd(void)
		{ 
			cout << "Route : " << route() << "\t Direction : " << dir() ;
			cout << "\tTrip Period: " << tripPeriod() ;
			cout << "\tNum. Stops : " << stopCount() << endl;
		};

			void stopKey::serialize(ofstream& pts);
			void stopKey::serializetext(ofstream& pts);
			void stopKey::serializetexthdr(ofstream& pts);
			void stopKey::deserialize(ifstream& pts);
			void stopKey::deserializetext(ifstream& pts);
//			string stopKey::makey( stopKey& kt);

		~stopKey(){
				//cout << "Deleting stop: " << *this << endl;
		}

	protected:
				string mRoute;  // Route
				string mschlName; // Trip Scheldule name (AM,PM, etc)
				string mtripPeriod; // Trip period
				string mdir; // Trip direction
				string mtripSTime; // Trip Start time 
				double mtripTime; // Trip Start time 
				long mtripId; // trip id (key) 
				long mtripNumber; // trip sequence number 
				long mstopCount; // number of stope in trip 
				double mtimeNexTrip; // time till next trip (headway) 
				double mhOnSum; // Historic Sum of the Boardings 
				double mhOffSum; // Historic Sum of the Alightings 
	};

	void stopKey::serialize(ofstream& pts)
	{
	 streamsize sizet=mRoute.size();// store Route ID's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mRoute.c_str(), sizet+1); // write Route Id including '\0' too

	 sizet=mschlName.size();// store schedule Name's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mschlName.c_str(), sizet+1); // write schedule Name including '\0' too

	 sizet=mdir.size();// store Route direction string length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mdir.c_str(), sizet+1); // write Route direction including '\0'

	 sizet=mtripPeriod.size();// store trip period's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mtripPeriod.c_str(), sizet+1); // write trip period including '\0' too

	 sizet=mtripSTime.size();// store start time's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mtripSTime.c_str(), sizet+1); // write start time including '\0' too

	 pts.write(reinterpret_cast<char *>(&mtripTime),sizeof(mtripTime));
	 pts.write(reinterpret_cast<char *>(&mtripId),sizeof(mtripId));
	 pts.write(reinterpret_cast<char *>(&mtripNumber), sizeof(mtripNumber));
	 pts.write(reinterpret_cast<char *>(&mstopCount), sizeof(mstopCount));
	 pts.write(reinterpret_cast<char *>(&mtimeNexTrip), sizeof(mtimeNexTrip));
	 pts.write(reinterpret_cast<char *>(&mhOnSum), sizeof(mhOnSum));
	 pts.write(reinterpret_cast<char *>(&mhOffSum), sizeof(mhOffSum));

	}


	void stopKey::deserialize(ifstream& pts)
	{
	 int len=0;
	 char *p=0;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mRoute=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for schedule Name 
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mschlName=p; // copy temp to data member (schedule Name )
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mdir=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mtripPeriod=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mtripSTime=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&mtripTime), sizeof(mtripTime));
	 pts.read(reinterpret_cast<char *>(&mtripId), sizeof(mtripId));
	 pts.read(reinterpret_cast<char *>(&mtripNumber),sizeof(mtripNumber));
	 pts.read(reinterpret_cast<char *>(&mstopCount),sizeof(mstopCount));
	 pts.read(reinterpret_cast<char *>(&mtimeNexTrip),sizeof(mtimeNexTrip));
	 pts.read(reinterpret_cast<char *>(&mhOnSum),sizeof(mhOnSum));
	 pts.read(reinterpret_cast<char *>(&mhOffSum),sizeof(mhOffSum));
	}


	stopKey::stopKey(string route1="",
				 string schlName1="", 
				 string dir1="", 
				 string tripPeriod1="", 
				 string tripSTime1="", 
				 double tripTime1=0, 
				 long tripId1=-1, 
				 long tripNumber1=-1, 
				 long stopCount1=-1,
				 double timeNexTrip1 = 0,
				 double hOnSum1 = 0,
				 double hOffSum1 = 0) 
	{
		stopKey::mRoute = route1;
		stopKey::mschlName = schlName1;
		stopKey::mdir = dir1;
		stopKey::mtripPeriod = tripPeriod1;
		stopKey::mtripSTime = tripSTime1;
		stopKey::mtripTime = tripTime1;
		stopKey::mtripPeriod = tripPeriod1;
		stopKey::mtripNumber = tripNumber1;
		stopKey::mstopCount=stopCount1; 
		stopKey::mtimeNexTrip=timeNexTrip1; 
		stopKey::mhOnSum= hOnSum1; 
		stopKey::mhOffSum= hOffSum1; 
	}


	void stopKey::serializetext(ofstream& pts)
	{
	 pts <<mRoute<<"\t"<<mdir<<"\t"<<mtripPeriod<<"\t"<<mtripId<<"\t"<<mtripNumber<<"\t"<<mstopCount<<endl;
	}
	void stopKey::serializetexthdr(ofstream& pts)
	{
	 pts <<"Route"<<"\t"<<"Direction"<<"\t"<<"TripPd"<<"\t"<<"TripId"<<"\t"<<"TripNumber"<<"\t"<<"StopCount"<<"\t"<<endl;
	}

	void stopKey::deserializetext(ifstream& pts)
	{
	 //pts >>mRoute>>mtripPeriod>>mtripId>>mtripNumber>>mstopCount>>endl;
	}

// Light Stop class for the purpose of mathcing stopId and StOrdr  to Stop Index numbers 

class ltStop 
{
private:
	long i,j,k;
	string dpkey;
	double dFract ;
public:
	ltStop() : i(-1),j(-1),k(-1) {};
	ltStop(const ltStop& ltstp )  {
		ltStop::i = ltstp.i;
		ltStop::j = ltstp.j;
		ltStop::k = ltstp.k;
		ltStop::dpkey = ltstp.dpkey;
		ltStop::dFract = ltstp.dFract;
	}
//  ltStop(long i1=-1, long j1=-1, long k1=-1)
//    : dpkey(makey(i1,j1,k1)),i(i1),j(j1), k(k1)  {
      //cout << "Creating DPStop Object from a stop object: " << *this << endl;
//  }


 /* ltStop(long i1=-1, long j1=-1, long k1=-1,long l1=-1,long m1=-1)
    : i(i1), j(j1), k(k1), l(l1), m(m1) {
		//cout << "Creating DPStop Object with out a stop object: " << *this << endl;
  }
  */
 ltStop& operator=(ltStop& ltstp) {
    ltStop::i = ltstp.get_i();
    ltStop::j = ltstp.get_j();
    ltStop::k = ltstp.get_k();
    ltStop::dpkey = ltstp.get_dpkey();
    ltStop::dFract = ltstp.get_dFract();

    return *this;
  }
  friend ostream&
  operator<<(ostream& os, const ltStop& ts) {
    return os <<ts.i << "\t" << ts.j << "\t" << ts.k<< "\t" << ts.dFract<< "\t" << ts.dpkey<<endl;
  }


 void ltStop::serializetext(ofstream& pts)
{
 pts <<i<<"\t"<<j<<"\t"<<k<<"\t"<<dFract<<"\t"<<dpkey;
}

void ltStop::serializetext3d(ofstream& pts)
{
 pts <<i<<"\t"<<j<<"\t"<<k<<"\t"<<dFract<<"\t"<<dpkey;
}

void ltStop::serializetext2(ofstream& pts)
{
 pts <<i<<"\t"<<j<<"\t"<<k<<"\t"<<dFract<<"\t"<<dpkey;
}

void ltStop::serializetext2d3(ofstream& pts)
{
 pts <<i<<"\t"<<j<<"\t"<<k<<"\t"<<dFract<<"\t"<<dpkey;
}

void ltStop::serializetexthdr(ofstream& pts)
{
 pts<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"dFract"<<"\t"<<"dpkey"<<"\t";
}

void ltStop::serializetext3dhdr(ofstream& pts)
{
 pts<<"i"<<"\t"<<"j"<<"\t"<<"dFract"<<"\t"<<"dpkey";
}

void ltStop::deserializetext(ifstream& pts)
{
 pts >>i>>j>>k>>dFract>>dpkey;
}

void ltStop::deserializetext3d(ifstream& pts)
{
 pts >>i>>j>>k>>dFract>>dpkey;
}

void ltStop::serialize(ofstream& pts)
{
 pts.write(reinterpret_cast<char *>(&i), sizeof(i));
 pts.write(reinterpret_cast<char *>(&j), sizeof(j));
 pts.write(reinterpret_cast<char *>(&k), sizeof(k));
 pts.write(reinterpret_cast<char *>(&dFract), sizeof(dFract));

 streamsize sizet=dpkey.size();// store dp key's length
 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
 pts.write(dpkey.c_str(), sizet+1); // write string dpkey and final '\0'
}


void ltStop::deserialize(ifstream& pts)
{
 int len=0;
 char *p=0;

 pts.read(reinterpret_cast<char *>(&i), sizeof(i));
 pts.read(reinterpret_cast<char *>(&j), sizeof(j));
 pts.read(reinterpret_cast<char *>(&k), sizeof(k));
 pts.read(reinterpret_cast<char *>(&dFract), sizeof(dFract));

 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
 p=new char [len+1]; // allocate temp buffer for stop Label
 pts.read(p, len+1); // copy string to temp, including '\0'
 dpkey=p; // copy temp to data member
 delete[] p;
}


virtual ~ltStop() { 
//	  delete tsem;
    //cout << "Deleting dpstop: " << *this << endl;
  }

	void set_i (long i1) {i=i1;}
	long get_i () {return i;}
	void set_j (long j1) {j=j1;}
	long get_j () {return j;}
	void set_k (long k1) {k=k1;}
	long get_k () {return k;}
	void set_dFract (double dF1) {dFract=dF1;}
	double get_dFract () {return dFract;}
	void set_dpkey(string dpkey1) {dpkey=dpkey1;}
	string get_dpkey () {return dpkey;}

	string makey(long i1,long j1,long k1,string sp="|")
	{
		return (to_string(i1) + sp + to_string(j1) + sp + to_string(k1) );
	}
	string makey(string i1,string j1,string k1,string sp="|")
	{
		return ((i1) + sp + (j1) + sp + (k1) );
	}

	string makey(ltStop& ltstp,char *sp="|")
	{
		return (to_string(ltstp.get_i()) + *sp + to_string(ltstp.get_j()) + *sp + to_string(ltstp.get_k()) );
	}

};


