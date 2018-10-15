using namespace std;

class sTrip {
	public:
//SELECT Id, SchlName, RteName, DirName, TimePeriod, TripSTime, TripKey, TripNumber, SchlTime, 
//AcTimeArr, AcTimeDep, SchlRunTm, StOrdr, Stop_Id, StrMainCrs, Ons, Offs, Lat, Lon, Geometry
//FROM RTD_Trips
	

		sTrip(long id1, string schlName1, string routeName1, string dirName1, string timePeriod1,
			string tripSTime1, long inTripKey, int tripNumber1,string schlTime1, string AcTimeArr,
			string acTimeDep1, int schlRunTm1 ,  int stOrdr1, long stopCode1, string strMainCrs1, 
				 double ons1, double offs1, double xc1 , double yc1);

		sTrip(long id1, string schlName1, string routeName1, string dirName1, string timePeriod1,
			string tripSTime1="", long inTripKey=-1, int tripNumber1=-1,string schlTime1="", 
			string AcTimeArr="", string acTimeDep1="", int schlRunTm1=-1 ,  string strMainCrs1="",
			int stOrdr1=-1, long stopCode1=-1, double ons1=0, double offs1=0, double xc1=-1 , 
			double yc1=-1, double zc1=-1);

		sTrip() : mid(0),   mschlName(""), mrouteName(""),mdirName(""),
			mtimePeriod(""), mtripSTime(""),mtripKey(-1), mtripNumber(-1),mschlTime(""),
			macTimeArr(""),macTimeDep(""),mschlRunTm(0),mstOrdr(-1), mstopCode(-1),
			mstopName(""), mstrMainCrs(""), mons (0), moffs(0), mxc(0),myc(0),mzc(0)
			 {};
		sTrip(const sTrip& ts )  
		{
			sTrip::mid = ts.mid;
			sTrip::mrouteName = ts.mrouteName;
			sTrip::mtripKey= ts.mtripKey; 
			sTrip::mtripNumber= ts.mtripNumber; 
			sTrip::mtimePeriod= ts.mtimePeriod; 
			sTrip::mtripSTime= ts.mtripSTime; 
			sTrip::mschlTime= ts.mschlTime; 
			sTrip::macTimeArr= ts.macTimeArr; 
			sTrip::macTimeDep= ts.macTimeDep; 
			sTrip::mschlRunTm= ts.mschlRunTm; 
			sTrip::mstrMainCrs= ts.mstrMainCrs; 
			sTrip::mstOrdr = ts.mstOrdr;
			sTrip::mstopCode = ts.mstopCode;
			sTrip::mons= ts.mons;
			sTrip::moffs= ts.moffs;
			sTrip::mxc = ts.mxc; 
			sTrip::myc = ts.myc; 
			sTrip::mzc = ts.mzc; 
			//  cout << "Ts[" << id << "]" << endl;
			//    ++copycons;
	  }
	  sTrip& operator=(const sTrip& ts) 
	  {
		//cout << "(" << id << ")=[" << ts.mid << "]" << endl;
			sTrip::mid = ts.mid;
			sTrip::mrouteName = ts.mrouteName;
			sTrip::mdirName = ts.mdirName;
			sTrip::mtripKey= ts.mtripKey; 
			sTrip::mtripNumber= ts.mtripNumber; 
			sTrip::mtimePeriod= ts.mtimePeriod; 
			sTrip::mtripSTime= ts.mtripSTime; 
			sTrip::mschlTime= ts.mschlTime; 
			sTrip::macTimeArr= ts.macTimeArr; 
			sTrip::macTimeDep= ts.macTimeDep; 
			sTrip::mschlRunTm= ts.mschlRunTm; 
			sTrip::mstOrdr = ts.mstOrdr;
			sTrip::mstopCode = ts.mstopCode;
			sTrip::mstrMainCrs= ts.mstrMainCrs; 
			sTrip::mons= ts.mons;
			sTrip::moffs= ts.moffs;
			sTrip::mxc = ts.mxc; 
			sTrip::myc = ts.myc; 
			sTrip::mzc = ts.mzc; 
		return *this;
	  }

		bool operator==(sTrip ts);

	  // Create a stop from a stop pointer:
	  sTrip(const sTrip* ts, const int& ix) 
		: mid(ts->mid+ix) {
		cout << "Copied stop " << *this << " from "
			 << *ts << endl;
	  }

	 
	   void setId (long id1) {mid=id1;}

	   long id () {return mid;}

		void setStOrdr (int StOrdr1) {mstOrdr=StOrdr1;}

		int stOrdr () {return mstOrdr;}  // stop order

		void setStopCode (long StopCode1) {mstopCode=StopCode1;}

		long stopCode () {return mstopCode;}

		void settripKey (long tripKey1) {mtripKey=tripKey1;}

		long tripKey () {return mtripKey;}

		void settripNumber (int tripNumber1) {mtripNumber=tripNumber1;}

		int tripNumber () {return mtripNumber;}

		void setschlRunTm (int schlRunTm1) {mschlRunTm=schlRunTm1;}

		int schlRunTm () {return mschlRunTm;}

		void setschlTime (string schlTime1) {mschlTime=schlTime1;}

		string schlTime () {return mschlTime;}

		void setstrMainCrs (string strMainCrs1) {mstrMainCrs=strMainCrs1;}

		string strMainCrs () {return mstrMainCrs;}

		void setacTimeArr (string acTimeArr1) {macTimeArr=acTimeArr1;}

		string acTimeArr () {return macTimeArr;}

		void setacTimeDep (string acTimeDep1) {macTimeArr=acTimeDep1;}

		string acTimeDep () {return macTimeDep;}

		void setons (double ons1) {mons=ons1;}

		double ons () {return mons;}

		void setoffs (double offs1) {mons=offs1;}

		double offs () {return moffs;}

		void setrouteName (string mrouteName1) {mrouteName=mrouteName1;}

		string routeName () {return mrouteName;}

		void setstopName (string stopName1) {mstopName = stopName1;}

		string stopName () {return mstopName;}

		void setdirName (string dirName1) {mdirName = dirName1;}

		string dirName () {return mdirName;}

		void setxc (double xc1) {mxc=xc1;}

		double xc () const {return mxc;}

		void setyc (double yc1) {myc=yc1;}

		double yc () const {return myc;}

		void setzc (double zc1) {mzc=zc1;}

		double zc () const {return mzc;}

		friend std::ostream& operator<<(
			std::ostream& os, const sTrip& ts) {
				return os << ts.mid << "\t" << ts.mstOrdr<< "\t"<< ts.mrouteName << "\t"<< ts.mstopCode
					<<"\t"<< ts.mstopName <<"\t"<< ts.mdirName<<"\t"<< ts.mtripKey <<"\t" 
					<< ts.mtripNumber <<  "\t"<< ts.macTimeArr<<"\t" << ts.macTimeDep<<"\t"
					<< ts.mschlRunTm<<"\t"<< ts.mstrMainCrs<<"\t" << ts.mons<<"\t" << ts.moffs 
					<<"\t"<<ts.mxc<<"\t"<<ts.myc<<"\t"<<ts.mzc<<"\t"<<endl ;
	  }

		void ts_hdr(ostream& out)
		{ 
		 out << "id"<< "\t" << "StopCode"<<"\t"<< "Stop Name" <<"\t"
		  << "Route"<<"\t" << "Parent Stn" <<"\t"<< "Dirn"<<"\t"<< "Trip Id " <<"\t"<< "Trip Shape"
		  <<"\t"<< "Trip Block"<<"\t" << "Trip Head"<<"\t"<<"XC"<<"\t"<<"YC"<<"\t"<<"ZC"<<"\t"<<endl ;
	  }
	  
		void show_tsf(ostream& out)
		{ 
			out << id() << "\t" << stOrdr()<< "\t"<< routeName()<< "\t"<< stopCode()
					<<"\t"<< stopName() <<"\t"<< dirName()<<"\t"<< tripKey() <<"\t" 
					<< tripNumber() <<  "\t"<< acTimeArr()<<"\t" << acTimeDep()<<"\t"
					<< schlRunTm()<<"\t"<< schlRunTm()<<"\t" << strMainCrs()<<"\t"
					<< ons()<<"\t" << offs() <<"\t"<<xc()<<"\t"<<yc()<<"\t"<<zc()<<"\t"<<endl ;
		};


		void show_tsdh(void)
		{ 
			cout << "Stop Id: " << id() << endl;
			cout << "Stop Name: " << stopCode() << endl;
			cout << "Direction: " << dirName() << endl;
		};

		void show_tsd(void)
		{ 
			cout << "Stop Id: " << id() << endl;
			cout << "Stop Order: " << stOrdr() << endl;
			cout << "Stope Code: " << stopCode() << endl;
			cout << "Route : " << routeName() << endl;
			cout << "Trip Id: " << tripKey() << endl;
			cout << "Trip Number: " << tripNumber() << endl;
			cout << "Trip Actual Arr. Time: " << acTimeArr() << endl;
			cout << "Trip Actual Dep. Time: " << acTimeDep() << endl;
			cout << "Scheduled Run Time : " << schlRunTm() << endl;
			cout << "Main Crossing Street  : " << strMainCrs() << endl;
			cout << "Boardings : " << ons() << endl;
			cout << "Alightings : " << offs() << endl;
			cout << "Xc: " << xc() << endl;
			cout << "Yc: " << yc() << endl;
			cout << "Zc: " << zc() << endl;
		};

			void sTrip::serialize(ofstream& pts);
			void sTrip::serializetext(ofstream& pts);
			void sTrip::serializetexthdr(ofstream& pts);
			void sTrip::deserialize(ifstream& pts);
			void sTrip::deserializetext(ifstream& pts);

		~sTrip(){
				//cout << "Deleting stop: " << *this << endl;
		}

	protected:
				long mid; // Stop id j
				string mschlName; // Schedule name
				string mrouteName; //Route Name
				string mdirName; // direction name
				string mtimePeriod; // Time Period
				string mtripSTime; // trip Start Time
				long mtripKey;  // Trip Key
				int mtripNumber; // Trip Number
				string mschlTime; // Schedule Time
				string macTimeArr; //Actual Arrival Time
				string macTimeDep; // Actual Departure Time
				int mschlRunTm; // Schedule Run Time Minutes
				int mstOrdr; // Stop order stop sequence
				long mstopCode; // Stop Code
				string mstopName; // Stop Name
				string mstrMainCrs; // Main crossing street
         		double mons; // Boardings 
         		double moffs; // Alightings 
				double mxc; // X Coordinate
         		double myc; // Y Coordinate
         		double mzc; // z Coordinate
	};

	void sTrip::serialize(ofstream& pts)
	{
	 pts.write(reinterpret_cast<char *>(&mid), sizeof(mid));
	 streamsize sizet=mschlName.size();// store Schedule Name's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mschlName.c_str(), sizet+1); // write Schedule Name including '\0' too
	 pts.write(reinterpret_cast<char *>(&mstOrdr), sizeof(mstOrdr));
	 sizet=mrouteName.size();// store Route name's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mrouteName.c_str(), sizet+1); // write Route name including '\0' too
	 sizet=mdirName.size();// store Direction name's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mdirName.c_str(), sizet+1); // write Direction name including '\0' too
	 sizet=mtimePeriod.size();// store Time Period's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mtimePeriod.c_str(), sizet+1); // write time period including '\0' too
	 sizet=mtripSTime.size();// store Trip Start Time's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mtripSTime.c_str(), sizet+1); // write Trip Start time including '\0' too
	 pts.write(reinterpret_cast<char *>(&mtripKey),sizeof(mtripKey));// Write Trip Key's length
	 pts.write(reinterpret_cast<char *>(&mtripNumber),sizeof(mtripNumber));// Write Trip Number's length
	 sizet=mschlTime.size();// store Scheduled Time's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mschlTime.c_str(), sizet+1); // write Scheduled Time including '\0' too
	 sizet=macTimeArr.size();// store Actual Arrival Time's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(macTimeArr.c_str(), sizet+1); // write Actual Arrival Time including '\0' too
	 sizet=macTimeDep.size();// store Actual Departure Time's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(macTimeDep.c_str(), sizet+1); // write Actual Departure Time including '\0' too
	 pts.write(reinterpret_cast<char *>(&mschlRunTm),sizeof(mschlRunTm));// store Schedule Run Time Minutes length
	 pts.write(reinterpret_cast<char *>(&mstOrdr),sizeof(mstOrdr));// store Stop Order length
	 pts.write(reinterpret_cast<char *>(&mstopCode),sizeof(mstopCode));// store Stop Code ID's length
	 sizet=mstopName.size();// store Trip ID's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mstopName.c_str(), sizet+1); // write Stop Name including '\0' too
	 sizet=mstrMainCrs.size();// store Main crossing street's length
	 pts.write(reinterpret_cast<char *>(&sizet),sizeof(sizet));
	 pts.write(mstrMainCrs.c_str(), sizet+1); // write Main crossing street's including '\0' too
	 pts.write(reinterpret_cast<char *>(&mons),sizeof(mons));
	 pts.write(reinterpret_cast<char *>(&moffs),sizeof(moffs));
	 pts.write(reinterpret_cast<char *>(&mxc),sizeof(mxc));
	 pts.write(reinterpret_cast<char *>(&myc),sizeof(myc));
	 pts.write(reinterpret_cast<char *>(&mzc),sizeof(mzc));

	}


	void sTrip::deserialize(ifstream& pts)
	{
	 int len=0;
	 char *p=0;
	 pts.read(reinterpret_cast<char *>(&mid), sizeof(mid));

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mschlName=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mrouteName=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mdirName=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mtimePeriod=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mtripSTime=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&mtripKey),sizeof(mtripKey)); // read stopcode 
	 pts.read(reinterpret_cast<char *>(&mtripNumber),sizeof(mtripNumber)); // read stopcode 

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mschlTime=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 macTimeArr=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 macTimeDep=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&mschlRunTm),sizeof(mschlRunTm)); // read schedule RunTm 
	 pts.read(reinterpret_cast<char *>(&mstOrdr),sizeof(mstOrdr)); // read stopcode 
	 pts.read(reinterpret_cast<char *>(&mstopCode),sizeof(mstopCode)); // read stopcode 

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mstopName=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&len), sizeof(len));
	 p=new char [len+1]; // allocate temp buffer for name
	 pts.read(p, len+1); // copy string to temp, including '\0'
	 mstrMainCrs=p; // copy temp to data member
	 delete[] p;

	 pts.read(reinterpret_cast<char *>(&mons),sizeof(mons));
	 pts.read(reinterpret_cast<char *>(&moffs),sizeof(moffs));
	 pts.read(reinterpret_cast<char *>(&mxc),sizeof(mxc));
	 pts.read(reinterpret_cast<char *>(&myc),sizeof(myc));
	 pts.read(reinterpret_cast<char *>(&mzc),sizeof(mzc));
	}


	sTrip::sTrip(long id1, string schlName1, string routeName1, string dirName1, string timePeriod1,
			string tripSTime1, long tripKey1, int tripNumber1,string schlTime1,	string acTimeArr1, 
			string acTimeDep1, int schlRunTm1 , string strMainCrs1,	int stOrdr1, long stopCode1, 
			double ons1, double offs1, double xc1,double yc1,double zc1 ) 
	{
			sTrip::mid = id1;
			sTrip::mrouteName = routeName1;
			sTrip::mtripKey= tripKey1; 
			sTrip::mtripNumber= tripNumber1; 
			sTrip::mtimePeriod= timePeriod1; 
			sTrip::mtripSTime= tripSTime1; 
			sTrip::mschlTime= schlTime1; 
			sTrip::macTimeArr= acTimeArr1; 
			sTrip::macTimeDep= acTimeDep1; 
			sTrip::mschlRunTm= schlRunTm1; 
			sTrip::mstrMainCrs= strMainCrs1; 
			sTrip::mstOrdr = stOrdr1;
			sTrip::mstopCode = stopCode1;
			sTrip::mons= ons1;
			sTrip::moffs= offs1;
			sTrip::mxc = xc1; 
			sTrip::myc = yc1; 
			sTrip::mzc = zc1; 
	}


	void sTrip::serializetext(ofstream& pts)
	{
	 pts << mid << "\t" << mstOrdr<< "\t"<< mrouteName << "\t"<< mstopCode
					<<"\t"<< mstopName <<"\t"<< mdirName<<"\t"<< mtripKey <<"\t" 
					<< mtripNumber <<  "\t"<< macTimeArr<<"\t" << macTimeDep<<"\t"
					<< mschlRunTm<<"\t"<< mstrMainCrs<<"\t" << mons<<"\t" << moffs 
					<<"\t"<<mxc<<"\t"<<myc<<"\t"<<mzc<<endl;
	}
	void sTrip::serializetexthdr(ofstream& pts)
	{
	 pts  << " mid " <<  "\t"  << " mstOrdr" <<  "\t" << " mrouteName " << "\t" << "mstopCode" << " \t" 
					<< " mstopName " << " \t"  << " mdirName" << " \t"  << " mtripKey " << " \t" 
					<< " mtripNumber " <<  "\t"  << " macTimeArr" << " \t"  << " macTimeDep" << " \t"
					<< " mschlRunTm" << " \t"  << " mstrMainCrs" << " \t"  << " mons" << " \t"  
					<< " moffs " << " \t"  << "mxc" << " \t"  << "myc" << " \t"<< "mzc"<< endl;
	}


	void sTrip::deserializetext(ifstream& pts)
	{
	 //pts >>mid>>mstOrdr>>mstopCode>>mtripKey>>mtripNumber>>mschlRunTm>>mons>>moffs>>mxc>>myc>>mzc>>endl;
	}


