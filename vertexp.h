#if !defined( _VERTEXP_H_ )
#define _VERTEXP_H_ 

// define the vertex object 

class vertexp {
public:
    
//    edgevect edgevert;

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
    return *this;
  }

	bool operator==(vertexp vx);

	vertexp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1,short toStop1, long index1,long lowlink1);

	vertexp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1, double tcost1,short toStop1, long index1,long lowlink1);

	void vxp(long pid1,long pidp1, short lbl1, double cost1, long fid1, long orig1,short toStop1, long index1,long lowlink1);

	void set_id (long pid1) {pid=pid1;}

	long get_id () const {return pid;}

	void set_idp (long pidp1) {pidp=pidp1;}

	long get_idp () const {return pidp;}  // predecessor

	void set_lbl (short lbl1) {lbl=lbl1;}

	short get_lbl () const {return lbl;}

	void set_toStop (short toStop1) {toStop=toStop1;}

	short get_toStop () const {return toStop;}

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

	void vertexp::serialize(ofstream& pfVx);
	void vertexp::deserialize(ifstream& pfVx);
	void vertexp::serializetext(ofstream& pVx);

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
//	         vertexp* pv; // predecessor vertex pointer in shortest path tree 
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

void vertexp::serializetext(ofstream& pVx)
{
 pVx <<pid<<"|"<<pidp<<"|"<<lbl<<"|"<<cost<<"|"<<tcost<<"|"<<fid<<"|"<<orig<<
	 "|"<<toStop<<"|"<<index<<"|"<<lowlink<<endl;
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