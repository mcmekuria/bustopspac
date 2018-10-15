maplngstop&  stopspacing_analysisx (mmaplng& mmapeidid,mmaplngpar& mmaparced,maplngvx& mapvert0,
						   maplngvx& mapvert1, mmaplnged& mmaped, mmaplngpar& mmaparStop, 
						   mmapdblng& mmapSver, mmaplng& mmapVxStop, mmapdblng& mmapVx0, 
						   mmapdblng& mmapVx1, bool blnHist,char edgeinfname[],
						   char vertinfname[], char pariname[],char tstopiname[], globcost gc1, 
						   pdhdway phway, lucodes luci,int q);
void vxv1 (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,
		  mmaplng& mmapved,mmaplng& mmapEdgeStop,short onoff); //vertex voronoi
void vxv2 (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,
		  mmaplng& mmapved,mmaplng& mmapEdgeStop,short onoff); //vertex voronoi

void vxv1 (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,mmaplng& mmapved,mmaplng& mmapVxStop,short onoff) {
// vertex voronoi analysis
	// variables
	double ecost,vcost,tcost,xcost=0;
    long k=0, vid=0,eid=0,ip=0,pid=0,origin=0,v1=0,v2=0,vcnt=0,o1=0; // vid - vertex id, eid - edge id, pid - parcel id
    vertexp* vx1,*vx2;
    edgev* evp;
	bool forwd=true,rev=true,fixedv2=false,bndEdge1=false,bndEdge2=false;
/* start the search from the sorted vertex origin collection  
 and process the vertex voronoi partition */
mmapSverit = mmapSver.begin();
while(mmapSverit != mmapSver.end()) 
{  
	  // get the minimum of the origin list
      vcost =(*mmapSverit).first;
      vid =(*mmapSverit).second;
	  mapverit = mapvert.find(vid);
	  vx1 = &(mapverit->second);

//      vx1->show_vertof(outfile);
	  origin = vx1->get_orig();
	  //	cout<<"Vert scan "<< (long)&(*vx1)<<" id "<<vx1->get_id()<<" orig "<<vx1->get_orig()<<" cost "<<vx1->get_cost()<<endl;
      v1 = vx1->get_id();
		//vid = (*mumved_AIter).first;
		mumVSit = mmapVxStop.find(vid);
		if (mumVSit!=mmapVxStop.end())
		{
			o1=mumVSit->second;
			vx1->set_orig(o1);
		} 
	  // erase/remove the fixed vertex from the sorted vertex list (the heap) 
		mmapSver.erase(vcost);
	    outfile<<"Vertex Fixed  "<<v1<<" orig "<<origin<<" cost "<<vcost<<endl;

	  // search through the edges coming out of this vertex - vid
      pair<mumved_AIter, mumved_AIter> vedrange = mmapved.equal_range(vid);
	  //size_t j = distance(vedrange.first,vedrange.second);
	  for (mumavedit = vedrange.first; mumavedit!=vedrange.second;mumavedit++)
		  //while (mumavedit->first == vedrange.first->first )  //mmapved2.end()))
	  {
		//vid = (*mumved_AIter).first;
		eid = (*mumavedit).second;
		// find the edge cost attribute from the object
//        mmapedit = mmaped.find(eid);

		pair<edfmap_Iter, edfmap_Iter> edrange = mmaped.equal_range(eid);
		//size_t j = distance(edrange.first,edrange.second);
//  find the edge details of eid  

		for (mmapedi = edrange.first; mmapedi!=edrange.second;mmapedi++)
	    {
  			 evp = &((*mmapedi).second);
			 if(vx1->get_orig()<=0) 
			 {
				vx1->set_orig(o1);
			 }
			 
	    	 forwd=(evp->get_frid()==vid);  // if the vertex being scanned is the edge's 1st vertex (forward direction)
			 if (forwd) {
				mapverit = mapvert.find(evp->get_toid()); // get the end vertex
			 } else 
			 {
				mapverit = mapvert.find(evp->get_frid()); // get the end vertex
			 }
// fix vertex by setting label = -1.0
			  if (vx1->get_orig()>0) {
				vx1->set_lbl(-1);
			  } else
			  {
				vx1->set_orig(o1);
			  }
			 vx2 = &((*mapverit).second);
			 v2 = vx2->get_id();
             outfile << "EID " <<eid<< " i "<<v1<<" Lbl "<<vx1->get_lbl()<< " j "<<v2<<" Lbl "<<vx2->get_lbl()<<endl;
		  fixedv2=(vx2->get_lbl()==-1); // if v2 has already been fixed
		  if (fixedv2)
		  { // skip if this vertex is already labelled (edge can not be labeled here)
				bndEdge1 = (vx1->get_cost()<(vx2->get_cost()+evp->get_cost()));
				bndEdge2 = (vx2->get_cost()<(vx1->get_cost()+evp->get_cost()));

    		if(bndEdge2)
			{
				//This edge may need to be fixed in reverse or it may be a boundary.
				//Check if vx1's cost is less than the cost of vx2+edge cost
				if (bndEdge1) 
				{// this is a boundary edge that will not benefit from the two vertices.
					evp->set_dirn(-2);
					if (evp->get_lbl()!=-1)
					{	
						if (forwd) {
							if(evp->get_scost()>vx1->get_cost()) {evp->set_scost(vx1->get_cost());}
							if(evp->get_tcost()>vx2->get_cost()) {evp->set_tcost(vx2->get_cost());}
						} else {
							if(evp->get_scost()>vx1->get_cost()) {evp->set_scost(vx2->get_cost());}
							if(evp->get_tcost()>vx2->get_cost()) {evp->set_tcost(vx1->get_cost());}
						}
					}// edge is already fixed skip updating the costs
				}
				else {// reverse fix will be applied
					vx1->set_cost(vx2->get_cost()+evp->get_cost());
					vx1->set_orig(vx2->get_orig());
					evp->set_orig(vx2->get_orig());
					evp->set_scost(vx2->get_cost());
					evp->set_tcost(vx1->get_cost());
				}
			}
			else  
			{ // vx2 > vx1+ecost hence update vertex 2's cost and origin
					vx2->set_cost(vx1->get_cost()+evp->get_cost());
					vx2->set_orig(vx1->get_orig());
					evp->set_orig(vx1->get_orig());
					evp->set_scost(vx1->get_cost());
					evp->set_tcost(vx2->get_cost());
			}
			  // if the begining edge is the same as vid (the edge is in forward orientation)

		  }	  // if v2 is already fixed
		  else 
		  { // v2 is not fixed yet
			 vcnt++;
			// cout <<"Visting Edge "<<eid<< " , End Vertex = "<<v2<<" ,total visited nodes "<<vcnt<<endl;
			// get the end vertex id
			// v2 = evp->get_toid();
		        ecost = evp->get_cost();
		        tcost = vcost + ecost;
            // update the relevant attibutes of the edge using start vertex attributes
				if (vcost < evp->get_scost())
				{
				  evp->set_lbl(-1);
      			  evp->set_orig(o1);
                  evp->set_scost(vcost);
                  evp->set_tcost(tcost);
				  if (forwd) 
				  {evp->set_dirn(1);}
				  else
				  {evp->set_dirn(-1);}
				}
                //evp->show_edge(outedgefile);
    			outfile<<" eid "<<evp->get_id()<<" i "<<evp->get_frid()<<" j "<<
					evp->get_toid()<<" cost "<<evp->get_cost()<<" Scost "<<evp->get_scost()<<endl;
				if (vx2->get_cost()== inf) // if the successor vertex cost = inf then it is not visited yet 
				{
    				vx2->set_cost(tcost); // update cost
					vx2->set_orig(origin); // update origin
	    			vx2->set_idp(vx1->get_id()); //update pred id
		    		vx2->set_tcost(tcost); // update tcost
					vx2->show_vertof(outfile);
				    mmapSver.insert(dblng_Pair(tcost,v2)); // insert it into the
				}
				else
				{  if (tcost < vx2->get_cost())  // the new cost is less than old cost
				   {
    				vx2->set_cost(tcost); // update cost
	    			vx2->set_orig(origin);  // update origin
	    			vx2->set_idp(vx1->get_id());  // update pred id
		    		vx2->set_tcost(tcost);  // update tcost
					vx2->show_vertof(outfile);
				   } 
				   else
				   {
					   if (tcost == vx2->get_cost())
				      {
						// tie breaker 
						  vx2->set_fid(vx1->get_id());
					  }
				   }
				}// end vertex is visited i.e. if (cost=INF)
		  }  // if vertex is already labeled then skip it.
		  // end vertex vertex (for current edge) is being scanned/fixed , i.e. vid = evp->get_toid 
				  if(vx1->get_cost()>(vx2->get_cost()+evp->get_cost())) 
				  { // update vx1 cost and update edge origin
					vx1->set_cost(vx2->get_cost()+evp->get_cost());
					vx1->set_orig(vx2->get_orig());
					evp->set_orig(vx2->get_orig());
					evp->set_scost(vx2->get_cost());
					evp->set_tcost(vx1->get_cost());
				  }
			  outfile << "EID " <<eid<<" vid "<<vid<< " & i "<<evp->get_frid()<<" mismatch. E-Orig : "<<evp->get_orig()
				  <<" Lbl "<<vx1->get_lbl()<< " j "<<evp->get_toid()<<" E-Cost "<<evp->get_cost()<<" Es-Cost "
				  <<evp->get_scost()<<" E-Cost "<<evp->get_cost()<<" Et-Cost "<<vx1->get_tcost()<<endl;
		}; // while edge match is present
	  }; // while vertex/edge match is present	// advance to the next edge
        mmapSverit = mmapSver.begin();  // pick the next lowest cost vertex 
} // while map heap list is present

cout << "End Vertex Voronoi for onoff = "<<onoff<<endl;
mapverit = mapvert.begin(); // get the end vertex
o1=0;
		(mapverit->second).show_verthdr(cout);
while (mapverit!=mapvert.end()) {
	if ((mapverit->second).get_orig()<=0) {
		o1++;
		(mapverit->second).show_vertof(cout);
	}
	mapverit++;
}
cout << "No. of unassigned Vertices = "<<o1<<" for travel direction = "<<onoff<<endl;


}; //End vxv routine


void vxv2 (mmapdblng& mmapSver,maplngvx& mapvert,mmaplnged& mmaped,mmaplng& mmapved,mmaplng& mmapVxStop,short onoff) 
{
// vertex voronoi analysis
	// variables
	double ecost,vcost,tcost,xcost=0;
    long k=0, vid=0,eid=0,ip=0,pid=0,origin=0,v1=0,v2=0,vcnt=0,o1=0; // vid - vertex id, eid - edge id, pid - parcel id
    vertexp* vx1,*vx2;
    edgev* evp;
	bool forwd=true,rev=true,fixedv2=false,bndEdge1=false,bndEdge2=false;
/* start the search from the sorted vertex origin collection  
 and process the vertex voronoi partition */
mmapSverit = mmapSver.begin();
while(mmapSverit != mmapSver.end()) 
{  
	  // get the minimum of the origin list
      vcost =(*mmapSverit).first;
      vid =(*mmapSverit).second;
	  mapverit = mapvert.find(vid);
	  vx1 = &(mapverit->second);

//      vx1->show_vertof(outfile);
	  origin = vx1->get_orig();
	  //	cout<<"Vert scan "<< (long)&(*vx1)<<" id "<<vx1->get_id()<<" orig "<<vx1->get_orig()<<" cost "<<vx1->get_cost()<<endl;
      v1 = vx1->get_id();
		//vid = (*mumved_AIter).first;
		mumVSit = mmapVxStop.find(vid);
		// find the stop id related to this vertex
		if (mumVSit!=mmapVxStop.end())
		{
			o1=mumVSit->second;
			// label the vertex with the stop 
			vx1->set_orig(o1);
		  // erase/remove the fixed vertex from the sorted vertex list (the heap) 
			mmapSver.erase(vcost);
		    outfile<<"Vertex Fixed  "<<v1<<" orig "<<o1<<" cost "<<vcost<<endl;
		} 

	  // search through the edges coming out of this vertex - vid
      pair<mumved_AIter, mumved_AIter> vedrange = mmapved.equal_range(vid);
	  //size_t j = distance(vedrange.first,vedrange.second);
	  for (mumavedit = vedrange.first; mumavedit!=vedrange.second;mumavedit++)
	  {
		eid = (*mumavedit).second;
		// find the edge cost attribute from the object
		pair<edfmap_Iter, edfmap_Iter> edrange = mmaped.equal_range(eid);
		//size_t j = distance(edrange.first,edrange.second);
//  find the edge details of eid  
		for (mmapedi = edrange.first; mmapedi!=edrange.second;mmapedi++)
		{
  			 evp = &((*mmapedi).second);
			 if (vx1->get_id()==evp->get_frid() || vx1->get_id()==evp->get_toid()) {

		    	 forwd=(evp->get_frid()==vid);  // if the vertex being scanned is the edge's 1st vertex (forward direction)
				 if (forwd) {
					evp->set_dirn(1);
					mapverit = mapvert.find(evp->get_toid()); // get the end vertex
				 } else 
				 {
					evp->set_dirn(-1);
					mapverit = mapvert.find(evp->get_frid()); // get the end vertex
				 }
// fix vertex by setting label = -1.0
				  if (vx1->get_orig()<=0) {
					vx1->set_orig(o1);
				  }
				 vx1->set_lbl(-1);
				 vx2 = &((*mapverit).second);
				 v2 = vx2->get_id();
	             outfile << "EID " <<eid<< " i "<<v1<<" Lbl "<<vx1->get_lbl()<< " j "<<v2<<" Lbl "<<vx2->get_lbl()<<endl;
				fixedv2=(vx2->get_lbl()==-1); // if v2 has already been fixed
			  if (fixedv2)
			  { // skip if this vertex is already labelled (edge can not be labeled here)
				bndEdge1 = roundownx((vx1->get_cost()),3)<roundownx((vx2->get_cost()+evp->get_cost()),3);
				bndEdge2 = roundownx((vx2->get_cost()),3)<roundownx((vx1->get_cost()+evp->get_cost()),3);

	    		if(bndEdge2)
				{
				//This edge may need to be fixed in reverse or it may be a boundary.
				//Check if vx1's cost is less than the cost of vx2+edge cost
					if (bndEdge1) 
					{// this is a boundary edge that will not benefit from the two vertices.
						evp->set_dirn(-2);
						if (evp->get_lbl()!=-1)
						{	
							if (forwd) {
								if(evp->get_scost()>vx1->get_cost()) {evp->set_scost(vx1->get_cost());}
								if(evp->get_tcost()>vx2->get_cost()) {evp->set_tcost(vx2->get_cost());}
							} else {
								if(evp->get_scost()>vx1->get_cost()) {evp->set_scost(vx2->get_cost());}
								if(evp->get_tcost()>vx2->get_cost()) {evp->set_tcost(vx1->get_cost());}
							}
						}// edge is already fixed skip updating the costs
					}
					else {// reverse fix will be applied
						vx1->set_cost(vx2->get_cost()+evp->get_cost());
						vx1->set_orig(vx2->get_orig());
						evp->set_orig(vx2->get_orig());
						evp->set_scost(vx2->get_cost());
						evp->set_tcost(vx1->get_cost());
					}
				}
				else  
				{ // vx2 > vx1+ecost hence update vertex 2's cost and origin
					vx2->set_cost(vx1->get_cost()+evp->get_cost());
					vx2->set_orig(vx1->get_orig());
					evp->set_orig(vx1->get_orig());
					evp->set_scost(vx1->get_cost());
					evp->set_tcost(vx2->get_cost());
				}
				  // if the begining edge is the same as vid (the edge is in forward orientation)
			  }	  // if v2 is already fixed
		  else 
		  { // v2 is not fixed yet
			 vcnt++;
			// cout <<"Visting Edge "<<eid<< " , End Vertex = "<<v2<<" ,total visited nodes "<<vcnt<<endl;
			// get the end vertex id
			// v2 = evp->get_toid();
		        ecost = evp->get_cost();
		        tcost = vcost + ecost;
            // update the relevant attibutes of the edge using start vertex attributes
				if (vcost < evp->get_scost())
				{
				  evp->set_lbl(-1);
      			  evp->set_orig(o1);
                  evp->set_scost(vcost);
                  evp->set_tcost(tcost);
				  if (forwd) 
				  {evp->set_dirn(1);}
				  else
				  {evp->set_dirn(-1);}
				}
                //evp->show_edge(outedgefile);
    			outfile<<" eid "<<evp->get_id()<<" i "<<evp->get_frid()<<" j "<<
					evp->get_toid()<<" cost "<<evp->get_cost()<<" Scost "<<evp->get_scost()<<endl;
				if (vx2->get_cost()== inf) // if the successor vertex cost = inf then it is not visited yet 
				{
    				vx2->set_cost(tcost); // update cost
					vx2->set_orig(o1); // update origin
	    			vx2->set_idp(vx1->get_id()); //update pred id
		    		vx2->set_tcost(tcost); // update tcost
					vx2->show_vertof(outfile);
				    mmapSver.insert(dblng_Pair(tcost,v2)); // insert it into the
				}
				else
				{  if (tcost < vx2->get_cost())  // the new cost is less than old cost
				   {
    				vx2->set_cost(tcost); // update cost
	    			vx2->set_orig(origin);  // update origin
	    			vx2->set_idp(vx1->get_id());  // update pred id
		    		vx2->set_tcost(tcost);  // update tcost
					vx2->show_vertof(outfile);
				   } 
				   else
				   {
					   if (tcost == vx2->get_cost())
				      {
						// tie breaker 
						  vx2->set_fid(vx1->get_id());
					  }
				   }
				}// end vertex is visited i.e. if (cost=INF)

		  }  // if vertex is already labeled then skip it.
		  // end vertex vertex (for current edge) is being scanned/fixed , i.e. vid = evp->get_toid 
				  if(vx1->get_cost()>(vx2->get_cost()+evp->get_cost())) 
				  { // update vx1 cost and update edge origin
					vx1->set_cost(vx2->get_cost()+evp->get_cost());
					vx1->set_orig(vx2->get_orig());
					evp->set_orig(vx2->get_orig());
					evp->set_scost(vx2->get_cost());
					evp->set_tcost(vx1->get_cost());
				  }
			  outfile << "EID " <<eid<<" vid "<<vid<< " & i "<<evp->get_frid()<<" mismatch. E-Orig : "<<evp->get_orig()
				  <<" Lbl "<<vx1->get_lbl()<< " j "<<evp->get_toid()<<" E-Cost "<<evp->get_cost()<<" Es-Cost "
				  <<evp->get_scost()<<" E-Cost "<<evp->get_cost()<<" Et-Cost "<<vx1->get_tcost()<<endl;
			} /* if vertex vid either the beginning or end of the current edge object data 
			 this happens due to duplicate edge ids that are made up of several edge pieces.  */
		}; // while edge match is present
	  }; // while vertex/edge match is present	// advance to the next edge
        mmapSverit = mmapSver.begin();  // pick the next lowest cost vertex 
} // while map heap list is present

cout << "End Vertex Voronoi for onoff = "<<onoff<<endl;
mapverit = mapvert.begin(); // get the end vertex
o1=0;

if (mapverit != mapvert.end()) {
		(mapverit->second).show_verthdr(cout);

	while (mapverit!=mapvert.end()) {
		if ((mapverit->second).get_orig()<=0) {
		o1++;
		(mapverit->second).show_vertof(cout);
		}
		mapverit++;
	}
}
cout << "No. of unassigned Vertices = "<<o1<<" for travel direction = "<<onoff<<endl;


}; //End vxv2 routine

/*
//: C08:ConstInitialization.cpp
// Initializing const in classes
#include <iostream>
using namespace std;

class Fred {
  const int size;
public:
  Fred(int sz);
  void print();
};

Fred::Fred(int sz) : size(sz) {}
void Fred::print() { cout << size << endl; }

    //: C08:BuiltInTypeConstructors.cpp
    #include <iostream>
    using namespace std;

    class B {
      int i;
    public:
      B(int ii);
      void print();
    };

    B::B(int ii) : i(ii) {}
    void B::print() { cout << i << endl; }

    int main() {
      B a(1), b(2);
      float pi(3.14159);
      a.print(); b.print();
      cout << pi << endl;
    } ///:~


*/
/*
 tstop& readtstop1(string& rec1, tstop& pstop, map<int , string>::iterator maphdrit, char *seps ) // *seps = "," 
{
   static int rno = 0; //record number

	int i=0; int ib=10; int j=0;
//OBJECTID,ORDER,LABEL,STOPNAME,CUM_DIST,CUM_TIME,UNDCRTM,CUM_TIME_P,HISTORIC,INBOUND,HISTONS,HISTOFFS,	
//HISTDEPVOL,ONS,OFFS,DEPVOL,PROBSTOP,DEPDELAY,ARRDELAY,DELAY,PVAL,AVAL,INCLUDE,CUM_TIME2	
//WKTMONS,WKTMOFFS,EXTERNAL,WALKCOST,RIDECOST,OPERCOST,TCOST,EXTERNAL,ELIMINATED,STOPLOC,EDGEID,POSALONG
	unsigned long id; // Stop id
//			 long Stopidp; // predecessor Stop id 
			long Edgeid;
//	         tstop* stv; // predecessor vertex pointer in shortest path tree 
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
			double totDelay; 
			double PVal;
			double AVal;
			double CRdTmE; 
			double WkTmOns; 
			double WkTmOffs;
			double WalkCost;
			double RideCost;
			double OperCost;
			double TCost;
			double posalong;


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
	i = mapfldsit->first;
	fldval = mapfldsit->second;
//	val1 = fldval.c_str();
	fldhdr = maphdrit->second;

	f1 = (stringUpper<string>(fldhdr)); 
	if (f1 == "OBJECTID") // tstop id
	{  id = fromString<long>(fldval); //strtol(val1,&stop1,ib);
	  if (id>0) {
		  pstop.set_id(id);
	  }
	  else
	  {
		  return pstop;
	  }
	}
	else if (f1 == "ORDER") // Order
	{
		  StOrdr = fromString<int>(fldval); //strtol(val1,&stop1,ib);
		  pstop.set_StOrdr(StOrdr);
	}
	else if (f1 == "LABEL") // stop label
	{
		StopLbl = fldval;
		pstop.set_StopLbl(StopLbl);
	}
	else if (f1 == "STOPNAME")  //stop name
	{
		StopName = fldval;
		pstop.set_StopName(StopName);
	  }
	else if (f1 == "CUM_DIST") // cum. distance
	{
		CumDist = fromString<double>(fldval); // strtod(val1,&stop1);
		pstop.set_CumDist(CumDist);
	} 
	else if (f1 == "CUM_TIME") // cum. ride time
	{ 
		CRdTm = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_CRdTm(CRdTm);
	}
	else if (f1 == "HISTORIC")
	{
		blnHist = fromString<short>(fldval); // strtol(val1,&stop1,ib);
		pstop.set_blnHist(blnHist);
	} 
	else if (f1 == "INBOUND") // Inbound
	{
		blnInbd = fromString<short>(fldval); //strtol(val1,&stop1,ib);
		pstop.set_blnInbd(blnInbd);
	}
	else if (f1 == "INCLUDE") // Include
	{
		blnIncl = fromString<short>(fldval); //strtol(val1,&stop1,ib);
		pstop.set_blnIncl(blnIncl);
	}
	else if (f1 == "HISTONS") // HistOns
	{
		HistOns = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_HistOns(HistOns);} 
	else if (f1 == "HISTOFFS") // HistOffs
	{	HistOffs = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_HistOffs(HistOffs);
	} 
	else if (f1 == "HISTDEPVOL") // HistDepVol
	{	
		HistDepVol = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_HistDepVol(HistDepVol);
	} 
	else if (f1 == "ONS") // Ons
	{
		Ons = fromString<double>(fldval); 
		pstop.set_Ons(Ons);} 
	else if (f1 == "OFFS") // Offs
	{	
		Offs = fromString<double>(fldval);
		pstop.set_Offs(Offs);
	} 
	else if (f1 == "DEPVOL") // DepVol
	{	
		DepVol = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_DepVol(DepVol);
	} 
	else if (f1 == "PROBSTOP") // DepVol
	{	
		probStop = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_probStop(probStop);
	}
	else if (f1 == "DEPDELAY") // Departure Delay
	{	
		depDelay = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_depDelay(depDelay);
	}
	else if (f1 == "ARRDELAY") // Arrival Delay
	{	
		arrDelay = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_arrDelay(arrDelay);
	}
	else if (f1 == "DELAY") // total Delay
	{	
		totDelay = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_totDelay(totDelay);
	}
	else if (f1 == "AVAL") // Attraction sum 
	{	
		AVal = fromString<double>(fldval); 
		pstop.set_AVal(AVal);
	}
	else if (f1 == "PVAL") // Production sum 
	{	
		PVal = fromString<double>(fldval); 
		pstop.set_PVal(PVal);
	}
	else if (f1 == "CUM_TIME2") // Ride time from stop to the end of the line (for ons) WkTmOns 
	{	
		CRdTmE = fromString<double>(fldval); 
		pstop.set_CRdTmE(CRdTmE);
	}
	else if (f1 == "CUM_TIME_P") // cum. ride time
	{ 
		CRdTmC = fromString<double>(fldval); //strtod(val1,&stop1);
		pstop.set_CRdTmC(CRdTmC);
	}
	else if (f1 == "UNDCRTM") // Ride time from stop to the end of the line (for ons) WkTmOns 
	{	
		undCRdTm = fromString<double>(fldval); 
		pstop.set_undCRdTm(undCRdTm);
	}
	else if (f1 == "WKTMONS") // Walk time summary for ons
	{	
		WkTmOns = fromString<double>(fldval); 
		pstop.set_WkTmOns(WkTmOns);
	}
	else if (f1 == "WKTMOFFS") // Walk time summary for ons
	{	
		WkTmOffs = fromString<double>(fldval); 
		pstop.set_WkTmOffs(WkTmOffs);
	}
	else if (f1 == "WALKCOST") // Walk Cost
	{	
		WalkCost = fromString<double>(fldval); 
		pstop.set_WalkCost(WalkCost);
	}
	else if (f1 == "RIDECOST") // Ride Cost
	{	
		RideCost = fromString<double>(fldval); 
		pstop.set_RideCost(RideCost);
	}
	else if (f1 == "OPERCOST") // Operating Cost
	{	
		OperCost = fromString<double>(fldval); 
		pstop.set_OperCost(OperCost);
	}
	else if (f1 == "TCOST") // Total Cost
	{	
		TCost = fromString<double>(fldval); 
		pstop.set_TCost(TCost);
	}
	else if (f1 == "EXTERNAL") // External stop
	{	
		blnExtr = fromString<short>(fldval); 
		pstop.set_blnExtr(blnExtr);
	} 
	else if (f1 == "ELIMINATED") // ELIMINATED stop
	{	
		blnElim = fromString<short>(fldval); 
		pstop.set_blnElim(blnElim);
	} 
	else if (f1 == "EDGEID") // Edgeid for stop
	{	
		Edgeid = fromString<long>(fldval); 
		pstop.set_Edgeid(Edgeid);
	} 
	else if (f1 == "POSALONG" || f1 == "PALONG") // Position along edge 
	{	
		posalong = fromString<double>(fldval); 
		pstop.set_posalong(posalong);
	} 

   // get the next fld
	mapfldsit++;
	maphdrit++;
}
   if (pstop.get_id()>0)
   {
	cout<<"Stop id "<<id<<", Name : "<<StopName<<", Cum Rd Tm "<<CRdTm<<", Historic "<<blnHist<<", Hist Ons "<<HistOns<<", Hist Offs "<<HistOffs<<endl;
   }

   return pstop;
}
*/

/*
	for (nit1=objrange.first; nit1!=objrange.second;nit1++)
	{
		 // insert the begining and ending vertices of this edge into the heap set
		eoid = nit1->second;			
		mit = m1.find(eoid);
		if (mit !=m1.end()) 
		{
		// insert the begining and ending vertices of this edge into the heap set
				o1= &(mit->second);
		 // get the vertices from the vert list to append it to the current edge
				v1 = o1->get_frid();
			    v2 = o1->get_toid();
				vcost = p1.get_CRdTm();
				xcost = p1.get_CRdTmE();
			// get the stop from the stop list and update the edge and vertex properties
				o1->set_orig(id);
				o1->set_lbl(-1);
				o1->set_dirn(1);
				mmapStopEdge.insert(lng_Pair(id,eoid));
				mmapEdgeStop.insert(lng_Pair(eid,eoid));
			// check to see if the direction of the edge is forward or reverse using the start cost (scost) & end cost (tcost)
			
				posalong = p1.get_posalong() ;
				o1->set_palong(posalong);
				vxmap_Iter = mapvert.find(v1);
				xdist = posalong * o1->get_cost();
			if ( vxmap_Iter != mapvert.end() )
			{ // get the vertex from the vertex list and update the vertex properties
				vx1 =  &(vxmap_Iter->second);
				o1->set_scost(vcost+xdist);
				o1->set_tcost((vcost+o1->get_cost()-xdist));
				vx1->set_origin_vals(-1,o1->get_lbl(),o1->get_orig(),o1->get_scost(),o1->get_scost());
				mmapSver.insert(dblng_Pair(vcost,v1));
				mmapVx0.insert(dblng_Pair((vcost+xdist),v1));
				mmapVx1.insert(dblng_Pair((xcost+xdist),v1));
				mmapVx0R.insert(lngdbl_Pair(v1,(vcost+xdist)));
				mmapVx1R.insert(lngdbl_Pair(v1,(xcost+xdist)));
				mmapStopVx.insert(lng_Pair(id,v1));
				mmapVxStop.insert(lng_Pair(v1,id));

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
					xdist = (1-posalong)*o1->get_cost();
					mmapVx0.insert(dblng_Pair((vcost+xdist),v2));
					mmapVx1.insert(dblng_Pair((xcost+xdist),v2));
					mmapVx0R.insert(lngdbl_Pair(v2,(vcost+xdist)));
					mmapVx1R.insert(lngdbl_Pair(v2,(xcost+xdist)));
					mmapStopVx.insert(lng_Pair(id,v2));
					mmapVxStop.insert(lng_Pair(v2,id));
	   				vx2->set_origin_vals(-1,o1->get_lbl(),o1->get_orig(),o1->get_tcost(),o1->get_tcost());
					mmapSver.insert(dblng_Pair(vx2->get_cost(),v2));
				} else		
				{
					cout << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<<endl;
					f1 << "The vertex set doesn't have an element "
					<< "with a key of " << v2 << " for Edge "<<eid<< endl;
				}  // 4i
			}// if the edge is found
			else		
			{
				cout << "There is no edge with a key of " << eid<< 
					" for the stop " << id << " in this set. "<<endl;
				f1 << "There is no edge with a key of " << eid<< 
					" for the stop " << id << " in this set. "<<endl;
			}  // 4i if the edge id is in the edge-object map
		} // if edge object id is found in the edgeid-objectid map
	}
*/

/*	
	edgev(const edgev& ev) {
        cout << "Copy Edge ev[" << id << "]" << endl;
       ++copycons;
	  }
	  edgev& operator=(const edgev& ev)  {
        cout << "Assign Edge (" << id << ")=[" << ev.id << "]" << endl;
        id = ev.id;
        ++assign;
        return *this;
	  }

  friend bool operator<(const edgev& evb, const edgev& evc ) {
    return evb.id < evc.id;
  }
  friend bool operator==(const edgev& evb,const edgev& evc) {
    return evb.id == evc.id;
  }
  virtual ~edgev() {
    std::cout << "~[" << id << "]" << endl;
    ++destroy;
  }

  friend ostream& operator<<(ostream& os,const edgev& evc) {
    return os << evc.id;
  }
//  friend class edrpt;

*/

/*

	outparcfile.open(outfilename);
	while (!outparcfile) 
	{ //i1
		cout << "Parcel Output file could not be created! "<<endl;
		cout << "Please re-enter file name again : "<<endl; 
		   cin>>outfilename;
    	outparcfile.clear();
    	outparcfile.open(outfilename);
	    if (++j>3) 
		{ //2i 
		   cout << "More than three trials!"<<endl;
		   cout << "Press enter to exit!"<<endl;
		   cin>>outfilename;
		   exit (0);
		} //2i
	} //1i

	if (outparcfile.is_open()) {
		outparcfile.close();
	}
	outparcfile.clear();

*/

/*
	/*		
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
// write out the vertex results
	fileName(vertinfname,".trv",outfilename);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

	outfile.open(outfilename);
	str1 = "pid|pidp|lbl|cost|tcost|fid|orig|toStop|Index|LowIndex";

	writeTextObjectData(mapvert,pVx1,i,outfile,str1);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

// aggregate the vertex results
// first remap the key of the vertices by stop
	mmapvert1 = remapVertId2Orig(mapvert,mmapvert1,*pVx1,ip);	

	VertAccum vertAcc;
	deque<VertAccum> vertSumAcc;
	vertSumAcc = aggbystop (tsidmap,mmapvert1,*pstop,*pVx1,vertAcc);
	str1 = "Summary of Vertex attributes by id \nId\tStop\tcost\tTotCost\tPa.Cnt\tindex\tLowLink";
	cout << str1<<endl;
	writeObjData(vertSumAcc,vertAcc,str1,cout);

	// write out the vertex results
	fileName(vertinfname,".inx",outfilename);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

	outfile.open(outfilename);
	str1 = "pid|pidp|lbl|cost|tcost|fid|orig|toStop|Index|LowIndex";

	writeTextObjectData(mapvert,pVx1,i,outfile,str1);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

	fileName(vertinfname,".iny",outfilename);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

	outfile.open(outfilename);
	str1 = "pid|pidp|lbl|cost|tcost|fid|orig|toStop|Index|LowIndex";

	writeTextObjectData(mapvert2,pVx1,i,outfile,str1);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
	mapvert2.clear();

*/


//	mumavedit = mmapV1V2.begin();
// since the edge mapping mmapV1V2 & mapvert are sorted by id the lowest 
//	vertext id will be used to connect all other vertices to it.
//	vxmap_Iter = mapvert.find(mumavedit->first);

/*
	vxmap_Iter = mapvert.begin();
	if (vxmap_Iter!=mapvert.end())
	{
		pVx1 = &(vxmap_Iter->second);
		v1 = pVx1->get_id();
		vxmap_Iter++;
		while (vxmap_Iter!=mapvert.end())
		{
			pVx2 = &(vxmap_Iter->second);
			v2 = pVx2->get_id();
			mmapV1V2.insert (lng_Pair(v1,v2));
			vxmap_Iter++;
		}
	}
*/

/*
		vxmap_Iter = mapvert.begin();
		pVx1 = &(vxmap_Iter->second);
// run the Strongly Connected Components Tarjan's routine
//	it is not necessary to first create a link from the first point to all the others so that there will be continuity
		tarjanr(pVx1);
	mmapvert.clear();
	mmapvert = remapVertId2Index(mapvert,mmapvert,*pVx1,ip);	
	VxVorAccum VxVAcc;
	deque<VxVorAccum> VxVorSumAcc;
	ip=inf;
	VxVorSumAcc = vertaggbyMaxKey (ip,j1,mmapvert,*pVx1,VxVAcc);

// write out the vertex results
	fileName(vertinfname,".trj",outfilename);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();

	outfile.open(outfilename);
	str1 = "pid|pidp|lbl|cost|tcost|fid|orig|toStop|Index|LowIndex";

	writeTextObjectData(mmapvert,pVx1,i,outfile,str1);
	if (outfile.is_open()) {
		outfile.close();
	}
	outfile.clear();
*/

// run a Vertex Voronoi


//   vxvtarx (mmapV1V2EgId,mmapSver,mmapVxStop,mapvert,mapvert2,maped,onoff);

// update the costs for the vertices that are assigned to stops using the stop table
		//updateVert(mapvert,tsidmap,mmapSver);
