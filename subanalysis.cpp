/*
 // call the voronoi vertext routine
short onoff = 0;

fileName(vertinfname,".bin",outfilename,onoff);

originout.open(outfilename,  ios::out | ios::binary|ios::trunc);
writeBindblong(mmapVx0,originout);
if (originout.is_open()) {
originout.close();
}
originout.clear();  

 onoff = 1;
 fileName(vertinfname,".bin",outfilename,onoff);
originout.open(outfilename,  ios::out | ios::binary|ios::trunc);
writeBindblong(mmapVx1,originout);
 if (originout.is_open()) {
 originout.close();
 }
 originout.clear();  

// read the original edge data from the edge file
onoff = 0;
fileName(edgeinfname,".bin",outedgefilename,-1);
ifstream edgein(outedgefilename, ios::in | ios::binary);
mmaped0.clear();
readBinEdgeObjectData(mmaped0,ev,edgein);
if (edgein.is_open()) {
edgein.close();
}
edgein.clear();

 vxv2 (mmapVx0,mapvert0,mmaped0,mmapved,mmapVxStop,onoff);
// save the vertex voronoi result 
 fileName(vertinfname,"vorvx.bin",outvertfilename,onoff);

ofstream vertout0(outvertfilename, ios::out | ios::binary|ios::trunc);
writeBinObjectData(mapvert0,vx,vertout0);
 if (vertout0.is_open()) {
    vertout0.close();
 }
 vertout0.clear();  

ifstream vertin0(outvertfilename, ios::in | ios::binary);

// write the edge object data
    fileName(edgeinfname,"ed.bin",outedgefilename,onoff);
	outedgeobjectfile.clear();
	outedgeobjectfile.open(outedgefilename, ios::out | ios::binary|ios::trunc);
	writeBinObjectData(mmaped0,ev,outedgeobjectfile);

	if (outedgeobjectfile.is_open()) {
		outedgeobjectfile.close();
	}
	outedgeobjectfile.clear();

// do the ons vertex voronoi data i.e towards the stop using the ride time from the stop to the end
onoff = 1;
// read the original edge data from the edge file
fileName(edgeinfname,".bin",outedgefilename,-1);
edgein.open(outedgefilename, ios::in | ios::binary);
mmaped1.clear();
readBinEdgeObjectData(mmaped1,ev,edgein);
if (edgein.is_open()) {
edgein.close();
}
edgein.clear();

vxv2 (mmapVx1,mapvert1,mmaped1,mmapved,mmapVxStop,onoff);


fileName(vertinfname,"vorvx.bin",outvertfilename,onoff);
ofstream vertout1(outvertfilename, ios::out | ios::binary|ios::trunc);
writeBinObjectData(mapvert1,vx,vertout1);
if (vertout1.is_open()) {
vertout1.close();
}
vertout1.clear();

string strVertHdr= "pid,pidp,lbl,cost,tcost,fid,orig,toStop";
fileName(vertinfname,"vorvx.csv",outvertfilename,onoff);
vertout1.open(outvertfilename, ios::out | ios::trunc);
writeTextObjectData(mapvert1,vx,vertout1,strVertHdr);
if (vertout1.is_open()) {
vertout1.close();
}
vertout1.clear();
j=0;
il=0;

// export the processed vertices
mapverit = mapvert0.begin();
while(mapverit != mapvert0.end()) 
{ 
	  // get the vertex id & vertex object
     k =(*mapverit).first;
	vxp =  &(mapverit->second);
	if (vxp->get_lbl()==0) { il++;}
	if (vxp->get_orig()>0) { j++;}
	// show the vertex details 
//	vxp->show_vertof(outvertfile);
//    cout<<" Vertex Origin pointer "<< (long) &(*vxp)<<endl;
	mapverit++;
} // while map origin is present
outvertfile.close();
cout << " End Vertex Export. "<<j<<" of " <<mapvert.size()<< " veritices assigned "<<il<<" Were not assigned to any stop."<<endl;
evp->show_edgehdr(outfile);


onoff = 1;
fileName(edgeinfname,"ed.bin",outedgefilename,onoff);
ofstream edgeout(outedgefilename, ios::out | ios::binary|ios::trunc);
writeBinObjectData(mmaped1,ev,edgeout);
if (edgeout.is_open()) {
  edgeout.close();
}
edgeout.clear();

// read the edge data
onoff = 0;
mmaped0.clear();
fileName(edgeinfname,"ed.bin",outedgefilename,onoff);
edgein.open(outedgefilename, ios::in | ios::binary);
readBinEdgeObjectData(mmaped0,ev,edgein);
if (edgein.is_open()) {
edgein.close();
}
edgein.clear();

// edge voronoi process
onoff=0; // off - from stop to parcel

edv (mapvert0, mmaped0,mmapved,mmapEdgeStop,onoff);

fileName(edgeinfname,"vored.bin",outedgefilename,onoff);
edgeout.open(outedgefilename, ios::out | ios::binary|ios::trunc);
writeBinObjectData(mmaped0,ev,edgeout);

if (edgeout.is_open()) {
edgeout.close();
}
edgeout.clear();

 string strEdgeHdr = "id,eid,eoid,evp,lbl,toStop,dirn,";
	 strEdgeHdr.append("ecost,scost,tcost,slen,palong,frid,toid,");
	 strEdgeHdr.append("orig,tway,efc,frfc,tofc,orfc,ornm,enote");

fileName(edgeinfname,"vored.csv",outedgefilename,onoff);
edgeout.open(outedgefilename, ios::out |ios::trunc);
writeTextObjectData(mmaped0,ev,edgeout,strEdgeHdr);

if (edgeout.is_open()) {
edgeout.close();
}
edgeout.clear();


// read the edge data
onoff = 1;
mmaped1.clear();
fileName(edgeinfname,"ed.bin",outedgefilename,onoff);
edgein.open(outedgefilename, ios::in | ios::binary);
readBinEdgeObjectData(mmaped1,ev,edgein);
if (edgein.is_open()) {
edgein.close();
}
edgein.clear();

// process the ons edge voronoi routine
onoff=1;
edv (mapvert1, mmaped1,mmapved,mmapEdgeStop,onoff);

fileName(edgeinfname,"vored.bin",outedgefilename,onoff);
edgeout.open(outedgefilename, ios::out | ios::binary|ios::trunc);
writeBinObjectData(mmaped1,ev,edgeout);

if (edgeout.is_open()) {
edgeout.close();
}
edgeout.clear();

fileName(edgeinfname,"vored.csv",outedgefilename,onoff);
edgeout.open(outedgefilename, ios::out | ios::trunc);
writeTextObjectData(mmaped0,ev,edgeout,strEdgeHdr);

cout << " End of street partitioning..."<<endl<<" Start Assigning parcels to stops..."<<endl;

 //open parcel data output file - empty file if it already exists before writing new data
// ofstream outparcfile(strcat_s(strcat_s("C:\\NEU\\bustops\\Boston\\parcel",to_string(onoff)),".bin"), ios::out | ios::binary|ios::trunc);
// print header to parcel file
	parp->show_parcelhdr(outparcfile);
// Assign parcels to stops and update walk cost
	blnHist = true;

// read the parcel data from the parcel file
onoff = -1;
fileName(pariname,".bin",outfilename,onoff);
ifstream parcin(outfilename, ios::in | ios::binary);
mmaparced0.clear();
readBinParcObjectData(mmaparced0,par1,parcin);
if (parcin.is_open()) {
parcin.close();
}
parcin.clear();
mmaparStop.clear();
onoff=0;
parcstopwalk(mmapeidvid,mmaparced0,mapvert0, mmaped0,mmaparStop,onoff,blnHist);

//calculate the summation for the onoff=0 (offs)
ParcelOffAccum paoff;
deque<ParcelOffAccum> vpaoff;
deque<ParcelOffAccum>::iterator vpait;
 // calculate summary of parcel impact by tranist stop - (aggregate parcel walk time impact by transit stop)
//vpaoff = pagwalkstopoff (tsidmap,mmaparStop,luci,gcost,onoff,paoff,blnHist);

mmaparStop = pacalc (tsidmap,mmaparStop,luci,gcost,onoff,blnHist);
vpaoff = pagstop (tsidmap,mmaparStop,paoff);

// for historic run update the offs, hist. walk time offs 
updStopDmndWlkTm(vpaoff, tsidmap, onoff,blnHist);

// calculate the parcel ons and offs using the aggregate 
if (blnHist) {
	calcParcelDemand0 ( tsidmap, mmaparStop, onoff, blnHist);	
}

// aggregate the parcel cost by transit stops
vpaoff = pagwalkstop<maplngstop,mmaplngpar,ParcelOffAccum> (tsidmap,mmaparStop,paoff);

if (outparcfile.is_open()) {
	outparcfile.close();
}
outparcfile.clear();

fileName(pariname,"stopid.bin",outfilename,onoff);
outparcfile.open(outfilename, ios::out | ios::binary | ios::trunc);
writeBinIdplusObjectData(mmaparStop,par1,outparcfile);
if (outparcfile.is_open()) {
outparcfile.close();
}
outparcfile.clear();

//onoff = 0;
//fileName(pariname,"out.bin",outfilename,onoff);
parcin.clear();
parcin.open(outfilename, ios::in | ios::binary);
mmaparStop.clear();
readBinIdParcObjectData(mmaparStop,par1,parcin);
if (parcin.is_open()) {
parcin.close();
}
parcin.clear();

fileName(pariname,"out.bin",outfilename,onoff);
outparcfile.open(outfilename, ios::in | ios::binary | ios::trunc);
writeBinObjectData(mmaparStop,par1,outparcfile);
if (outparcfile.is_open()) {
outparcfile.close();
}
outparcfile.clear();

//vpaoff.clear();

//onoff = 0;
//fileName(pariname,"out.bin",outfilename,onoff);
parcin.open(outfilename, ios::in | ios::binary);
mmaparced1.clear();
readBinParcObjectData(mmaparced1,par1,parcin);
if (parcin.is_open()) {
parcin.close();
}
parcin.clear();

onoff = 1;
mmaparStop.clear();
parcstopwalk(mmapeidvid,mmaparced1,mapvert1, mmaped1,mmaparStop,onoff,blnHist);

fileName(pariname,"out.bin",outfilename,onoff);
outparcfile.open(outfilename, ios::out | ios::binary | ios::trunc);
writeBinObjectData(mmaparStop,par1,outparcfile);
if (outparcfile.is_open()) {
outparcfile.close();
}
outparcfile.clear();

fileName(pariname,"stopid.bin",outfilename,onoff);
outparcfile.open(outfilename, ios::out | ios::binary | ios::trunc);
writeBinIdplusObjectData(mmaparStop,par1,outparcfile);
if (outparcfile.is_open()) {
	outparcfile.close();
}
outparcfile.clear();

outparcfile.clear();
fileName(pariname,"pass.csv",outparfilename,onoff);
outparcfile.open(outparfilename, ios::out | ios::trunc);
par1.texthdr(outparcfile);
writeTextObjectData(mmaparStop,par1,outparcfile,"");
 if (outparcfile.is_open()) {
	outparcfile.close();
 }
outparcfile.clear();

// calculate summary of parcel impact by tranist stop - (aggregate parcel walk time impact by transit stop)
ParcelOnAccum paon;
deque<ParcelOnAccum> vpaon;
deque<ParcelOnAccum>::iterator vpaonit;
if (blnHist) {
	mmaparStop = pacalc (tsidmap,mmaparStop,luci,gcost,onoff,blnHist);
}
vpaon = pagstop (tsidmap,mmaparStop,paon);

updStopDmndWlkTm(vpaon, tsidmap, onoff,blnHist);

// calculate the parcel ons and offs using the aggregate 
if (blnHist) {
	calcParcelDemand0 ( tsidmap, mmaparStop, onoff, blnHist);	
}

// reaggregate the parcel data as a check against the historic demand per stop 
//vpaon = pagwalkstop<maplngstop,mmaplngpar,ParcelOnAccum> (tsidmap,mmaparStop,paon);
vpaon = pagstop (tsidmap,mmaparStop,paon);

//open parcel data output file - empty file if it already exists before writing new data
 onoff=-1;
 fileName(pariname,"out.bin",outparfilename,onoff);
 outparcfile.clear();
 outparcfile.open(outparfilename, ios::out | ios::binary|ios::trunc);
// write the parcel data to file
 writeBinObjectData(mmaparStop,par1,outparcfile);
// close parcel output file
 if (outparcfile.is_open()) {
	 outparcfile.close();
 }
 outparcfile.clear();

// write the parcel data to text file format file
onoff = 2;
fileName(pariname,"parcass.csv",outparfilename,onoff);
outparcfile.open(outparfilename, ios::out | ios::trunc);
par1.texthdr(outparcfile);
writeTextObjectData(mmaparStop,par1,outparcfile,"");
 if (outparcfile.is_open()) {
	outparcfile.close();
 }
outparcfile.clear();


//vpaon.clear();

 // calculate the undelayed run time
 tsidmit = tsidmap.begin();
 tsidmit2 = tsidmap.end();
    undelTm<maplngstop::iterator>(tsidmit, tsidmit2,gcost->get_unitontm(),
			gcost->get_unitofftm(), phway->get_hdway());

//    ofstream paout("C:\\NEU\\bustops\\Boston\\parcel.bin", ios::out | ios::binary);

// loop over the stops and calculate the riding and operating details
	tsidmit = tsidmap.begin();
	tsidmit2 = tsidmap.end();


	// for historic run update the offs, hist. walk time offs 
	onoff = 0;
	mmaparStop0.clear();
	SortOnOffIdParcObject(mmaparStop,mmaparStop0,onoff);
	vpaoff = pagstop (tsidmap,mmaparStop0,paoff);
	updStopDmndWlkTm(vpaoff, tsidmap, onoff,blnHist);

	// for historic run update the offs, hist. walk time ons 
	onoff = 1;
	mmaparStop1.clear();
	SortOnOffIdParcObject(mmaparStop,mmaparStop1,onoff);
	vpaon = pagstop (tsidmap,mmaparStop1,paon);
	updStopDmndWlkTm(vpaon, tsidmap, onoff,blnHist);

	ImpactCalc<maplngstop::iterator>(tsidmit,tsidmit2,*gcost,*phway);

//open parcel data output file - empty file if it already exists before writing new data
 onoff=-1;
 fileName(tstopiname,"out.csv",outfilename,onoff);

 ofstream stopout(outfilename, ios::out |ios::trunc);
 string txthdr =  "Stopid,Stopidp,Edgeid,palong,StOrdr,StopLbl,StopName,lbl,blnHist,blnInbd,blnExtr,blnIncl,blnElim,";
	txthdr.append("CumDist,CRdTm,undCRdTm,CRdTmC,HistOns,HistOffs,HistDepVol,Ons,Offs,DepVol,probStoph,probStop,");
	txthdr.append("depDelay,arrDelay,dwlDelay,totDelay,PVal,AVal,CRdTmE,hWkTmOns,hWkTmOffs,WkTmOns,WkTmOffs,");
	txthdr.append("WalkCost,RideCost,OperCost,TCost");

// write the parcel data to file
 writeTextObjectData(tsidmap,stop0,stopout,txthdr);
// close parcel output file
 if (stopout.is_open()) {
	 stopout.close();
 }
 stopout.clear();
*/
