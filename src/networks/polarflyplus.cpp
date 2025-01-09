/*
  Copyright (c) 2007-2015, Trustees of The Leland Stanford Junior University
  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this list
  of conditions and the following disclaimer.
  Redistributions in binary form must reproduce the above copyright notice, this 
  list of conditions and the following disclaimer in the documentation and/or 
  other materials provided with the distribution.
  Neither the name of the Stanford University nor the names of its contributors 
  may be used to endorse or promote products derived from this software without 
  specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "booksim.hpp"
#include <vector>
#include <sstream>

#include "polarflyplus.hpp"
#include "random_utils.hpp"
#include "misc_utils.hpp"
#include "globals.hpp"

#define POLAR_LATENCY 
#define Hypercubeport 7
#define Polarflyport 8

int gP, gA, gG;

//Hypercube : Local (group)
//Polarfly  : Global

int polarfly_table[57][8]={{9,10,11,12,13,14,15,1}
,{16,17,18,12,19,20,21,2}
,{22,23,24,25,13,20,26,3}
,{27,28,29,25,19,14,30,4}
,{27,23,18,31,32,33,15,5}
,{22,28,11,34,35,33,21,6}
,{16,10,24,34,32,36,30,7}
,{9,17,29,31,35,36,26,8}
,{27,37,24,1,38,39,8,21}
,{40,23,41,1,35,19,7,42}
,{43,44,29,1,32,6,20,45}
,{46,25,34,47,1,2,48,31}
,{49,17,50,1,3,51,33,30}
,{22,52,18,1,4,36,53,54}
,{16,28,55,1,56,5,57,26}
,{22,37,29,2,56,51,7,15}
,{40,28,50,2,32,13,8,54}
,{49,52,24,2,35,5,14,45}
,{43,10,41,2,4,39,33,26}
,{27,44,11,2,3,36,57,42}
,{9,23,55,2,38,6,53,30}
,{16,37,41,31,3,6,14,54}
,{40,10,29,47,3,5,53,21}
,{9,28,18,48,3,39,7,45}
,{46,12,38,3,35,32,4,56}
,{43,52,55,34,3,19,8,15}
,{9,37,50,34,4,5,20,42}
,{40,17,24,48,4,6,57,15}
,{16,23,11,47,4,51,8,45}
,{49,44,55,31,4,13,7,21}
,{22,44,41,48,12,5,8,30}
,{43,17,11,25,38,5,7,54}
,{46,36,51,5,19,13,6,39}
,{27,52,50,47,12,6,7,26}
,{49,10,18,25,56,6,8,42}
,{46,33,14,53,7,8,57,20}
,{46,40,43,22,16,9,27,49}
,{9,52,41,25,32,51,57,21}
,{9,44,24,47,56,19,33,54}
,{46,37,44,10,23,28,17,52}
,{22,10,50,31,38,19,57,45}
,{27,10,55,48,35,51,20,54}
,{49,37,11,48,32,19,53,26}
,{40,52,11,31,56,39,20,30}
,{46,55,18,29,41,50,24,11}
,{40,37,55,25,12,36,33,45}
,{49,23,29,34,12,39,57,54}
,{43,28,24,31,12,51,53,42}
,{43,37,18,47,35,13,57,30}
,{27,17,41,34,56,13,53,45}
,{16,52,29,48,38,13,33,42}
,{40,44,18,34,38,51,14,26}
,{43,23,50,48,56,36,14,21}
,{22,17,55,47,32,39,14,42}
,{46,45,26,21,42,54,15,30}
,{16,44,50,25,35,39,53,15}
,{49,28,41,47,38,36,20,15}
};


//calculate the hop count between src and destination
int polarflyplusnew_hopcnt(int src, int dest) 
{
  int hopcnt;
  int dest_grp_ID, src_grp_ID; 
  int src_hopcnt, dest_hopcnt;
  int src_intm, dest_intm;
  int grp_output, dest_grp_output;
  int grp_output_RID;

  int _grp_num_routers= gA;
  int _grp_num_nodes =_grp_num_routers*gP;
  
  dest_grp_ID = int(dest/_grp_num_nodes);
  src_grp_ID = int(src / _grp_num_nodes);
  
  //source and dest are in the same group, either 0-1 hop
  if (dest_grp_ID == src_grp_ID) {
    if ((int)(dest / gP) == (int)(src /gP))
      hopcnt = 0;
    else
      hopcnt = 1;
    
  } else {
    //source and dest are in the same group
    //find the number of hops in the source group
    //find the number of hops in the dest group
    if (src_grp_ID > dest_grp_ID)  {
      grp_output = dest_grp_ID;
      dest_grp_output = src_grp_ID - 1;
    }
    else {
      grp_output = dest_grp_ID - 1;
      dest_grp_output = src_grp_ID;
    }
    grp_output_RID = ((int) (grp_output / (gP))) + src_grp_ID * _grp_num_routers;
    src_intm = grp_output_RID * gP;

    grp_output_RID = ((int) (dest_grp_output / (gP))) + dest_grp_ID * _grp_num_routers;
    dest_intm = grp_output_RID * gP;

    //hop count in source group
    if ((int)( src_intm / gP) == (int)( src / gP ) )
      src_hopcnt = 0;
    else
      src_hopcnt = 1; 

    //hop count in destination group
    if ((int)( dest_intm / gP) == (int)( dest / gP ) ){
      dest_hopcnt = 0;
    }else{
      dest_hopcnt = 1;
    }

    //tally
    hopcnt = src_hopcnt + 1 + dest_hopcnt;
  }

  return hopcnt;  
}


//packet output port based on the source, destination and current location
int polarflyplus_port(int rID, int source, int dest){
  int _grp_num_routers= gA;
  int _grp_num_nodes =_grp_num_routers*gP;

  int out_port = -1;
  int grp_ID = int(rID / _grp_num_routers); 
  int dest_grp_ID = int(dest/_grp_num_nodes);
  int grp_output=-1;
  int grp_RID=-1;
  
  //which router within this group the packet needs to go to
  if (dest_grp_ID == grp_ID) {
    grp_RID = int(dest / gP);
  } else {
    if (grp_ID > dest_grp_ID) {
      grp_output = dest_grp_ID;
    } else {
      grp_output = dest_grp_ID - 1;
    }
    grp_RID = int(grp_output /gP) + grp_ID * _grp_num_routers;
  }

  //At the last hop
  if (dest >= rID*gP && dest < (rID+1)*gP) {    
    out_port = dest%gP;
  } else if (grp_RID == rID) {
    //At the optical link
    out_port = gP + (gA-1) + grp_output %(gP);
  } else {
    //need to route within a group
    assert(grp_RID!=-1);

    if (rID < grp_RID){
      out_port = (grp_RID % _grp_num_routers) - 1 + gP;
    }else{
      out_port = (grp_RID % _grp_num_routers) + gP;
    }
  }  
 
  assert(out_port!=-1);
  return out_port;
}


PolarFlyplusNew::PolarFlyplusNew( const Configuration &config, const string & name ) :
  Network( config, name )
{

  _ComputeSize( config );
  _Alloc( );
  _BuildNet( config );
}

void PolarFlyplusNew::_ComputeSize( const Configuration &config )
{

  // LIMITATION
  // _n == # of dimensions within a group
  // _p == # of processors within a router
  _p = config.GetInt( "k" );// # of nodes in each switch=1
  _n = config.GetInt( "n" );
  _p=1;
  _n=1;
  assert(_n==1);

  _k = Polarflyport+Hypercubeport+1; // Polarfly + Hyoercube+  CPU  

  // FIX...
  gK = _p; gN = _n;

  //group : Hypercube
  _a = powi(2,Hypercubeport);
  _nodes = _a*powi(2, Polarflyport); 
  _num_of_switch = _nodes / _p;
  _channels = _num_of_switch * (Polarflyport + Hypercubeport); 
  _size = _num_of_switch;
  
  gG = _g;
  gP = _p;
  gA = _a;
  _grp_num_routers = gA;
  _grp_num_nodes =_grp_num_routers*gP;

}

void PolarFlyplusNew::_BuildNet( const Configuration &config )
{

  int _output=-1;
  int _input=-1;
  int _dim_ID=-1;
  int _num_ports_per_switch=Polarflyport+Hypercubeport+1;
  int c;

  ostringstream router_name;

  cout << " Polarflyplus " << endl;
  cout << " p = " << _p << " n = " << _n << endl;
  cout << " each switch - total radix =  "<< _k << endl;
  cout << " # of switches = "<<  _num_of_switch << endl;
  cout << " # of channels = "<<  _channels << endl;
  cout << " # of nodes ( size of network ) = " << _nodes << endl;
  cout << " # of groups (_g) = " << _g << endl;
  cout << " # of routers per group (_a) = " << _a << endl;

  for ( int node = 0; node < _num_of_switch; ++node ) {
    // ID of the group
    int grp_ID;
    grp_ID = (int) (node/_a);
    router_name << "router";
    
    router_name << "_" <<  node ;

    _routers[node] = Router::NewRouter( config, this, router_name.str( ), 
					node, _k, _k );
    _timed_modules.push_back(_routers[node]);

    router_name.str("");

    for ( int cnt = 0; cnt < _p; ++cnt ) {
      c = _p * node +  cnt;
      _routers[node]->AddInputChannel( _inject[c], _inject_cred[c] );

    }

    for ( int cnt = 0; cnt < _p; ++cnt ) {
      c = _p * node +  cnt;
      _routers[node]->AddOutputChannel( _eject[c], _eject_cred[c] );

    }

    //

    if (_n > 1 )  { cout << " ERROR: n>1 dimension NOT supported yet... " << endl; exit(-1); }

    //********************************************
    //   connect OUTPUT channels
    //********************************************
    // add hypercube output channel
    //
    //_chan[output] : {{src hypercube ports}, {src polarfly ports}}
    
    for ( int dim = 0; dim < _n; ++dim ) {
      for ( int cnt = 0; cnt < Hypercubeport; ++cnt ) {
	_output = (Polarflyport+Hypercubeport+1) * node + cnt;

	_routers[node]->AddOutputChannel( _chan[_output], _chan_cred[_output] );

#ifdef POLAR_LATENCY
	_chan[_output]->SetLatency(10);
	_chan_cred[_output]->SetLatency(10);
#endif
      }
    }

    //add polarfly output channel
    for ( int cnt = 0; cnt < Polarflyport; ++cnt ) {
      _output = (Polarflyport+Hypercubeport+1) * node + Hypercubeport + cnt;

      //_chan[_output].global = true;
      _routers[node]->AddOutputChannel( _chan[_output], _chan_cred[_output] );
#ifdef POLAR_LATENCY
      _chan[_output]->SetLatency(10);
      _chan_cred[_output]->SetLatency(10);
#endif
    }


    //********************************************
    //   connect INPUT channels
    //********************************************
    //_chan[input] :  {{dest hypercube ports}, {dest polarfly ports}}
    _num_ports_per_switch = Polarflyport + Hypercubeport;

    //hypercube calculation
    for ( int cnt = 0; cnt < Hypercubeport; ++cnt ) {
	_input = (node ^ powi(2, Hypercubeport)) * _num_ports_per_switch + cnt;
	_routers[node]->AddInputChannel( _chan[_input], _chan_cred[_input] );
      }

    //Polarfly  table refer
    int dest_polarport;
    for ( int cnt = 0; cnt < Polarflyport; ++cnt ) {
      for(int i = 0 ; i < 8; i++){
          if (polarfly_table[polarfly_table[grp_ID][cnt]][i]==grp_ID){
                dest_polarport=i;
	  }
      }
      _input = polarfly_table[_grpID][cnt]*powi(2,Hypercubeport)*_num_ports_per_switch + dest_polarport;
	  
      _routers[node]->AddInputChannel( _chan[_input], _chan_cred[_input] );
    }

  }

  cout<<"Done links"<<endl;
}


int PolarFlyplusNew::GetN( ) const
{
  return _n;
}

int PolarFlyplusNew::GetK( ) const
{
  return _k;
}

void PolarFlyplusNew::InsertRandomFaults( const Configuration &config )
{
 
}

double PolarFlyplusNew::Capacity( ) const
{
  return (double)_k / 8.0;
}

void PolarFlyplusNew::RegisterRoutingFunctions(){

  gRoutingFunctionMap["min_polarflynew"] = &min_polarflyplusnew;
  gRoutingFunctionMap["ugal_polarflynew"] = &ugal_polarflyplusnew;
}


void min_polarflyplusnew( const Router *r, const Flit *f, int in_channel, 
		       OutputSet *outputs, bool inject )
{
  outputs->Clear( );

  if(inject) {
    int inject_vc= RandomInt(gNumVCs-1);
    outputs->AddRange(-1, inject_vc, inject_vc);
    return;
  }

  int _grp_num_routers= gA;

  int dest  = f->dest;
  int rID =  r->GetID(); 

  int grp_ID = int(rID / _grp_num_routers); 
  int debug = f->watch;
  int out_port = -1;
  int out_vc = 0;
  int dest_grp_ID=-1;

  if ( in_channel < gP ) {
    out_vc = 0;
    f->ph = 0;
    if (dest_grp_ID == grp_ID) {
      f->ph = 1;
    }
  } 


  out_port = polarflyplus_port(rID, f->src, dest);

  //optical dateline
  if (out_port >=gP + (gA-1)) {
    f->ph = 1;
  }  
  
  out_vc = f->ph;
  if (debug)
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
	       << "	through output port : " << out_port 
	       << " out vc: " << out_vc << endl;
  outputs->AddRange( out_port, out_vc, out_vc );
}


//Basic adaptive routign algorithm for the polarfly
void ugal_polarflyplusnew( const Router *r, const Flit *f, int in_channel, 
			OutputSet *outputs, bool inject )
{
  //need 3 VCs for deadlock freedom

  assert(gNumVCs==3);
  outputs->Clear( );
  if(inject) {
    int inject_vc= RandomInt(gNumVCs-1);
    outputs->AddRange(-1, inject_vc, inject_vc);
    return;
  }
  
  //this constant biases the adaptive decision toward minimum routing
  //negative value woudl biases it towards nonminimum routing
  int adaptive_threshold = 30;

  int _grp_num_routers= gA;
  int _grp_num_nodes =_grp_num_routers*gP;
  int _network_size =  gA * gP * gG;

 
  int dest  = f->dest;
  int rID =  r->GetID(); 
  int grp_ID = (int) (rID / _grp_num_routers);
  int dest_grp_ID = int(dest/_grp_num_nodes);

  int debug = f->watch;
  int out_port = -1;
  int out_vc = 0;
  int min_queue_size;
  int nonmin_queue_size;
  int intm_grp_ID;
  int intm_rID;

  if(debug){
    cout<<"At router "<<rID<<endl;
  }
  int min_router_output, nonmin_router_output;
  
  //at the source router, make the adaptive routing decision
  if ( in_channel < gP )   {
    //dest are in the same group, only use minimum routing
    if (dest_grp_ID == grp_ID) {
      f->ph = 2;
    } else {
      //select a random node
      f->intm =RandomInt(_network_size - 1);
      intm_grp_ID = (int)(f->intm/_grp_num_nodes);
      if (debug){
	cout<<"Intermediate node "<<f->intm<<" grp id "<<intm_grp_ID<<endl;
      }
      
      //random intermediate are in the same group, use minimum routing
      if(grp_ID == intm_grp_ID){
	f->ph = 1;
      } else {
	//congestion metrics using queue length, obtained by GetUsedCredit()
	min_router_output = polarflyplus_port(rID, f->src, f->dest); 
      	min_queue_size = max(r->GetUsedCredit(min_router_output), 0) ; 

      
	nonmin_router_output = polarflyplus_port(rID, f->src, f->intm);
	nonmin_queue_size = max(r->GetUsedCredit(nonmin_router_output), 0);

	//congestion comparison, could use hopcnt instead of 1 and 2
	if ((1 * min_queue_size ) <= (2 * nonmin_queue_size)+adaptive_threshold ) {	  
	  if (debug)  cout << " MINIMAL routing " << endl;
	  f->ph = 1;
	} else {
	  f->ph = 0;
	}
      }
    }
  }

  //transition from nonminimal phase to minimal
  if(f->ph==0){
    intm_rID= (int)(f->intm/gP);
    if( rID == intm_rID){
      f->ph = 1;
    }
  }

  //port assignement based on the phase
  if(f->ph == 0){
    out_port = polarflyplus_port(rID, f->src, f->intm);
  } else if(f->ph == 1){
    out_port = polarflyplus_port(rID, f->src, f->dest);
  } else if(f->ph == 2){
    out_port = polarflyplus_port(rID, f->src, f->dest);
  } else {
    assert(false);
  }

  } else if(f->ph == 2){
    out_port = polarflyplus_port(rID, f->src, f->dest);
  } else {
    assert(false);
  }

  //optical dateline
  if (f->ph == 1 && out_port >=gP + (gA-1)) {
    f->ph = 2;
  }  

  //vc assignemnt based on phase
  out_vc = f->ph;

  outputs->AddRange( out_port, out_vc, out_vc );
}
