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
#define Hypercubeport 1
#define Polarflyport 3

int gP_polar, gA_polar, gG_polar;

//Hypercube : Local (group)
//Polarfly  : Global
vector<vector<string>> dbg;
// node0-p0, node1-p0, hypercube
/*
int polarfly_table[57][8]=
{
{8,9,10,11,12,13,14,0},
{15,16,17,11,18,19,20,1},
{21,22,23,24,12,19,25,2},
{26,27,28,24,18,13,29,3},
{26,22,17,30,31,32,14,4},
{21,27,10,33,34,32,20,5},
{15,9,23,33,31,35,29,6},
{8,16,28,30,34,35,25,7},
{26,36,23,0,37,38,7,20},
{39,22,40,0,34,18,6,41},
{42,43,28,0,31,5,19,44},
{45,24,33,46,0,1,47,30},
{48,16,49,0,2,50,32,29},
{21,51,17,0,3,35,52,53},
{15,27,54,0,55,4,56,25},
{21,36,28,1,55,50,6,14},
{39,27,49,1,31,12,7,53},
{48,51,23,1,34,4,13,44},
{42,9,40,1,3,38,32,25},
{26,43,10,1,2,35,56,41},
{8,22,54,1,37,5,52,29},
{15,36,40,30,2,5,13,53},
{39,9,28,46,2,4,52,20},
{8,27,17,47,2,38,6,44},
{45,11,37,2,34,31,3,55},
{42,51,54,33,2,18,7,14},
{8,36,49,33,3,4,19,41},
{39,16,23,47,3,5,56,14},
{15,22,10,46,3,50,7,44},
{48,43,54,30,3,12,6,20},
{21,43,40,47,11,4,7,29},
{42,16,10,24,37,4,6,53},
{45,35,50,4,18,12,5,38},
{26,51,49,46,11,5,6,25},
{48,9,17,24,55,5,7,41},
{45,32,13,52,6,7,56,19},
{45,39,42,21,15,8,26,48},
{8,51,40,24,31,50,56,20},
{8,43,23,46,55,18,32,53},
{45,36,43,9,22,27,16,51},
{21,9,49,30,37,18,56,44},
{26,9,54,47,34,50,19,53},
{48,36,10,47,31,18,52,25},
{39,51,10,30,55,38,19,29},
{45,54,17,28,40,49,23,10},
{39,36,54,24,11,35,32,44},
{48,22,28,33,11,38,56,53},
{42,27,23,30,11,50,52,41},
{42,36,17,46,34,12,56,29},
{26,16,40,33,55,12,52,44},
{15,51,28,47,37,12,32,41},
{39,43,17,33,37,50,13,25},
{42,22,49,47,55,35,13,20},
{21,16,54,46,31,38,13,41},
{45,44,25,20,41,53,14,29},
{15,43,49,24,34,38,52,14},
{48,27,40,46,37,35,19,14}
};
*/
/*
int polarfly_table[13][4]=
{
{4,5,6,0},
{7,8,6,1},
{7,5,9,2},
{4,8,9,3},
{7,10,0,3},
{11,8,0,2},
{12,9,0,1},
{4,10,1,2},
{11,5,1,3},
{12,6,2,3},
{12,11,7,4},
{12,10,8,5},
{11,10,9,6}
};
*/

int polarfly_table[7][3]=
{
{3,4,0},
{5,4,1},
{6,4,2},
{6,5,0},
{0,1,2},
{6,3,1},
{5,3,2}
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

  int _grp_num_routers= gA_polar;
  int _grp_num_nodes =_grp_num_routers*gP_polar;
  
  dest_grp_ID = int(dest/_grp_num_nodes);
  src_grp_ID = int(src / _grp_num_nodes);
  
  //source and dest are in the same group, either 0-1 hop
  if (dest_grp_ID == src_grp_ID) {
    if ((int)(dest / gP_polar) == (int)(src /gP_polar))
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
    grp_output_RID = ((int) (grp_output / (gP_polar))) + src_grp_ID * _grp_num_routers;
    src_intm = grp_output_RID * gP_polar;

    grp_output_RID = ((int) (dest_grp_output / (gP_polar))) + dest_grp_ID * _grp_num_routers;
    dest_intm = grp_output_RID * gP_polar;

    //hop count in source group
    if ((int)( src_intm / gP_polar) == (int)( src / gP_polar ) )
      src_hopcnt = 0;
    else
      src_hopcnt = 1; 

    //hop count in destination group
    if ((int)( dest_intm / gP_polar) == (int)( dest / gP_polar ) ){
      dest_hopcnt = 0;
    }else{
      dest_hopcnt = 1;
    }

    //tally
    hopcnt = src_hopcnt + 1 + dest_hopcnt;
  }

  return hopcnt;  
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
  //_p = config.GetInt( "k" );// # of nodes in each switch=1
  //_n = config.GetInt( "n" );
  _p=1;
  _n=1;
  assert(_n==1);

  _k = Polarflyport + Hypercubeport + 1; // Polarfly + Hyoercube+  CPU  

  // FIX...
  gK = _p; gN = _n;

  //group : Hypercube
  _a = powi(2,Hypercubeport);
  //_g = 57; //F7 Polarfly
  //_g = 13; //F3 Polarfly
  _g = 7; //F2 Polarfly	    
  _nodes   = _a * _p * _g;
 
  _num_of_switch = _nodes / _p;
  _channels = _num_of_switch * (Polarflyport + Hypercubeport); 
  _size = _num_of_switch;
  
  gG_polar = _g;
  gP_polar = _p;
  gA_polar = _a;
  _grp_num_routers = gA_polar;
  _grp_num_nodes =_grp_num_routers*gP_polar;

}

void PolarFlyplusNew::_BuildNet( const Configuration &config )
{

  int _output=-1;
  int _input=-1;
  //int _dim_ID=-1;
  //int _num_ports_per_switch=Polarflyport+Hypercubeport+1;
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
  int ch_count=0;
  for ( int node = 0; node < _num_of_switch; ++node ) {
    // ID of the group
    int grp_ID;
    grp_ID = (int) (node/_a);
    router_name << "router";
    
    //router_name << "_" <<  node ;
    for(int i = 0 ; i < Hypercubeport; i++){
     router_name << "_" << ((node >> i)%2) ;
    }
    router_name << "_" << (node >> Hypercubeport);
    //cout << router_name.str( ) << endl;
    _routers[node] = Router::NewRouter( config, this, router_name.str(), 
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

    if (_n > 1 )  { cout << " ERROR: n>1 dimension NOT supported yet... " << endl; exit(-1); }

    //********************************************
    //   connect OUTPUT channels
    //********************************************
    // add hypercube output channel
    //
    //_chan[output] : {{src hypercube ports}, {src polarfly ports}}
    
      for ( int cnt = 0; cnt < Hypercubeport; ++cnt ) {
	_output = (Polarflyport+Hypercubeport) * node + cnt;
        dbg.push_back({ to_string(_output), "Hypercube","node"+to_string(node)+"-port"+to_string(cnt) });
	_routers[node]->AddOutputChannel( _chan[_output], _chan_cred[_output] );

#ifdef POLAR_LATENCY
	_chan[_output]->SetLatency(10);
	_chan_cred[_output]->SetLatency(10);
#endif
      }
    //add polarfly output channel
    for ( int cnt = 0; cnt < Polarflyport; ++cnt ) {
      _output = (Polarflyport+Hypercubeport) * node + Hypercubeport + cnt;
      //_chan[_output].global = true;
      //if(grp_ID < Polarflyport  && cnt==Polarflyport-1) continue; //red group : no self-connection
      dbg.push_back({ to_string(_output), "Polarfly","node"+to_string(node)+"-port"+to_string(cnt+Hypercubeport) });
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
    //_num_ports_per_switch = Polarflyport + Hypercubeport;
   
    //hypercube calculation
    for ( int cnt = 0; cnt < Hypercubeport; ++cnt ) {
	_input = (node ^ (1<<cnt)) * (Polarflyport+Hypercubeport) + cnt;
       dbg[ch_count].push_back (to_string(_input)+" node"+to_string(node ^ (1<<cnt))+"-port"+to_string(cnt));
       ch_count++;
	_routers[node]->AddInputChannel( _chan[_input], _chan_cred[_input] );
      }
    //Polarfly  table refer
    int dest_polarport;
    for ( int cnt = 0; cnt < Polarflyport; ++cnt ) {
      for(int i = 0 ; i < 8; i++){
          if (polarfly_table[polarfly_table[grp_ID][cnt]][i]==grp_ID){
                dest_polarport=i+Hypercubeport;
		break;
	  }
      }
      int hyperadd = node%(powi(2,Hypercubeport));
      int dest_node_add = polarfly_table[grp_ID][cnt]*powi(2,Hypercubeport)+hyperadd;
      //cout << hyperadd << " " << polarfly_table[grp_ID][cnt] << endl; 
      _input = dest_node_add * (Polarflyport+Hypercubeport) + dest_polarport;
      //if(grp_ID < Polarflyport && cnt== Polarflyport-1) continue; //red group : no self-connection
       dbg[ch_count].push_back (to_string(_input)+" node"+to_string(dest_node_add)+"-port"+to_string(dest_polarport));
       ch_count++;
      _routers[node]->AddInputChannel( _chan[_input], _chan_cred[_input] );
    }

  }

  for(int i = 0; i < ch_count; i++){
     for(int j = 0 ; j < 4; j++){
         cout << dbg[i][j] << " ";
     } cout << endl;
  }
  cout<<"Done links"<<endl;
}

void PolarFlyplusNew::InsertRandomFaults( const Configuration &config )
{
 
}

double PolarFlyplusNew::Capacity( ) const
{
  return (double)_k / 8.0;
}

void PolarFlyplusNew::RegisterRoutingFunctions(){

}
