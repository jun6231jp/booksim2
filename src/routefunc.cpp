// $Id$

/*
 Copyright (c) 2007-2015, Trustees of The Leland Stanford Junior University
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice, this 
 list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice, this
 list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

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

/*routefunc.cpp
 *
 *This is where most of the routing functions reside. Some of the topologies
 *has their own "register routing functions" which must be called to access
 *those routing functions. 
 *
 *After writing a routing function, don't forget to register it. The reg 
 *format is rfname_topologyname. 
 *
 */
#include <bitset>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <cassert>
#include <typeinfo>
#include "booksim.hpp"
#include "routefunc.hpp"
#include "kncube.hpp"
#include "hcube.hpp"
#include "polarflyplus.hpp"
#include "random_utils.hpp"
#include "misc_utils.hpp"
#include "fattree.hpp"
#include "tree4.hpp"
#include "qtree.hpp"
#include "cmesh.hpp"
#include "polarfly_tables.hpp"

#define VCNUM 3 //Total VC: VCNUM*2 (Request+Reply)
map<string, tRoutingFunction> gRoutingFunctionMap;
int  Hypercubeport,Polarflyport;
/* Global information used by routing functions */
int Faultescape=0;
int gNumVCs;

/* Add more functions here
 *
 */

// ============================================================
//  Balfour-Schultz
int gReadReqBeginVC, gReadReqEndVC;
int gWriteReqBeginVC, gWriteReqEndVC;
int gReadReplyBeginVC, gReadReplyEndVC;
int gWriteReplyBeginVC, gWriteReplyEndVC;

// ============================================================
//  QTree: Nearest Common Ancestor
// ===
void qtree_nca( const Router *r, const Flit *f,
		int in_channel, OutputSet* outputs, bool inject)
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int height = QTree::HeightFromID( r->GetID() );
    int pos    = QTree::PosFromID( r->GetID() );
    
    int dest   = f->dest;
    
    for (int i = height+1; i < gN; i++) 
      dest /= gK;
    if ( pos == dest / gK ) 
      // Route down to child
      out_port = dest % gK ; 
    else
      // Route up to parent
      out_port = gK;        

  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

// ============================================================
//  Tree4: Nearest Common Ancestor w/ Adaptive Routing Up
// ===
void tree4_anca( const Router *r, const Flit *f,
		 int in_channel, OutputSet* outputs, bool inject)
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int range = 1;
  
  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int dest = f->dest;
    
    const int NPOS = 16;
    
    int rH = r->GetID( ) / NPOS;
    int rP = r->GetID( ) % NPOS;
    
    if ( rH == 0 ) {
      dest /= 16;
      out_port = 2 * dest + RandomInt(1);
    } else if ( rH == 1 ) {
      dest /= 4;
      if ( dest / 4 == rP / 2 )
	out_port = dest % 4;
      else {
	out_port = gK;
	range = gK;
      }
    } else {
      if ( dest/4 == rP )
	out_port = dest % 4;
      else {
	out_port = gK;
	range = 2;
      }
    }
    
    //  cout << "Router("<<rH<<","<<rP<<"): id= " << f->id << " dest= " << f->dest << " out_port = "
    //       << out_port << endl;

  }

  outputs->Clear( );

  for (int i = 0; i < range; ++i) 
    outputs->AddRange( out_port + i, vcBegin, vcEnd );
}

// ============================================================
//  Tree4: Nearest Common Ancestor w/ Random Routing Up
// ===
void tree4_nca( const Router *r, const Flit *f,
		int in_channel, OutputSet* outputs, bool inject)
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int dest = f->dest;
    
    const int NPOS = 16;
    
    int rH = r->GetID( ) / NPOS;
    int rP = r->GetID( ) % NPOS;
    
    if ( rH == 0 ) {
      dest /= 16;
      out_port = 2 * dest + RandomInt(1);
    } else if ( rH == 1 ) {
      dest /= 4;
      if ( dest / 4 == rP / 2 )
	out_port = dest % 4;
      else
	out_port = gK + RandomInt(gK-1);
    } else {
      if ( dest/4 == rP )
	out_port = dest % 4;
      else
	out_port = gK + RandomInt(1);
    }
    
    //  cout << "Router("<<rH<<","<<rP<<"): id= " << f->id << " dest= " << f->dest << " out_port = "
    //       << out_port << endl;

  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

// ============================================================
//  FATTREE: Nearest Common Ancestor w/ Random  Routing Up
// ===
void fattree_nca( const Router *r, const Flit *f,
               int in_channel, OutputSet* outputs, bool inject)
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {
    
    int dest = f->dest;
    int router_id = r->GetID(); //routers are numbered with smallest at the top level
    int routers_per_level = powi(gK, gN-1);
    int pos = router_id%routers_per_level;
    int router_depth  = router_id/ routers_per_level; //which level
    int routers_per_neighborhood = powi(gK,gN-router_depth-1);
    int router_neighborhood = pos/routers_per_neighborhood; //coverage of this tree
    int router_coverage = powi(gK, gN-router_depth);  //span of the tree from this router
    

    //NCA reached going down
    if(dest <(router_neighborhood+1)* router_coverage && 
       dest >=router_neighborhood* router_coverage){
      //down ports are numbered first

      //ejection
      if(router_depth == gN-1){
	out_port = dest%gK;
      } else {	
	//find the down port for the destination
	int router_branch_coverage = powi(gK, gN-(router_depth+1)); 
	out_port = (dest-router_neighborhood* router_coverage)/router_branch_coverage;
      }
    } else {
      //up ports are numbered last
      assert(in_channel<gK);//came from a up channel
      out_port = gK+RandomInt(gK-1);
    }
  }  
  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

// ============================================================
//  FATTREE: Nearest Common Ancestor w/ Adaptive Routing Up
// ===
void fattree_anca( const Router *r, const Flit *f,
                int in_channel, OutputSet* outputs, bool inject)
{

  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));


  int out_port;

  if(inject) {

    out_port = -1;

  } else {


    int dest = f->dest;
    int router_id = r->GetID(); //routers are numbered with smallest at the top level
    int routers_per_level = powi(gK, gN-1);
    int pos = router_id%routers_per_level;
    int router_depth  = router_id/ routers_per_level; //which level
    int routers_per_neighborhood = powi(gK,gN-router_depth-1);
    int router_neighborhood = pos/routers_per_neighborhood; //coverage of this tree
    int router_coverage = powi(gK, gN-router_depth);  //span of the tree from this router
    

    //NCA reached going down
    if(dest <(router_neighborhood+1)* router_coverage && 
       dest >=router_neighborhood* router_coverage){
      //down ports are numbered first

      //ejection
      if(router_depth == gN-1){
	out_port = dest%gK;
      } else {	
	//find the down port for the destination
	int router_branch_coverage = powi(gK, gN-(router_depth+1)); 
	out_port = (dest-router_neighborhood* router_coverage)/router_branch_coverage;
      }
    } else {
      //up ports are numbered last
      assert(in_channel<gK);//came from a up channel
      out_port = gK;
      int random1 = RandomInt(gK-1); // Chose two ports out of the possible at random, compare loads, choose one.
      int random2 = RandomInt(gK-1);
      if (r->GetUsedCredit(out_port + random1) > r->GetUsedCredit(out_port + random2)){
	out_port = out_port + random2;
      }else{
	out_port =  out_port + random1;
      }
    }
  }  
  outputs->Clear( );
  
  outputs->AddRange( out_port, vcBegin, vcEnd );
}




// ============================================================
//  Mesh - adatpive XY,YX Routing 
//         pick xy or yx min routing adaptively at the source router
// ===

int dor_next_mesh( int cur, int dest, bool descending = false );

void adaptive_xy_yx_mesh( const Router *r, const Flit *f, 
		 int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else if(r->GetID() == f->dest) {

    // at destination router, we don't need to separate VCs by dim order
    out_port = 2*gN;

  } else {

    //each class must have at least 2 vcs assigned or else xy_yx will deadlock
    int const available_vcs = (vcEnd - vcBegin + 1) / 2;
    assert(available_vcs > 0);
    
    int out_port_xy = dor_next_mesh( r->GetID(), f->dest, false );
    int out_port_yx = dor_next_mesh( r->GetID(), f->dest, true );

    // Route order (XY or YX) determined when packet is injected
    //  into the network, adaptively
    bool x_then_y;
    if(in_channel < 2*gN){
      x_then_y =  (f->vc < (vcBegin + available_vcs));
    } else {
      int credit_xy = r->GetUsedCredit(out_port_xy);
      int credit_yx = r->GetUsedCredit(out_port_yx);
      if(credit_xy > credit_yx) {
	x_then_y = false;
      } else if(credit_xy < credit_yx) {
	x_then_y = true;
      } else {
	x_then_y = (RandomInt(1) > 0);
      }
    }
    
    if(x_then_y) {
      out_port = out_port_xy;
      vcEnd -= available_vcs;
    } else {
      out_port = out_port_yx;
      vcBegin += available_vcs;
    }

  }

  outputs->Clear();

  outputs->AddRange( out_port , vcBegin, vcEnd );
  
}

void xy_yx_mesh( const Router *r, const Flit *f, 
		 int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else if(r->GetID() == f->dest) {

    // at destination router, we don't need to separate VCs by dim order
    out_port = 2*gN;

  } else {

    //each class must have at least 2 vcs assigned or else xy_yx will deadlock
    int const available_vcs = (vcEnd - vcBegin + 1) / 2;
    assert(available_vcs > 0);

    // Route order (XY or YX) determined when packet is injected
    //  into the network
    bool x_then_y = ((in_channel < 2*gN) ?
		     (f->vc < (vcBegin + available_vcs)) :
		     (RandomInt(1) > 0));

    if(x_then_y) {
      out_port = dor_next_mesh( r->GetID(), f->dest, false );
      vcEnd -= available_vcs;
    } else {
      out_port = dor_next_mesh( r->GetID(), f->dest, true );
      vcBegin += available_vcs;
    }

  }

  outputs->Clear();

  outputs->AddRange( out_port , vcBegin, vcEnd );
  
}

//
// End Balfour-Schultz
//=============================================================

//=============================================================

int dor_next_mesh( int cur, int dest, bool descending )
{
  if ( cur == dest ) {
    return 2*gN;  // Eject
  }

  int dim_left;

  if(descending) {
    for ( dim_left = ( gN - 1 ); dim_left > 0; --dim_left ) {
      if ( ( cur * gK / gNodes ) != ( dest * gK / gNodes ) ) { break; }
      cur = (cur * gK) % gNodes; dest = (dest * gK) % gNodes;
    }
    cur = (cur * gK) / gNodes;
    dest = (dest * gK) / gNodes;
  } else {
    for ( dim_left = 0; dim_left < ( gN - 1 ); ++dim_left ) {
      if ( ( cur % gK ) != ( dest % gK ) ) { break; }
      cur /= gK; dest /= gK;
    }
    cur %= gK;
    dest %= gK;
  }

  if ( cur < dest ) {
    return 2*dim_left;     // Right
  } else {
    return 2*dim_left + 1; // Left
  }
}

//=============================================================

void dor_next_torus( int cur, int dest, int in_port,
		     int *out_port, int *partition,
		     bool balance = false )
{
  int dim_left;
  int dir;
  int dist2;

  for ( dim_left = 0; dim_left < gN; ++dim_left ) {
    if ( ( cur % gK ) != ( dest % gK ) ) { break; }
    cur /= gK; dest /= gK;
  }
  
  if ( dim_left < gN ) {

    if ( (in_port/2) != dim_left ) {
      // Turning into a new dimension

      cur %= gK; dest %= gK;
      dist2 = gK - 2 * ( ( dest - cur + gK ) % gK );
      
      if ( ( dist2 > 0 ) || 
	   ( ( dist2 == 0 ) && ( RandomInt( 1 ) ) ) ) {
	*out_port = 2*dim_left;     // Right
	dir = 0;
      } else {
	*out_port = 2*dim_left + 1; // Left
	dir = 1;
      }
      
      if ( partition ) {
	if ( balance ) {
	  // Cray's "Partition" allocation
	  // Two datelines: one between k-1 and 0 which forces VC 1
	  //                another between ((k-1)/2) and ((k-1)/2 + 1) which 
	  //                forces VC 0 otherwise any VC can be used
	  
	  if ( ( ( dir == 0 ) && ( cur > dest ) ) ||
	       ( ( dir == 1 ) && ( cur < dest ) ) ) {
	    *partition = 1;
	  } else if ( ( ( dir == 0 ) && ( cur <= (gK-1)/2 ) && ( dest >  (gK-1)/2 ) ) ||
		      ( ( dir == 1 ) && ( cur >  (gK-1)/2 ) && ( dest <= (gK-1)/2 ) ) ) {
	    *partition = 0;
	  } else {
	    *partition = RandomInt( 1 ); // use either VC set
	  }
	} else {
	  // Deterministic, fixed dateline between nodes k-1 and 0
	  
	  if ( ( ( dir == 0 ) && ( cur > dest ) ) ||
	       ( ( dir == 1 ) && ( dest < cur ) ) ) {
	    *partition = 1;
	  } else {
	    *partition = 0;
	  }
	}
      }
    } else {
      // Inverting the least significant bit keeps
      // the packet moving in the same direction
      *out_port = in_port ^ 0x1;
    }    

  } else {
    *out_port = 2*gN;  // Eject
  }
}

//=============================================================

void dim_order_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
	cout << "     routefunc mesh routing start : cur " << r << " dest " << f->dest << " in_port " << in_channel << " in_vc " << f->vc << endl;
  int out_port = inject ? -1 : dor_next_mesh( r->GetID( ), f->dest );
  
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  if ( !inject && f->watch ) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
	       << "Adding VC range [" 
	       << vcBegin << "," 
	       << vcEnd << "]"
	       << " at output port " << out_port
	       << " for flit " << f->id
	       << " (input port " << in_channel
	       << ", destination " << f->dest << ")"
	       << "." << endl;
  }
  
  outputs->Clear();

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void dim_order_hcube( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int out_port = inject ? -1 : dor_next_mesh( r->GetID( ), f->dest );

  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  if ( !inject && f->watch ) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
               << "Adding VC range ["
               << vcBegin << ","
               << vcEnd << "]"
               << " at output port " << out_port
               << " for flit " << f->id
               << " (input port " << in_channel
               << ", destination " << f->dest << ")"
               << "." << endl;
  }
    outputs->Clear();

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void dim_order_ni_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int out_port = inject ? -1 : dor_next_mesh( r->GetID( ), f->dest );
  
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  // at the destination router, we don't need to separate VCs by destination
  if(inject || (r->GetID() != f->dest)) {

    int const vcs_per_dest = (vcEnd - vcBegin + 1) / gNodes;
    assert(vcs_per_dest > 0);

    vcBegin += f->dest * vcs_per_dest;
    vcEnd = vcBegin + vcs_per_dest - 1;

  }
  
  if( !inject && f->watch ) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
	       << "Adding VC range [" 
	       << vcBegin << "," 
	       << vcEnd << "]"
	       << " at output port " << out_port
	       << " for flit " << f->id
	       << " (input port " << in_channel
	       << ", destination " << f->dest << ")"
	       << "." << endl;
  }
  
  outputs->Clear( );
  
  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void dim_order_pni_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int out_port = inject ? -1 : dor_next_mesh( r->GetID(), f->dest );
  
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  if(inject || (r->GetID() != f->dest)) {
    int next_coord = f->dest;
    if(!inject) {
      int out_dim = out_port / 2;
      for(int d = 0; d < out_dim; ++d) {
	next_coord /= gK;
      }
    }
    next_coord %= gK;
    assert(next_coord >= 0 && next_coord < gK);
    int vcs_per_dest = (vcEnd - vcBegin + 1) / gK;
    assert(vcs_per_dest > 0);
    vcBegin += next_coord * vcs_per_dest;
    vcEnd = vcBegin + vcs_per_dest - 1;
  }

  if( !inject && f->watch ) {
    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
	       << "Adding VC range [" 
	       << vcBegin << "," 
	       << vcEnd << "]"
	       << " at output port " << out_port
	       << " for flit " << f->id
	       << " (input port " << in_channel
	       << ", destination " << f->dest << ")"
	       << "." << endl;
  }
  
  outputs->Clear( );
  
  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

// Random intermediate in the minimal quadrant defined
// by the source and destination
int rand_min_intr_mesh( int src, int dest )
{
  int dist;

  int intm = 0;
  int offset = 1;

  for ( int n = 0; n < gN; ++n ) {
    dist = ( dest % gK ) - ( src % gK );

    if ( dist > 0 ) {
      intm += offset * ( ( src % gK ) + RandomInt( dist ) );
    } else {
      intm += offset * ( ( dest % gK ) + RandomInt( -dist ) );
    }

    offset *= gK;
    dest /= gK; src /= gK;
  }

  return intm;
}

//=============================================================

void romm_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    if ( in_channel == 2*gN ) {
      f->ph   = 0;  // Phase 0
      f->intm = rand_min_intr_mesh( f->src, f->dest );
    } 

    if ( ( f->ph == 0 ) && ( r->GetID( ) == f->intm ) ) {
      f->ph = 1; // Go to phase 1
    }

    out_port = dor_next_mesh( r->GetID( ), (f->ph == 0) ? f->intm : f->dest );

    // at the destination router, we don't need to separate VCs by phase
    if(r->GetID() != f->dest) {

      //each class must have at least 2 vcs assigned or else valiant valiant will deadlock
      int available_vcs = (vcEnd - vcBegin + 1) / 2;
      assert(available_vcs > 0);

      if(f->ph == 0) {
	vcEnd -= available_vcs;
      } else {
	assert(f->ph == 1);
	vcBegin += available_vcs;
      }
    }

  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void romm_ni_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  // at the destination router, we don't need to separate VCs by destination
  if(inject || (r->GetID() != f->dest)) {

    int const vcs_per_dest = (vcEnd - vcBegin + 1) / gNodes;
    assert(vcs_per_dest > 0);

    vcBegin += f->dest * vcs_per_dest;
    vcEnd = vcBegin + vcs_per_dest - 1;

  }

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    if ( in_channel == 2*gN ) {
      f->ph   = 0;  // Phase 0
      f->intm = rand_min_intr_mesh( f->src, f->dest );
    } 

    if ( ( f->ph == 0 ) && ( r->GetID( ) == f->intm ) ) {
      f->ph = 1; // Go to phase 1
    }

    out_port = dor_next_mesh( r->GetID( ), (f->ph == 0) ? f->intm : f->dest );

  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void min_adapt_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  outputs->Clear( );
  
  if(inject) {
    // injection can use all VCs
    outputs->AddRange(-1, vcBegin, vcEnd);
    return;
  } else if(r->GetID() == f->dest) {
    // ejection can also use all VCs
    outputs->AddRange(2*gN, vcBegin, vcEnd);
    return;
  }

  int in_vc;

  if ( in_channel == 2*gN ) {
    in_vc = vcEnd; // ignore the injection VC
  } else {
    in_vc = f->vc;
  }
  
  // DOR for the escape channel (VC 0), low priority 
  int out_port = dor_next_mesh( r->GetID( ), f->dest );    
  outputs->AddRange( out_port, 0, vcBegin, vcBegin );
  
  if ( f->watch ) {
      *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		  << "Adding VC range [" 
		  << vcBegin << "," 
		  << vcBegin << "]"
		  << " at output port " << out_port
		  << " for flit " << f->id
		  << " (input port " << in_channel
		  << ", destination " << f->dest << ")"
		  << "." << endl;
   }
  
  if ( in_vc != vcBegin ) { // If not in the escape VC
    // Minimal adaptive for all other channels
    int cur = r->GetID( );
    int dest = f->dest;
    
    for ( int n = 0; n < gN; ++n ) {
      if ( ( cur % gK ) != ( dest % gK ) ) { 
	// Add minimal direction in dimension 'n'
	if ( ( cur % gK ) < ( dest % gK ) ) { // Right
	  if ( f->watch ) {
	    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
			<< "Adding VC range [" 
		       << (vcBegin+1) << "," 
			<< vcEnd << "]"
			<< " at output port " << 2*n
			<< " with priority " << 1
			<< " for flit " << f->id
			<< " (input port " << in_channel
			<< ", destination " << f->dest << ")"
			<< "." << endl;
	  }
	  outputs->AddRange( 2*n, vcBegin+1, vcEnd, 1 ); 
	} else { // Left
	  if ( f->watch ) {
	    *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
			<< "Adding VC range [" 
		       << (vcBegin+1) << "," 
			<< vcEnd << "]"
			<< " at output port " << 2*n+1
			<< " with priority " << 1
			<< " for flit " << f->id
			<< " (input port " << in_channel
			<< ", destination " << f->dest << ")"
			<< "." << endl;
	  }
	  outputs->AddRange( 2*n + 1, vcBegin+1, vcEnd, 1 ); 
	}
      }
      cur  /= gK;
      dest /= gK;
    }
  } 
}

//=============================================================

void planar_adapt_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  outputs->Clear( );
  
  if(inject) {
    // injection can use all VCs
    outputs->AddRange(-1, vcBegin, vcEnd);
    return;
  }

  int cur     = r->GetID( ); 
  int dest    = f->dest;

  if ( cur != dest ) {
   
    int in_vc   = f->vc;
    int vc_mult = (vcEnd - vcBegin + 1) / 3;

    // Find the first unmatched dimension -- except
    // for when we're in the first dimension because
    // of misrouting in the last adaptive plane.
    // In this case, go to the last dimension instead.

    int n;
    for ( n = 0; n < gN; ++n ) {
      if ( ( ( cur % gK ) != ( dest % gK ) ) &&
	   !( ( in_channel/2 == 0 ) &&
	      ( n == 0 ) &&
	      ( in_vc < vcBegin+2*vc_mult ) ) ) {
	break;
      }

      cur  /= gK;
      dest /= gK;
    }

    assert( n < gN );

    if ( f->watch ) {
      *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		  << "PLANAR ADAPTIVE: flit " << f->id 
		  << " in adaptive plane " << n << "." << endl;
    }

    // We're in adaptive plane n

    // Can route productively in d_{i,2}
    bool increase;
    bool fault;
    if ( ( cur % gK ) < ( dest % gK ) ) { // Increasing
      increase = true;
      if ( !r->IsFaultyOutput( 2*n ) ) {
	outputs->AddRange( 2*n, vcBegin+2*vc_mult, vcEnd );
	fault = false;

	if ( f->watch ) {
	  *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		      << "PLANAR ADAPTIVE: increasing in dimension " << n
		      << "." << endl;
	}
      } else {
	fault = true;
      }
    } else { // Decreasing
      increase = false;
      if ( !r->IsFaultyOutput( 2*n + 1 ) ) {
	outputs->AddRange( 2*n + 1, vcBegin+2*vc_mult, vcEnd ); 
	fault = false;

	if ( f->watch ) {
	  *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		      << "PLANAR ADAPTIVE: decreasing in dimension " << n
		      << "." << endl;
	}
      } else {
	fault = true;
      }
    }
      
    n = ( n + 1 ) % gN;
    cur  /= gK;
    dest /= gK;
      
    if ( !increase ) {
      vcBegin += vc_mult;
    }
    vcEnd = vcBegin + vc_mult - 1;
      
    int d1_min_c;
    if ( ( cur % gK ) < ( dest % gK ) ) { // Increasing in d_{i+1}
      d1_min_c = 2*n;
    } else if ( ( cur % gK ) != ( dest % gK ) ) {  // Decreasing in d_{i+1}
      d1_min_c = 2*n + 1;
    } else {
      d1_min_c = -1;
    }
      
    // do we want to 180?  if so, the last
    // route was a misroute in this dimension,
    // if there is no fault in d_i, just ignore
    // this dimension, otherwise continue to misroute
    if ( d1_min_c == in_channel ) { 
      if ( fault ) {
	d1_min_c = in_channel ^ 1;
      } else {
	d1_min_c = -1;
      }

      if ( f->watch ) {
	*gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		    << "PLANAR ADAPTIVE: avoiding 180 in dimension " << n
		    << "." << endl;
      }
    }
      
    if ( d1_min_c != -1 ) {
      if ( !r->IsFaultyOutput( d1_min_c ) ) {
	outputs->AddRange( d1_min_c, vcBegin, vcEnd );
      } else if ( fault ) {
	// major problem ... fault in d_i and d_{i+1}
	r->Error( "There seem to be faults in d_i and d_{i+1}" );
      }
    } else if ( fault ) { // need to misroute!
      bool atedge;
      if ( cur % gK == 0 ) {
	d1_min_c = 2*n;
	atedge = true;
      } else if ( cur % gK == gK - 1 ) {
	d1_min_c = 2*n + 1;
	atedge = true;
      } else {
	d1_min_c = 2*n + RandomInt( 1 ); // random misroute

	if ( d1_min_c  == in_channel ) { // don't 180
	  d1_min_c = in_channel ^ 1;
	}
	atedge = false;
      }
      
      if ( !r->IsFaultyOutput( d1_min_c ) ) {
	outputs->AddRange( d1_min_c, vcBegin, vcEnd );
      } else if ( !atedge && !r->IsFaultyOutput( d1_min_c ^ 1 ) ) {
	outputs->AddRange( d1_min_c ^ 1, vcBegin, vcEnd );
      } else {
	// major problem ... fault in d_i and d_{i+1}
	r->Error( "There seem to be faults in d_i and d_{i+1}" );
      }
    }
  } else {
    outputs->AddRange( 2*gN, vcBegin, vcEnd ); 
  }
}

//=============================================================
/*
  FIXME: This is broken (note that f->dr is never actually modified).
  Even if it were, this should really use f->ph instead of introducing a single-
  use field.

void limited_adapt_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  outputs->Clear( );

  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  if ( inject ) {
    outputs->AddRange( -1, vcBegin, vcEnd - 1 );
    f->dr = 0; // zero dimension reversals
    return;
  }

  int cur = r->GetID( );
  int dest = f->dest;
  
  if ( cur != dest ) {
    if ( ( f->vc != vcEnd ) && 
	 ( f->dr != vcEnd - 1 ) ) {
      
      for ( int n = 0; n < gN; ++n ) {
	if ( ( cur % gK ) != ( dest % gK ) ) { 
	  int min_port;
	  if ( ( cur % gK ) < ( dest % gK ) ) { 
	    min_port = 2*n; // Right
	  } else {
	    min_port = 2*n + 1; // Left
	  }
	  
	  // Go in a productive direction with high priority
	  outputs->AddRange( min_port, vcBegin, vcEnd - 1, 2 );
	  
	  // Go in the non-productive direction with low priority
	  outputs->AddRange( min_port ^ 0x1, vcBegin, vcEnd - 1, 1 );
	} else {
	  // Both directions are non-productive
	  outputs->AddRange( 2*n, vcBegin, vcEnd - 1, 1 );
	  outputs->AddRange( 2*n+1, vcBegin, vcEnd - 1, 1 );
	}
	
	cur  /= gK;
	dest /= gK;
      }
      
    } else {
      outputs->AddRange( dor_next_mesh( cur, dest ),
			 vcEnd, vcEnd, 0 );
    }
    
  } else { // at destination
    outputs->AddRange( 2*gN, vcBegin, vcEnd ); 
  }
}
*/
//=============================================================

void valiant_mesh( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    if ( in_channel == 2*gN ) {
      f->ph   = 0;  // Phase 0
      f->intm = RandomInt( gNodes - 1 );
    }

    if ( ( f->ph == 0 ) && ( r->GetID( ) == f->intm ) ) {
      f->ph = 1; // Go to phase 1
    }

    out_port = dor_next_mesh( r->GetID( ), (f->ph == 0) ? f->intm : f->dest );

    // at the destination router, we don't need to separate VCs by phase
    if(r->GetID() != f->dest) {

      //each class must have at least 2 vcs assigned or else valiant valiant will deadlock
      int const available_vcs = (vcEnd - vcBegin + 1) / 2;
      assert(available_vcs > 0);

      if(f->ph == 0) {
	vcEnd -= available_vcs;
      } else {
	assert(f->ph == 1);
	vcBegin += available_vcs;
      }
    }

  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void valiant_torus( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int phase;
    if ( in_channel == 2*gN ) {
      phase   = 0;  // Phase 0
      f->intm = RandomInt( gNodes - 1 );
    } else {
      phase = f->ph / 2;
    }

    if ( ( phase == 0 ) && ( r->GetID( ) == f->intm ) ) {
      phase = 1; // Go to phase 1
      in_channel = 2*gN; // ensures correct vc selection at the beginning of phase 2
    }
  
    int ring_part;
    dor_next_torus( r->GetID( ), (phase == 0) ? f->intm : f->dest, in_channel,
		    &out_port, &ring_part, false );

    f->ph = 2 * phase + ring_part;

    // at the destination router, we don't need to separate VCs by phase, etc.
    if(r->GetID() != f->dest) {

      int const ring_available_vcs = (vcEnd - vcBegin + 1) / 2;
      assert(ring_available_vcs > 0);

      if(ring_part == 0) {
	vcEnd -= ring_available_vcs;
      } else {
	assert(ring_part == 1);
	vcBegin += ring_available_vcs;
      }

      int const ph_available_vcs = ring_available_vcs / 2;
      assert(ph_available_vcs > 0);

      if(phase == 0) {
	vcEnd -= ph_available_vcs;
      } else {
	assert(phase == 1);
	vcBegin += ph_available_vcs;
      }
    }

  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void valiant_ni_torus( const Router *r, const Flit *f, int in_channel, 
		       OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  // at the destination router, we don't need to separate VCs by destination
  if(inject || (r->GetID() != f->dest)) {

    int const vcs_per_dest = (vcEnd - vcBegin + 1) / gNodes;
    assert(vcs_per_dest > 0);

    vcBegin += f->dest * vcs_per_dest;
    vcEnd = vcBegin + vcs_per_dest - 1;

  }

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int phase;
    if ( in_channel == 2*gN ) {
      phase   = 0;  // Phase 0
      f->intm = RandomInt( gNodes - 1 );
    } else {
      phase = f->ph / 2;
    }

    if ( ( f->ph == 0 ) && ( r->GetID( ) == f->intm ) ) {
      f->ph = 1; // Go to phase 1
      in_channel = 2*gN; // ensures correct vc selection at the beginning of phase 2
    }
  
    int ring_part;
    dor_next_torus( r->GetID( ), (f->ph == 0) ? f->intm : f->dest, in_channel,
		    &out_port, &ring_part, false );

    f->ph = 2 * phase + ring_part;

    // at the destination router, we don't need to separate VCs by phase, etc.
    if(r->GetID() != f->dest) {

      int const ring_available_vcs = (vcEnd - vcBegin + 1) / 2;
      assert(ring_available_vcs > 0);

      if(ring_part == 0) {
	vcEnd -= ring_available_vcs;
      } else {
	assert(ring_part == 1);
	vcBegin += ring_available_vcs;
      }

      int const ph_available_vcs = ring_available_vcs / 2;
      assert(ph_available_vcs > 0);

      if(phase == 0) {
	vcEnd -= ph_available_vcs;
      } else {
	assert(phase == 1);
	vcBegin += ph_available_vcs;
      }
    }

    if (f->watch) {
      *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		 << "Adding VC range [" 
		 << vcBegin << "," 
		 << vcEnd << "]"
		 << " at output port " << out_port
		 << " for flit " << f->id
		 << " (input port " << in_channel
		 << ", destination " << f->dest << ")"
		 << "." << endl;
    }

  }
  
  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void dim_order_torus( const Router *r, const Flit *f, int in_channel, 
		      OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {
    
    int cur  = r->GetID( );
    int dest = f->dest;

    dor_next_torus( cur, dest, in_channel,
		    &out_port, &f->ph, false );


    // at the destination router, we don't need to separate VCs by ring partition
    if(cur != dest) {

      int const available_vcs = (vcEnd - vcBegin + 1) / 2;
      assert(available_vcs > 0);

      if ( f->ph == 0 ) {
	vcEnd -= available_vcs;
      } else {
	vcBegin += available_vcs;
      } 
    }

    if ( f->watch ) {
      *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		 << "Adding VC range [" 
		 << vcBegin << "," 
		 << vcEnd << "]"
		 << " at output port " << out_port
		 << " for flit " << f->id
		 << " (input port " << in_channel
		 << ", destination " << f->dest << ")"
		 << "." << endl;
    }

  }
 
  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void dim_order_ni_torus( const Router *r, const Flit *f, int in_channel, 
			 OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {
    
    int cur  = r->GetID( );
    int dest = f->dest;

    dor_next_torus( cur, dest, in_channel,
		    &out_port, NULL, false );

    // at the destination router, we don't need to separate VCs by destination
    if(cur != dest) {

      int const vcs_per_dest = (vcEnd - vcBegin + 1) / gNodes;
      assert(vcs_per_dest);

      vcBegin += f->dest * vcs_per_dest;
      vcEnd = vcBegin + vcs_per_dest - 1;

    }

    if ( f->watch ) {
      *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		 << "Adding VC range [" 
		 << vcBegin << "," 
		 << vcEnd << "]"
		 << " at output port " << out_port
		 << " for flit " << f->id
		 << " (input port " << in_channel
		 << ", destination " << f->dest << ")"
		 << "." << endl;
    }

  }
  
  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void dim_order_bal_torus( const Router *r, const Flit *f, int in_channel, 
			  OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int cur  = r->GetID( );
    int dest = f->dest;

    dor_next_torus( cur, dest, in_channel,
		    &out_port, &f->ph, true );

    // at the destination router, we don't need to separate VCs by ring partition
    if(cur != dest) {

      int const available_vcs = (vcEnd - vcBegin + 1) / 2;
      assert(available_vcs > 0);

      if ( f->ph == 0 ) {
	vcEnd -= available_vcs;
      } else {
	assert(f->ph == 1);
	vcBegin += available_vcs;
      } 
    }

    if ( f->watch ) {
      *gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
		 << "Adding VC range [" 
		 << vcBegin << "," 
		 << vcEnd << "]"
		 << " at output port " << out_port
		 << " for flit " << f->id
		 << " (input port " << in_channel
		 << ", destination " << f->dest << ")"
		 << "." << endl;
    }

  }
  
  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}

//=============================================================

void min_adapt_torus( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  outputs->Clear( );

  if(inject) {
    // injection can use all VCs
    outputs->AddRange(-1, vcBegin, vcEnd);
    return;
  } else if(r->GetID() == f->dest) {
    // ejection can also use all VCs
    outputs->AddRange(2*gN, vcBegin, vcEnd);
  }

  int in_vc;
  if ( in_channel == 2*gN ) {
    in_vc = vcEnd; // ignore the injection VC
  } else {
    in_vc = f->vc;
  }
  
  int cur = r->GetID( );
  int dest = f->dest;

  int out_port;

  if ( in_vc > ( vcBegin + 1 ) ) { // If not in the escape VCs
    // Minimal adaptive for all other channels
    
    for ( int n = 0; n < gN; ++n ) {
      if ( ( cur % gK ) != ( dest % gK ) ) {
	int dist2 = gK - 2 * ( ( ( dest % gK ) - ( cur % gK ) + gK ) % gK );
	
	if ( dist2 > 0 ) { /*) || 
			     ( ( dist2 == 0 ) && ( RandomInt( 1 ) ) ) ) {*/
	  outputs->AddRange( 2*n, vcBegin+3, vcBegin+3, 1 ); // Right
	} else {
	  outputs->AddRange( 2*n + 1, vcBegin+3, vcBegin+3, 1 ); // Left
	}
      }

      cur  /= gK;
      dest /= gK;
    }
    
    // DOR for the escape channel (VCs 0-1), low priority --- 
    // trick the algorithm with the in channel.  want VC assignment
    // as if we had injected at this node
    dor_next_torus( r->GetID( ), f->dest, 2*gN,
		    &out_port, &f->ph, false );
  } else {
    // DOR for the escape channel (VCs 0-1), low priority 
    dor_next_torus( cur, dest, in_channel,
		    &out_port, &f->ph, false );
  }

  if ( f->ph == 0 ) {
    outputs->AddRange( out_port, vcBegin, vcBegin, 0 );
  } else  {
    outputs->AddRange( out_port, vcBegin+1, vcBegin+1, 0 );
  } 
}

//=============================================================

void dest_tag_fly( const Router *r, const Flit *f, int in_channel, 
		   OutputSet *outputs, bool inject )
{
  int vcBegin = 0, vcEnd = gNumVCs-1;
  if ( f->type == Flit::READ_REQUEST ) {
    vcBegin = gReadReqBeginVC;
    vcEnd = gReadReqEndVC;
  } else if ( f->type == Flit::WRITE_REQUEST ) {
    vcBegin = gWriteReqBeginVC;
    vcEnd = gWriteReqEndVC;
  } else if ( f->type ==  Flit::READ_REPLY ) {
    vcBegin = gReadReplyBeginVC;
    vcEnd = gReadReplyEndVC;
  } else if ( f->type ==  Flit::WRITE_REPLY ) {
    vcBegin = gWriteReplyBeginVC;
    vcEnd = gWriteReplyEndVC;
  }
  assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

  int out_port;

  if(inject) {

    out_port = -1;

  } else {

    int stage = ( r->GetID( ) * gK ) / gNodes;
    int dest  = f->dest;

    while( stage < ( gN - 1 ) ) {
      dest /= gK;
      ++stage;
    }

    out_port = dest % gK;
  }

  outputs->Clear( );

  outputs->AddRange( out_port, vcBegin, vcEnd );
}



//=============================================================

void chaos_torus( const Router *r, const Flit *f, 
		  int in_channel, OutputSet *outputs, bool inject )
{
  outputs->Clear( );

  if(inject) {
    outputs->AddRange(-1, 0, 0);
    return;
  }

  int cur = r->GetID( );
  int dest = f->dest;
  
  if ( cur != dest ) {
    for ( int n = 0; n < gN; ++n ) {

      if ( ( cur % gK ) != ( dest % gK ) ) { 
	int dist2 = gK - 2 * ( ( ( dest % gK ) - ( cur % gK ) + gK ) % gK );
      
	if ( dist2 >= 0 ) {
	  outputs->AddRange( 2*n, 0, 0 ); // Right
	} 
	
	if ( dist2 <= 0 ) {
	  outputs->AddRange( 2*n + 1, 0, 0 ); // Left
	}
      }

      cur  /= gK;
      dest /= gK;
    }
  } else {
    outputs->AddRange( 2*gN, 0, 0 ); 
  }
}

void chaos_mesh( const Router *r, const Flit *f, 
		  int in_channel, OutputSet *outputs, bool inject )
{
  outputs->Clear( );

  if(inject) {
    outputs->AddRange(-1, 0, 0);
    return;
  }

  int cur = r->GetID( );
  int dest = f->dest;
  
  if ( cur != dest ) {
    for ( int n = 0; n < gN; ++n ) {
      if ( ( cur % gK ) != ( dest % gK ) ) { 
	// Add minimal direction in dimension 'n'
	if ( ( cur % gK ) < ( dest % gK ) ) { // Right
	  outputs->AddRange( 2*n, 0, 0 ); 
	} else { // Left
	  outputs->AddRange( 2*n + 1, 0, 0 ); 
	}
      }
      cur  /= gK;
      dest /= gK;
    }
  } else {
    outputs->AddRange( 2*gN, 0, 0 ); 
  }
}
//=======================polarfly+ routing functions=========================
int polarport_cal(int src_grp, int dest_grp){
   
    int global_port=-1;

    // 1hop
    for(int i=0; i < Polarflyport; i++){
        if(polarfly_connection_table[src_grp][i]==dest_grp){
	    if (src_grp == dest_grp){//red group self connection
                global_port = (i + 1) % Polarflyport;
	    }
	    else {
                global_port = i;
            }
	    break;
        }
    }
    // 2hop
    if(global_port==-1){
       for(int i=0; i < Polarflyport; i++){
          for(int j=0; j < Polarflyport; j++){
              if(polarfly_connection_table[src_grp][i]==polarfly_connection_table[dest_grp][j]){
                  global_port = i;
                  break;
              }
          }
          if(global_port > -1){break;}
       }
    }
    assert(global_port>-1);
    global_port += Hypercubeport + 1; //local port offset + eject port offset
    return global_port;
}

//for request
int* order0;
int* order1;
int* order2;
//for reply
int* order3;
int* order4;
int* order5;

void make_order(){
  order0=(int*)malloc(sizeof(int)*Hypercubeport);
  order1=(int*)malloc(sizeof(int)*Hypercubeport);
  order2=(int*)malloc(sizeof(int)*Hypercubeport);
  order3=(int*)malloc(sizeof(int)*Hypercubeport);
  order4=(int*)malloc(sizeof(int)*Hypercubeport);
  order5=(int*)malloc(sizeof(int)*Hypercubeport);
  for(int i = 0 ; i < Hypercubeport; i++){
     order0[i]=i;
     order1[i]=Hypercubeport-1-i;
     order2[i]=i;
     order3[i]=Hypercubeport-1-i;
     order4[i]=i;
     order5[i]=Hypercubeport-1-i;
  }
}


tuple<int,int> hyperport_cal(int hypercube_mv, int in_dim_order, int in_vc){
     int out_port = -1;
     int out_dim_order = 0;
     for(int i = in_dim_order; i < Hypercubeport; i++){
        int dim = -1; 
     	if(in_vc==0){dim = order0[i];}
	if(in_vc==1){dim = order1[i];}
	if(in_vc==2){dim = order2[i];}
	if(in_vc==VCNUM){dim = order3[i];}
	if(in_vc==VCNUM+1){dim = order4[i];}
        if(in_vc==VCNUM+2){dim = order5[i];}
 	if(((hypercube_mv >> dim) & 1) == 1){
	     out_dim_order = i+1;
	     out_port = dim+1;
             break;
	}
     }
     return {out_dim_order,out_port};
}

/* ============== adaptive routing =====================
int polarfly_fault_escape(int src_grp, int in_port, int global_port){
     int out_port=-1;
     if (in_port > Hypercubeport) {//Global receive
        for (int i = 0 ; i < Polarflyport; i++){
           if(((i+Hypercubeport+1) != in_port) && ((i+Hypercubeport+1) != global_port)) out_port = i;
	   break;
	}
     }
     else{ //Local receive
        out_port = (global_port-Hypercubeport-1)+1; //increment
        if(out_port >= Polarflyport - 1 && src_grp < Polarflyport){
	     out_port = out_port % (Polarflyport - 1);
        } //red group 
        if(out_port >= Polarflyport && src_grp >= Polarflyport){
	     out_port = out_port % Polarflyport;
        } //green, blue group
     }
     out_port += (Hypercubeport + 1); //Local offset + inject offset
     return out_port;
}
*/

//===================source routing=====================
int traffic_table[total_node][node_port] = {0};
struct LocalMoveResult {
    int current;
    int local_move;
    bool routing_complete;
    bool fault_detected;
};

LocalMoveResult process_local_move(int current, int dest, int hypercube_mv, int in_vc, int id) {
    LocalMoveResult result = {current, 0, false, false};
    int local_port = 0;
    int dim_order=0;
    if (hypercube_mv == 0) {
        return result;
    }
    int mv = 0;
    for (int k = 0; k < Hypercubeport; k++) {
	auto [dim_order_temp, local_port_temp]= hyperport_cal(hypercube_mv, dim_order, in_vc);
	if (local_port_temp >= 0 && !fault_table[current][local_port_temp]) {
            local_port = local_port_temp;
	    dim_order = dim_order_temp;
	    cout << "source routing id:" << id << " local:" << (in_vc%VCNUM)+1 << " node" << current << " -> ";
            current ^= (1 << (local_port - 1));
            mv ^= (1 << (local_port - 1));
	    cout << "node" << current << endl;
        }
        else if (fault_table[current][local_port_temp]) {
            cout << "source routing id:" << id << " local:" << (in_vc%VCNUM)+1 << " node" << current << " port" << local_port_temp << " err" << endl;
            result.fault_detected = true;
            break;
        }
        else {
            break;
        }
    }
    if (!result.fault_detected) {
       result.local_move = mv;
    }
    if (bitmask(current, Hypercubeport) == bitmask(dest, Hypercubeport)) {
        result.routing_complete = true;
    }
    result.current = current;
    return result;
}

struct GlobalMoveResult {
    int current;
    int current_group;
    int global_port;
    bool routing_complete;
};

GlobalMoveResult process_global_move(int current, int current_group, int dest_group, int phase, int id) {
    GlobalMoveResult result = {current, current_group, 0, false};
    
    int global_port = polarport_cal(current_group, dest_group);
    if (!fault_table[current][global_port]) {
	cout << "source routing id:" << id << " grp" << current_group << " -> ";
        current_group = polarfly_connection_table[current_group][global_port - Hypercubeport - 1];
        cout << "grp" << current_group << endl;
	current = current_group * (1 << Hypercubeport) + bitmask(current, Hypercubeport);
        
        if (current_group == dest_group) {
            result.routing_complete = true;
        }
        result.global_port = global_port;
    }
    else {
        cout << "source routing id:" << id << " global" << phase << " node" 
             << current
             << " port" << global_port << " err" << endl;
    }
    
    result.current = current;
    result.current_group = current_group;
    return result;
}

void source_routing(const Flit *f, int current_node, int destination_node) {
    int local_move1 = 0, local_move2 = 0, local_move3 = 0;
    int global_port1 = 0, global_port2 = 0;
    int routing_result = 0;
    int min_ok_weight = 0;
    vector<vector<int>> oklist; //localmv1, localmv2, localmv3, global1, blobal2, weight
    int in_vc=f->vc;
    if ( f->type == Flit::READ_REPLY || f->type == Flit::WRITE_REPLY) {
        in_vc+=VCNUM; //reply inject
    }
    int hypercube_moves = bitmask(current_node, Hypercubeport) ^ bitmask(destination_node, Hypercubeport);
    bitset<8> bit_min(hypercube_moves);
    int min_weight = bit_min.count();

    make_order();
    for (int i = 0; i < (1<<Hypercubeport); i++) {
        bool local_routing_complete = false;
        bool global_routing_complete = false;
	for (int j = 0 ; j < (1<<Hypercubeport); j++) {
	    local_move1 = 0; local_move2 = 0; local_move3 = 0;
	    global_port1 = 0; global_port2 = 0;
            int current = current_node;
            const int dest = destination_node;
            int current_group = current >> Hypercubeport;
            const int dest_group = dest >> Hypercubeport;
            int hypercube_mv3 = j;
	    int hypercube_mv2 = i;
	    int hypercube_mv1 = hypercube_moves^i^j; 
	    bitset<8> bit1(hypercube_mv1);
	    bitset<8> bit2(hypercube_mv2);
	    bitset<8> bit3(hypercube_mv3);
	    int weight = bit1.count() + bit2.count() + bit3.count();
            cout << "source routing id:" << f->pid << " hypercube hops:" << bitset<8>(hypercube_mv1) << " " << bitset<8>(hypercube_mv2) << " " << bitset<8>(hypercube_mv3) << endl;

            local_routing_complete = (bitmask(current, Hypercubeport) == bitmask(dest, Hypercubeport));
            global_routing_complete = (current_group == dest_group);

            // local move 1
            if (hypercube_mv1 > 0) {
                auto result = process_local_move(current, dest, hypercube_mv1, in_vc, f->pid);
                if (result.fault_detected) continue;
                current = result.current;
                local_move1 = result.local_move;
                local_routing_complete = result.routing_complete;
            }
            if (local_routing_complete && global_routing_complete) {
	        vector<int> okpath = {local_move1,local_move2,local_move3,global_port1,global_port2,weight};
                oklist.push_back(okpath);
                continue;	
	    }
            // global move 1
            if (!global_routing_complete) {
                auto result = process_global_move(current, current_group, dest_group, 1, f->pid);
                if (result.global_port == 0) continue;
                current = result.current;
                current_group = result.current_group;
                global_port1 = result.global_port;
                global_routing_complete = result.routing_complete;
            }
	    if((global_port1 == 0) && (local_move1 == 0))continue;
	    if (local_routing_complete && global_routing_complete){
                vector<int> okpath = {local_move1,local_move2,local_move3,global_port1,global_port2,weight};
                oklist.push_back(okpath);
                continue;
	    }
            // local move 2
            if (hypercube_mv2 > 0) {
                auto result = process_local_move(current, dest, hypercube_mv2, in_vc+1, f->pid);
                if (result.fault_detected) continue;
                current = result.current;
                local_move2 = result.local_move;
                local_routing_complete = result.routing_complete;
            }
	    if (local_routing_complete && global_routing_complete) {
		vector<int> okpath = {local_move1,local_move2,local_move3,global_port1,global_port2,weight};
                oklist.push_back(okpath);
                continue;
	    }
            // global move 2
            if (!global_routing_complete) {
                auto result = process_global_move(current, current_group, dest_group, 2, f->pid);
                if (result.global_port == 0) continue;
                current = result.current;
                current_group = result.current_group;
                global_port2 = result.global_port;
                global_routing_complete = result.routing_complete;
            }
	    if((global_port2 == 0) && (local_move2 == 0))continue;
	    if (local_routing_complete && global_routing_complete){
	        vector<int> okpath = {local_move1,local_move2,local_move3,global_port1,global_port2,weight};
                oklist.push_back(okpath);
                continue;
	    }
            // local move 3
            if (hypercube_mv3 > 0) {
                auto result = process_local_move(current, dest, hypercube_mv3, in_vc+2, f->pid);
                if (result.fault_detected) continue;
                current = result.current;
                local_move3 = result.local_move;
                local_routing_complete = result.routing_complete;
            }
	    if (local_routing_complete && global_routing_complete){
	        vector<int> okpath = {local_move1,local_move2,local_move3,global_port1,global_port2,weight};
                oklist.push_back(okpath);
                continue;
	    }
        }
    }

    sort(oklist.begin(), oklist.end(), [](const vector<int>& a, const vector<int>& b) {
        if(a[5] != b[5])return a[5] < b[5];
        if(a[0] != b[0])return a[0] > b[0];
	if(a[1] != b[1])return a[1] > b[1];
	return a[2] > b[2];
    });

    if(oklist.size()==0 && (current_node>>Hypercubeport)==(destination_node>>Hypercubeport)){// polarfly escape
    	for (int i = 0; i < (1<<Hypercubeport); i++) {
         for (int j = 0 ; j < (1<<Hypercubeport); j++) {
          for(int k = 0; k < Polarflyport; k++){
	    bool local_routing_complete = false;
	    bool global_routing_complete = false;
	    int current = current_node;
	    const int dest = destination_node;
	    int current_group = current >> Hypercubeport;
	    const int dest_group = dest >> Hypercubeport;
            int hypercube_mv3 = j;
	    int hypercube_mv2 = i;
	    int hypercube_mv1 = hypercube_moves^i^j;
	    if (hypercube_mv2==0)continue;
	    bitset<8> bit1(hypercube_mv1);
	    bitset<8> bit2(hypercube_mv2);
	    bitset<8> bit3(hypercube_mv3);
	    int weight = bit1.count() + bit2.count() + bit3.count();
            cout << "source routing id:" << f->pid << " hypercube hops:" << bitset<8>(hypercube_mv1) << " " << bitset<8>(hypercube_mv2) << " " << bitset<8>(hypercube_mv3) << " polarfly external esc"<< endl;
            local_move1=0;local_move2=0;local_move3=0;
	    global_port1 = 0; global_port2 = 0;
            // local move 1
            if (hypercube_mv1 > 0) {
                auto result = process_local_move(current, dest, hypercube_mv1, in_vc, f->pid);
                if (result.fault_detected) continue;
                current = result.current;
                local_move1 = result.local_move;
            }
            int dest_group_esc = polarfly_connection_table[current_group][j];
	    if(dest_group_esc==current_group)continue; //red group last port
	    //global move 1
            auto result = process_global_move(current, current_group, dest_group_esc, 1, f->pid);
            if (result.global_port == 0) continue;
            current = result.current;
            current_group = result.current_group;
            global_port1 = result.global_port;
            // local move 2
            if (hypercube_mv2 > 0) {
                auto result = process_local_move(current, dest, hypercube_mv2, in_vc+1, f->pid);
                if (result.fault_detected) continue;
                current = result.current;
                local_move2 = result.local_move;
                local_routing_complete = result.routing_complete;
            }
	    //global back 2
            result = process_global_move(current, current_group, dest_group, 2, f->pid);
            if (result.global_port == 0) continue;
            current = result.current;
            current_group = result.current_group;
            global_port2 = result.global_port;
            global_routing_complete = result.routing_complete;
	    if(!global_routing_complete)continue;   
            //local move 3
	    if (hypercube_mv3 > 0) {
                auto result = process_local_move(current, dest, hypercube_mv3, in_vc+2, f->pid);
                if (result.fault_detected) continue;
                current = result.current;
                local_move3 = result.local_move;
                local_routing_complete = result.routing_complete;
            }
            if (local_routing_complete && global_routing_complete){
                vector<int> okpath = {local_move1,local_move2,local_move3,global_port1,global_port2,weight+2};
                oklist.push_back(okpath);
            }
	  }
         }
	}
    }

    if(oklist.size()==0){
      routing_result=1;
      local_move1 = 0;
      local_move2 = 0;
      local_move3 = 0;
      global_port1 = 0;
      global_port2 = 0;
      min_ok_weight = 0;
    } 
    else{
      local_move1 = oklist[0][0];
      local_move2 = oklist[0][1];
      local_move3 = oklist[0][2];
      global_port1 = oklist[0][3];
      global_port2 = oklist[0][4];
      min_ok_weight = oklist[0][5];
    }
    ((int*)f->data)[0] = local_move1;
    ((int*)f->data)[1] = local_move2;
    ((int*)f->data)[2] = local_move3;
    ((int*)f->data)[3] = global_port1;
    ((int*)f->data)[4] = global_port2;
    cout << "source routing id:" << f->id 
         << " src:" << current_node 
         << " dest:" << destination_node 
         << " mv:" << hypercube_moves 
         << " localmv1:" << local_move1 
         << " localmv2:" << local_move2 
         << " localmv3:" << local_move3 
         << " global1:" << global_port1 
         << " global2:" << global_port2;
    if(routing_result==0) {cout << " routing:OK";}
    if(routing_result==1) {cout << " routing:NG";}
    cout << " extrahops:"  << (min_ok_weight-min_weight);
    if(min_weight == min_ok_weight){cout << " minimal";}
    if(min_weight < min_ok_weight){cout << " non-minimal";}
	 cout << endl;
}


void dim_order_polarflyplus( const Router *r, const Flit *f, int in_channel,
                      OutputSet *outputs, bool inject )
{
   outputs->Clear( );
   int out_port=-1;
   int out_vc=0;
   if(inject) {
      outputs->AddRange( -1, 0, 0);
   } else {

    int in_vc=f->vc;
    const int in_port=in_channel;
    const int cur = r->GetID( );
    const int dest = f->dest;
    // ========for adaptive routing======= 
    //const int src_grp = cur>>Hypercubeport;
    //const int dest_grp = dest>>Hypercubeport;
    //int global_port;
    //int local_port;
    //int escape_flag=0;
    if(dest == cur) {
            cout << "routefunc polarfly+ id:" << f->pid << " eject" << endl;
            out_port = 0; // Eject
	    out_vc = in_vc;
    }
    else {
          //==================source routing==================== 
          out_port=-1;
	  if(in_port == 0){ //injection port
             source_routing(f, cur, dest); // write in flits
	  }

          const int local_mv1 = ((int*)f->data)[0];
          const int local_mv2 = ((int*)f->data)[1];
	  const int local_mv3 = ((int*)f->data)[2];
	  const int global_port1 = ((int*)f->data)[3];
	  const int global_port2 = ((int*)f->data)[4];
          int order=0;
          if(in_vc == 0){for(int i=0;i<Hypercubeport;i++){if(in_port==order0[i]+1){order=i;break;}}}
          if(in_vc == 1){for(int i=0;i<Hypercubeport;i++){if(in_port==order1[i]+1){order=i;break;}}}
          if(in_vc == 2){for(int i=0;i<Hypercubeport;i++){if(in_port==order2[i]+1){order=i;break;}}}
          if(in_vc == 3){for(int i=0;i<Hypercubeport;i++){if(in_port==order3[i]+1){order=i;break;}}}
          if(in_vc == 4){for(int i=0;i<Hypercubeport;i++){if(in_port==order4[i]+1){order=i;break;}}}
          if(in_vc == 5){for(int i=0;i<Hypercubeport;i++){if(in_port==order5[i]+1){order=i;break;}}}
          if(in_port==0){
		  order=-1;
                  if ( f->type == Flit::READ_REPLY || f->type == Flit::WRITE_REPLY) {
                     in_vc += VCNUM; //reply inject
                  }	  
	  }

	  if(in_port <= Hypercubeport && order < Hypercubeport) { //Local receive : local move
              for(int i = order+1; i < Hypercubeport; i++){
                if(in_vc == 0){
      		      if((local_mv1 >> order0[i])%2==1){
			out_port=order0[i]+1;
		        out_vc=in_vc;
			cout << "routefunc polarfly+ id:" << f->pid << " local move1 port" << out_port << " vc" << out_vc << endl;
			break;
		      }
		}
		else if (in_vc == 1){
                      if((local_mv2 >> order1[i])%2==1){
                        out_port=order1[i]+1;
                        out_vc=in_vc;
                        cout << "routefunc polarfly+ id:" << f->pid << " local move2 port" << out_port << " vc" << out_vc << endl;
                        break;
                      }
		}
		else if (in_vc == 2){
                      if((local_mv3 >> order2[i])%2==1){
                        out_port=order2[i]+1;
                        out_vc=in_vc;
                        cout << "routefunc polarfly+ id:" << f->pid << " local move3 port" << out_port << " vc" << out_vc << endl;
			break;
                      }
		}
		else if(in_vc == 3){
		      if((local_mv1 >> order3[i])%2==1){
		        out_port=order3[i]+1;
			out_vc=in_vc;
			cout << "routefunc polarfly+ id:" << f->pid << " local move1 port" << out_port << " vc" << out_vc << endl;
     			break;
		      } 
		}
		else if (in_vc == 4){
		      if((local_mv2 >> order4[i])%2==1){
			out_port=order4[i]+1; 
			out_vc=in_vc;
			cout << "routefunc polarfly+ id:" << f->pid << " local move2 port" << out_port << " vc" << out_vc << endl;
     			break;
		      }
		}
		else if (in_vc == 5){
		      if((local_mv3 >> order5[i])%2==1){
			out_port=order5[i]+1;
			out_vc=in_vc;
			cout << "routefunc polarfly+ id:" << f->pid << " local move3 port" << out_port << " vc" << out_vc << endl;
     			break;
		      } 
		}
	      }
	    }
            if((in_port <= Hypercubeport && order == Hypercubeport) || (in_port <= Hypercubeport && out_port == -1)) { 
	    //Local LSB receive or No local mv : local -> global
               if (in_vc == 0 || in_vc == 3){
                 if(global_port1>0){
	           out_port = global_port1;
	           out_vc = in_vc;
                   cout << "routefunc polarfly+ id:" << f->pid << " global move1 port" << out_port << " vc" << out_vc << endl;
		 }
	       }
	       else if (in_vc == 1 || in_vc == 4){
                 if(global_port2>0){
	           out_port = global_port2;
                   out_vc = in_vc;
                   cout << "routefunc polarfly+ id:" << f->pid << " global move2 port" << out_port << " vc" << out_vc << endl;
		 }
	       }
	       if(out_port == -1){//no global move: local -> local
	         for(int i = 0; i < Hypercubeport; i++){
                     if (in_vc == 0){
                        if((local_mv2 >> order1[i])%2==1){
	                  out_port=order1[i]+1;
	                  out_vc=in_vc+1;
	                  cout << "routefunc polarfly+ id:" << f->pid << " local move2 port" << out_port << " vc" << out_vc << " hypercube internal esc" << endl;
			  break;
		       	}
		     }
                     else if (in_vc == 1){
                        if((local_mv3 >> order2[i])%2==1){
			  out_port=order2[i]+1;
			  out_vc=in_vc+1;
			  cout << "routefunc polarfly+ id:" << f->pid << " local move3 port" << out_port << " vc" << out_vc << " hypercube internal esc"<< endl;
                          break;
                        }
		     }
		     else if (in_vc == 3){
           	        if((local_mv2 >> order4[i])%2==1){
			  out_port=order4[i]+1;
		          out_vc=in_vc+1;
	                  cout << "routefunc polarfly+ id:" << f->pid << " local move2 port" << out_port << " vc" << out_vc << " hypercube internal esc" << endl;
	                  break;
			}
		     }
		     else if (in_vc == 4){
			if((local_mv3 >> order5[i])%2==1){
		          out_port=order5[i]+1;
		          out_vc=in_vc+1;
		          cout << "routefunc polarfly+ id:" << f->pid << " local move3 port" << out_port << " vc" << out_vc << " hypercube internal esc" << endl;
                          break;
                        }
                     }
		  }
	       }
   	    }
            if(in_port > Hypercubeport) {  //Global receive: global -> local move 
              for(int i = 0;i < Hypercubeport ; i++){
                if (in_vc == 0){  
      		      if((local_mv2 >> order1[i])%2==1){
                        out_port=order1[i]+1;
                        out_vc=in_vc+1;
                        cout << "routefunc polarfly+ id:" << f->pid << " local move2 port" << out_port << " vc" << out_vc << endl;
                        break;
		      }
		 }
		 else if (in_vc == 1){
                      if((local_mv3 >> order2[i])%2==1){
                        out_port=order2[i]+1;
                        out_vc=in_vc+1;
			cout << "routefunc polarfly+ id:" << f->pid << " local move3 port" << out_port << " vc" << out_vc << endl;
                        break;
                      }
		 }
		 else if (in_vc == 3){
                      if((local_mv2 >> order4[i])%2==1){
                        out_port=order4[i]+1;
                        out_vc=in_vc+1;
                        cout << "routefunc polarfly+ id:" << f->pid << " local move2 port" << out_port << " vc" << out_vc << endl;
                        break;
                      }
                 }
                 else if (in_vc == 4){
                      if((local_mv3 >> order5[i])%2==1){
                        out_port=order5[i]+1;
                        out_vc=in_vc+1;
                        cout << "routefunc polarfly+ id:" << f->pid << " local move3 port" << out_port << " vc" << out_vc << endl;
                        break;
                      }
                 }
              }
	      if (out_port==-1){ //Global receive no hypercube mv: global -> global
                 out_port = global_port2;
                 out_vc = in_vc+1;
                 cout << "routefunc polarfly+ id:" << f->pid << " global move2 port" << out_port << " vc" << out_vc << endl;
	      }
	    }

/* ============adaptive routing=============
	 
            out_port=-1;
	    if(in_port == 0){ //injection port
   	       if ( local_port > -1 ) {
                  if (!r->IsFaultyOutput( local_port )){
                    cout << "routefunc polarfly+ id:" << f->pid << " local move port" << local_port << endl;
                    out_port = local_port;
                  }
		  else{ // Polarfly fault escape
                    if(escape_flag==0){Faultescape++; escape_flag=1;}
                    cout << "#routefunc polarfly+ id:" << f->pid << " global_port" << global_port << " failure" << endl;
                    out_port = polarfly_fault_escape(src_grp,in_port,global_port);
                  }
	       }
	       else{ //No local move: Gloval move vc0
		    if (!r->IsFaultyOutput( global_port )){
                       cout << "routefunc polarfly+ id:" << f->pid << " global move port" << global_port << endl;
                       if (in_port!=global_port) {out_port = global_port;}
                       else{out_port = polarfly_fault_escape(src_grp,in_port,global_port);}
                    }
                    else{ // Polarfly fault escape
                        if(escape_flag==0){Faultescape++; escape_flag=1;}
                        cout << "routefunc polarfly+ id:" << f->pid << " global_port" << global_port << " failure" << endl;
                        out_port = polarfly_fault_escape(src_grp,in_port,global_port);
                    }
	        }
	        out_vc = in_vc;
	   }
	   else if(in_port < Hypercubeport) { //Local receive : local move
	        if ( local_port > -1 ) {
                  if (!r->IsFaultyOutput( local_port )){ 
		    cout << "routefunc polarfly+ id:" << f->pid << " local move port" << local_port << endl;
    	            out_port = local_port;
		    out_vc = in_vc;
		  }	
		  else{ // Hypercube fault escape
		    if(escape_flag==0){Faultescape++; escape_flag=1;}
		    cout << "routefunc polarfly+ id:" << f->pid << " local_port" << local_port << " failure" << endl;
		  }
		}
           }
           if(in_port == Hypercubeport || out_port == -1) { //Local LSB receive or No local move : local -> global move
                  if (!r->IsFaultyOutput( global_port )){ 
 		    cout << "routefunc polarfly+ id:" << f->pid << " local->global move port" << global_port << endl;
	            out_port = global_port; 
		  }
                  else{ // Polarfly fault escape
		    if(escape_flag==0){Faultescape++; escape_flag=1;}
                    cout << "routefunc polarfly+ id:" << f->pid << " global_port" << global_port << " failure" << endl;
                    out_port = polarfly_fault_escape(src_grp,in_port,global_port);
		  }
		  out_vc = in_vc+1;
           }
           if(in_port > Hypercubeport) { 
	        if ( local_port > -1 ) { //Global receive: global -> local move
                  if (!r->IsFaultyOutput( local_port )){ 
                    cout << "routefunc polarfly+ id:" << f->pid << " global->local move port" << local_port << endl;
    	            out_port = local_port;
         	    out_vc = in_vc; 
		  }
		  else{ // Hypercube fault escape
                    if(escape_flag==0){Faultescape++; escape_flag=1;}
	            cout << "routefunc polarfly+ id:" << f->pid << " local_port" << local_port << " failure" << endl;
		  }
		}
		else { //Global receive and No local move: global move + vc increment
	            if (!r->IsFaultyOutput( global_port )){
                       cout << "routefunc polarfly+ id:" << f->pid << " global move vc+1 port" << global_port << endl;
		       if (in_port!=global_port) {out_port = global_port;}
		       else{out_port = polarfly_fault_escape(src_grp,in_port,global_port);}
		    }  
		    else{ // Polarfly fault escape 
			if(escape_flag==0){Faultescape++; escape_flag=1;}
			cout << "routefunc polarfly+ id:" << f->pid << " global_port" << global_port << " failure" << endl;
			out_port = polarfly_fault_escape(src_grp,in_port,global_port);
                    }
		    if (out_port > -1 ) out_vc = in_vc + 1;
		}
           }
           if(out_port == -1){ //Global receive and Hypercube escape: global move + vc increment
                  if (!r->IsFaultyOutput( global_port )){
		    cout << "routefunc polarfly+ id:" << f->pid << " global move vc+1 port" << global_port << endl;
              	    if (in_port!=global_port) {out_port = global_port;}
                    else{out_port = polarfly_fault_escape(src_grp,in_port,global_port);}
		  }
		  else{ // Polarfly fault escape
	            if(escape_flag==0){Faultescape++; escape_flag=1;}
	            cout << "routefunc polarfly+ id:" << f->pid << " global_port" << global_port << " failure" << endl;
		    out_port = polarfly_fault_escape(src_grp,in_port,global_port); 
                  }
		  if (out_port > -1 ) out_vc = in_vc + 1;
           }
*/	   
	   if(out_port==-1){cout << "routefunc polarfly+ router:" << r->FullName() << "#" << r->GetID( ) << " id:" << f->pid << " outporterr" << endl;}
           assert(out_port > -1);
    } 
    
    cout << "routefunc polarfly+ router:" << r->FullName() << "#" << r->GetID( ) << " id:" << f->pid << " src:" << f->src << " dest:" << f->dest << " in_p" << in_channel << "_vc" << f->vc << " out_p" << out_port << "_vc" << out_vc << " type:" << f->type << endl;
   //VC chk
   cout << "routefunc polarfly+ id:" << f->pid << " vc chk" << endl;
   if ( f->type == Flit::READ_REQUEST || f->type == Flit::WRITE_REQUEST) {
       if((cur!=dest) && ((out_vc < 0) || (out_vc > VCNUM-1)))cout << "routefunc polarfly+ id:" << f->pid << " vcerr:" << out_vc << endl; 
       assert((cur==dest) || ((out_vc >= 0) && (out_vc <= VCNUM-1))); //pending or VC0-2
   }
   else if ( f->type == Flit::READ_REPLY || f->type == Flit::WRITE_REPLY) {
       if((cur!=dest) && ((out_vc < VCNUM) || (out_vc > 2*VCNUM-1)))cout << "routefunc polarfly+ id:" << f->pid << " vcerr:" << out_vc << endl;
       assert((cur==dest) || ((out_vc >= VCNUM) && (out_vc <= 2*VCNUM-1))); //pending or VC3-5
   }
   else {
       assert((out_vc >= 0) && (out_vc <= 2*VCNUM-1));
   }
    traffic_table[ r->GetID( )][out_port]++;
    //cout << "traffic table node" <<  r->GetID( ) << " port" << out_port << " :" << traffic_table[ r->GetID( )][out_port] << endl;
    outputs->AddRange( out_port, out_vc, out_vc );
   }
}

//=============================================================

void InitializeRoutingMap( const Configuration & config )
{
  Hypercubeport = config.GetInt( "k" );
  Polarflyport = config.GetInt( "n" );
  gNumVCs = config.GetInt( "num_vcs" );

  //
  // traffic class partitions
  //
  gReadReqBeginVC    = config.GetInt("read_request_begin_vc");
  if(gReadReqBeginVC < 0) {
    gReadReqBeginVC = 0;
  }
  gReadReqEndVC      = config.GetInt("read_request_end_vc");
  if(gReadReqEndVC < 0) {
    gReadReqEndVC = gNumVCs / 2 - 1;
  }
  gWriteReqBeginVC   = config.GetInt("write_request_begin_vc");
  if(gWriteReqBeginVC < 0) {
    gWriteReqBeginVC = 0;
  }
  gWriteReqEndVC     = config.GetInt("write_request_end_vc");
  if(gWriteReqEndVC < 0) {
    gWriteReqEndVC = gNumVCs / 2 - 1;
  }
  gReadReplyBeginVC  = config.GetInt("read_reply_begin_vc");
  if(gReadReplyBeginVC < 0) {
    gReadReplyBeginVC = gNumVCs / 2;
  }
  gReadReplyEndVC    = config.GetInt("read_reply_end_vc");
  if(gReadReplyEndVC < 0) {
    gReadReplyEndVC = gNumVCs - 1;
  }
  gWriteReplyBeginVC = config.GetInt("write_reply_begin_vc");
  if(gWriteReplyBeginVC < 0) {
    gWriteReplyBeginVC = gNumVCs / 2;
  }
  gWriteReplyEndVC   = config.GetInt("write_reply_end_vc");
  if(gWriteReplyEndVC < 0) {
    gWriteReplyEndVC = gNumVCs - 1;
  }

  /* Register routing functions here */

  gRoutingFunctionMap["dim_order_polarflyplus"]         = &dim_order_polarflyplus;
  // ===================================================
  // Balfour-Schultz
  gRoutingFunctionMap["nca_fattree"]         = &fattree_nca;
  gRoutingFunctionMap["anca_fattree"]        = &fattree_anca;
  gRoutingFunctionMap["nca_qtree"]           = &qtree_nca;
  gRoutingFunctionMap["nca_tree4"]           = &tree4_nca;
  gRoutingFunctionMap["anca_tree4"]          = &tree4_anca;
  gRoutingFunctionMap["dor_mesh"]            = &dim_order_mesh;
  gRoutingFunctionMap["dor_hcube"]            = &dim_order_hcube;
  gRoutingFunctionMap["xy_yx_mesh"]          = &xy_yx_mesh;
  gRoutingFunctionMap["adaptive_xy_yx_mesh"]          = &adaptive_xy_yx_mesh;
  // End Balfour-Schultz
  // ===================================================

  gRoutingFunctionMap["dim_order_mesh"]  = &dim_order_mesh;
  gRoutingFunctionMap["dim_order_hcube"]  = &dim_order_hcube;
  gRoutingFunctionMap["dim_order_ni_mesh"]  = &dim_order_ni_mesh;
  gRoutingFunctionMap["dim_order_pni_mesh"]  = &dim_order_pni_mesh;
  gRoutingFunctionMap["dim_order_torus"] = &dim_order_torus;
  gRoutingFunctionMap["dim_order_ni_torus"] = &dim_order_ni_torus;
  gRoutingFunctionMap["dim_order_bal_torus"] = &dim_order_bal_torus;

  gRoutingFunctionMap["romm_mesh"]       = &romm_mesh; 
  gRoutingFunctionMap["romm_ni_mesh"]    = &romm_ni_mesh;

  gRoutingFunctionMap["min_adapt_mesh"]   = &min_adapt_mesh;
  gRoutingFunctionMap["min_adapt_torus"]  = &min_adapt_torus;

  gRoutingFunctionMap["planar_adapt_mesh"] = &planar_adapt_mesh;

  // FIXME: This is broken.
  //  gRoutingFunctionMap["limited_adapt_mesh"] = &limited_adapt_mesh;

  gRoutingFunctionMap["valiant_mesh"]  = &valiant_mesh;
  gRoutingFunctionMap["valiant_torus"] = &valiant_torus;
  gRoutingFunctionMap["valiant_ni_torus"] = &valiant_ni_torus;

  gRoutingFunctionMap["dest_tag_fly"] = &dest_tag_fly;

  gRoutingFunctionMap["chaos_mesh"]  = &chaos_mesh;
  gRoutingFunctionMap["chaos_torus"] = &chaos_torus;
}
