/*********************************************************************
#
#  Copyright (C) 2011, 2014 International Business Machines
#
#  Author:  Frank Liu, IBM
#
#  All Rights Reserved. This program and the accompanying materials
#  are made available under the terms of the Eclipse Public License v1.0
#  which accompanies this distribution, and is available at
#  http://www.eclipse.org/legal/epl-v10.html
#
*********************************************************************/

/*
  An example to demonstrate how to run SPRINT through API
*/

#include <stdio.h>

#include "base.h"
#include "waveforms.h"
#include "options.h"
#include "subcatch.h"

/* we need to initialize two global objects */
Options      OPT;
Stats        STAT;

int main() {
  Subcatchment *BASIN;  // each "subcatchment" forms a set of Saint Venant equations

  // create a Subcatchment object to hold river network
  BASIN = new Subcatchment(0); // index 0 is arbitrary
  
  // initialize graph and indicating that there will be less than 128 nodes, less than 128 edges, and
  // less than 32 junctions
  BASIN->InitGraph(128,128,32);

  // all units are assumed to be in SI. So you should convert them 
  // each node should have a unique index, which should be an unsigned interger
  // starting from 0

  // node 0: s0=0.0083, n=0.04, trapezoidal width=0.1, side slope=1.5  (required)
  // x-coord=0, y-coord=0, inital q=0 initial a=0 (these are optional)
  BASIN->MakeNode(0, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);

  // node 1, same specs (doesn't have to be. simply because I am lazy)
  BASIN->MakeNode(1, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);
  BASIN->MakeNode(2, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5); // node 2
  BASIN->MakeNode(3, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5); // 3
  BASIN->MakeNode(4, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5); // etc etc
  BASIN->MakeNode(5, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);
  BASIN->MakeNode(6, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);
  BASIN->MakeNode(7, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);
  BASIN->MakeNode(8, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);
  BASIN->MakeNode(9, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);
  BASIN->MakeNode(10, 0.0083, 0.04, 0, 0, 0, 0, 702.23, 0, TRAP, 0.1,1.5);

  // need to send a signal that the node assignments are done!
  // after this point, no new nodes can be added.
  BASIN->AssignNodes();

  /* we need to define a bunch of piece-wise-linear time series, for Q, lateral and A */
  double time[] = {0, 60, 120, 180, 240, 300, 360,420, 480};
  double junk1[] = {5, 5, 6, 7, 8, 7, 7, 6, 6};
  // we ask subcatchment to create a source for us, but we hold the returned index key because we
  // need to use it later
  source_id idx1 = BASIN->MakeSource(9, &time[0], &junk1[0]);
  
  // ditto
  double junk2[] = {4, 4, 4, 5, 3, 3, 3, 2, 2};
  source_id idx2 = BASIN->MakeSource(9, &time[0], &junk2[0]);

  double junk3[] = {8, 8, 8, 8.1, 8.2, 8.3, 8.4, 8.4, 8.3};
  source_id idx3 = BASIN->MakeSource(9, &time[0], &junk3[0]);

  double junk4[] = {0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.1,0}; 
  source_id idx4 = BASIN->MakeSource(9, &time[0], &junk4[0]);

  double junk5[] = {0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1,0.1}; 
  source_id idx5 = BASIN->MakeSource(9, &time[0], &junk5[0]);
  
  /* add edges to connect the nodes */
  // segment connection node 0 and 1, where 0 is the upstream node, at the length of 40m
  BASIN->MakeStvEquation(0, 1, 40);

  BASIN->MakeStvEquation(1, 2, 50); // 1 to 2
  BASIN->MakeStvEquation(2, 3, 60, idx4); // 2 to 3, also has a lateral source,
  // now you see the reason why we have asked for the source index keys

  BASIN->MakeStvEquation(4, 5, 80); // 4 to 5
  BASIN->MakeStvEquation(5, 6, 80); // 5 to 6

  BASIN->MakeStvEquation(7, 8, 80); 
  BASIN->MakeStvEquation(8, 9, 80, idx5); // with a lateral flow
  BASIN->MakeStvEquation(9, 10, 80); 

  /* we need to define a junction, node 3, 6, flow into node 7 */
  // assuming the flow mix is 50%-50% 
  int in_nodes[]={3, 6};
  double r[]={0.5,0.5};
  BASIN->MakeDepEquation(7, 2, &in_nodes[0], &r[0]);

  /* add the upstream Q sources at node 0 */
  // again use the source index key 
  BASIN->MakeQrSrcEquation(0, idx1 );

  // at node 4
  BASIN->MakeQrSrcEquation(4, idx2);

  // at node 10, this is the downstream boundary condition, 
  // but it's time varying so we define it as a "source"
  BASIN->MakeArSrcEquation(10, idx3 );

  // First we need to do a topology check to see if the network is properly connected:
  int rc = BASIN->TopologyCheck();
  if ( rc<0 ) { 
    // something is not right, let's see the diagnosis
    BASIN->TopoPrintErrMsg(stdout,"-- tchk error --", NULL);
    return (-1); // there is no need to proceed any further if the topology is wrong
  }

  /* 
     At this point, the stream network has been properly defined and checked. We are ready
     to start the simulation 
  */

  // Before going to the next step, we need to make a solver
  BASIN->MakeSolver( BASIN->GetNumNodes() );

  /* we should first compute the steady state solutions */
  // two step steady state
  // the first step is a hybrid method with graph+simplified SVE solve
  int rc2 = BASIN->SteadySolve(120, 120, 5*OPT.Tol());  
  // the second step is a pseudo-unsteady approach
  int rc3 = BASIN->SteadySolve(25.0, 5, 45, OPT.Tol());
  if ( rc2 < 0 || rc3 < 0 ) {
    // any of the two return values less than 0 means things are not right
    // since the topology of the stream net work is correct, most likely the boundary
    // conditions are out of whack, or some nodes have bizzar cross sections
    printf("hmmm, steady state solve couldn't be completed. Check your boundary conditions\n");
    return (-2);
  }
  
  /*
    Now we are ready to compute unsteady solutions
  */

  /* We need to allocate an object of waveform class to hold the unsteady results */
  Waveforms *WV = new Waveforms;
  WV->Init( BASIN->GetNumNodes() );
  
  // here we can also change some options, which controls the behavior
  // see options.h for details
  
  // the variables we are going to store in the waveform object
  int what2save = PRT_Q|PRT_A|PRT_D|PRT_Z; // everything

  // proceed to unsteady simulation 
  // the time interval is 0 to 150, store results in the waveform object WV
  BASIN->UnsteadySolve(150, 1, 25, OPT.Tol(), WV, NULL, what2save, 0, NULL);

  // temporary space for printing
  double rep_time[1024]; 
  double rep_val1[1024], rep_val2[1024];
  int    len, jj;

  // query what's inside the waveform object
  len = WV->Length();
  printf("Found %d records in the waveform container\n", len);
  WV->TimeSeq(&rep_time[0]);
  WV->GetQatNode(6, &rep_val1[0]); // get Q and Z at node 6 and print
  WV->GetZatNode(6, &rep_val2[0]);
  printf("Q and Z at node 6:\n");
  for (jj=0; jj<len; jj++) {
    printf("t = %8.3f  Q = %12.5e,  Z = %12.5e\n", rep_time[jj], rep_val1[jj], rep_val2[jj]);
  }

  /* 
     here we can do aother bunch of stuff. Example:
         further processing the results in waveform object;
	 invoke another model, e.g., run-off, 

	 or simply go to get a cup of joe ......
	 
   */

  // if you want to save the previous results, make a copy
  Waveforms NWV( *WV );

  len = NWV.Length();
  printf("Found %d records in the copy of the waveform container\n", len);
  NWV.TimeSeq(&rep_time[0]);
  NWV.GetQatNode(6, &rep_val1[0]); // again node 6, printing, should be the same 
  NWV.GetZatNode(6, &rep_val2[0]);
  printf("Q and Z at node 6:\n");
  for (jj=0; jj<len; jj++) {
    printf("t = %8.3f  Q = %12.5e,  Z = %12.5e\n", rep_time[jj], rep_val1[jj], rep_val2[jj]);
  }

  // now we are ready to proceed forward, this time to the time point 420 sec (from 150 sec)
  // note that the existing content in the waveform object will be overwritten
  BASIN->UnsteadySolve(420, 1, 25, OPT.Tol(), WV, NULL, what2save, 0, NULL);

  // let's see what we got
  len = WV->Length();
  printf("Found %d records in the waveform container\n", len);
  WV->TimeSeq(&rep_time[0]);
  WV->GetQatNode(6, &rep_val1[0]); // again node 6, printing
  WV->GetZatNode(6, &rep_val2[0]);
  printf("Q and Z at node 6:\n");
  for (jj=0; jj<len; jj++) {
    printf("t = %8.3f  Q = %12.5e,  Z = %12.5e\n", rep_time[jj], rep_val1[jj], rep_val2[jj]);
  }


  /*

    Again we can do other things here. For example, we just received the run-off from a
    hydrological model for the next few minutes. So we are going to "update" the Q sources
    before proceeding

  */

  // we have new data from 520 onward for two upstream Q sources, so we update them
  // the other three will be kept at the value of the last time point (i.e., 480 sec)

  double new_time[] = {520, 560, 600, 620, 650, 720};
  
  double junk9[] = {6, 6, 5, 5, 6, 7};
  int rr = BASIN->UpdateSource(idx1, 6, &new_time[0], &junk9[0]);

  printf("Update return code: %d\n", rr);  // this return code will tell us whether we succeeded, -1 means bad

  double junk8[] = {2, 3, 3, 4, 5, 4};
  BASIN->UpdateSource(idx2, 6, &new_time[0], &junk8[0]); // should check the return flag,
							 // but I am lazy ...
  
  // perform simulation to the 680 sec time point
  BASIN->UnsteadySolve(680, 1, 25, OPT.Tol(), WV, NULL, what2save, 0, NULL);

  // check ...
  len = WV->Length();
  printf("Found %d records in the waveform container\n", len);
  WV->TimeSeq(&rep_time[0]);
  WV->GetQatNode(6, &rep_val1[0]);  // again check node 6, but can be any valid node
  WV->GetZatNode(6, &rep_val2[0]);
  printf("Q and Z at node 6:\n");
  for (jj=0; jj<len; jj++) {
    printf("t = %8.3f  Q = %12.5e,  Z = %12.5e\n", rep_time[jj], rep_val1[jj], rep_val2[jj]);
  }
  
  /*
    etc, etc
    
  */

  // before exiting, we need to do some house cleaning
  if (WV) delete WV;
  if (BASIN) delete BASIN;

  printf("Bye!\n");

  return 0;
}

// Local Variables:
// mode: c++
// End:
