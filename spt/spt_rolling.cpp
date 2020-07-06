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
*********************************************************************
  top level to run spt netlist

  This is a special version of top lever SPRNT driver. It will *not* be
automatically built when "make" is invoked. If you want to build it, you have to
know how to modify the Makefile. If you don't know how, then probably you don't
want to touch (and use) this driver anyway.

  Purpose:
    This special SPRNT driver will run the given netlist, but will always load
    the initial conditions from a "state" file. After completing the unsteady
solve, it will also save the current state to the same file, overwritten its
content.

    The use scenario of this SPRNT drive is the following:

        - run SPRNT using netlist1 (and store the state to a state file)
        - run SPRNT using netlist2, load the initial conditions from state file,
          after completing, save the states to state again
        - run SPRNT using netlist3, load the initial conditios from state file,
etc

    It's clear that netlist1 and netlist2 should be built for the same river
network, otherwise all hell will break lose. Therefore, netlist1, netlist2,
netlist3,... should be identical, except for the unsteady portion of the forcing
terms, the boundary conditions and the lateral flows are different.

    The specifications of time are assumed to be additive. The forcing terms in
    netlist1 should be specified as [0,T1], the netlist2 [0,T2], and in
    netlist3 [0, T3], which correpsonds to the real time sequence as
    [0, T1, T1+T2, T1+T2+T3]

    The API-mode does the samething, but much cleaner.

*********************************************************************/

#include <ctype.h>
#include <stdio.h>

#include "build_defs.h"
#include "sim_parse.h"
#include "smap.h"
#include "spttimemeas.h"

#include "base.h"
#include "options.h"
#include "subcatch.h"

// OPT needs to be in the global namespace
Options OPT;
Stats STAT;

int read_spt_from_file(FILE *fF, Subcatchment *s, NameStore *NODE_NAMES,
                       SMap *hw, FILE *out);
char *get_basename(char *);

int main(int argc, char **argv) {
  const int buf_len = 512;
  int rc, docheck, print_flag;
  char fname[buf_len], outname[buf_len], statefname[buf_len];
  FILE *F, *FS;
  mytm TM;

  // These are permanent stores and we need to worry about memory management
  Subcatchment *SUB = NULL;
  NameStore *NODE_NAMES;
  SMap HASH;

  if (argc < 3) {
    printf("Usage: %s <netlistfile> <state file name>\n", argv[0]);
    printf("       %s -checkonly <netlistfile>\n", argv[0]);
    return (-1);
  }

  docheck = 0;
  if (argc == 3 && strcmp(argv[1], "-checkonly") == 0) {
    docheck = 1;
    strncpy(fname, argv[2], buf_len);
  } else {
    docheck = 0;
    strncpy(fname, argv[1], buf_len);
    strncpy(statefname, argv[2], buf_len);
  }

  if ((F = fopen(fname, "r")) == NULL) {
    printf("Bummer: unable to open netlist file %s\n", argv[1]);
    return (-1);
  }

  PrintHeader(stdout);

  /* allocate the subcatchment and the topograph */
  SUB = new Subcatchment(0);
  NODE_NAMES = new NameStore();

  // default option values
  OPT.Alpha() = 0.8;
  OPT.Beta() = 1.03;
  OPT.HT() = 1e-3;
  OPT.MinDT() = 0.001;
  OPT.SSTol() = 5e-8;
  OPT.Tol() = 1e-6;
  OPT.EpsilonA() = 1e-4;
  OPT.DebugLevel() = 0;

  rc = read_spt_from_file(F, SUB, NODE_NAMES, &HASH,
                          stdout); // also modifies OPT and
                                   // STAT in the global
                                   // namespace
  // save the netlist file
  STAT.SetInFile(get_basename(fname));

  if (rc < 0)
    goto out; // there is error!

  // do a topo check
  rc = SUB->TopologyCheck();

  if (rc < 0) {
    SUB->TopoPrintErrMsg(stdout, "== tchk ==", NODE_NAMES->Store());
    goto out;
  }

  // quit if docheck == 1
  if (OPT.CheckOnly() == 1 || docheck == 1) {
    printf("The connectivity of the given netlist appears to be correct\n");
    goto out;
  }

  // make solver
  SUB->MakeSolver(SUB->GetNumNodes());

  //  SUB->InitSolutions();
  FS = fopen(statefname, "r");
  rc = ERROR;

  if (FS) {
    const int perf_check = 0; // force to ignore the check. Could be dangerous!
    rc = SUB->LoadSteadyStateFromFile(FS, perf_check);
    SUB->CopyXpToXtm1();
    SUB->CopyXpToX();

    if (rc == OK)
      printf("[II]: Loaded states from \"%s\"\n", statefname);
    fclose(FS);
  }

  // steady solve, two phases
  if (rc != OK) { // either no state file, or state file is not good
    rc = SUB->SteadySolve(600, 600, 2 * OPT.Tol());
    if (rc < 0) {
      fprintf(stdout, "[EE]: first phase of steady-state solve failed\n");
      goto out;
    }

    rc = SUB->SteadySolve(25.0, 45, 45, OPT.Tol());
    if (rc < 0) {
      fprintf(stdout, "[EE]: second phase of steady-state solve failed\n");
      goto out;
    }
  }

  // store back if needed
  if (rc > 0 && STAT.SSFile()) {
    FS = fopen(STAT.SSFile(), "w");
    SUB->SaveSteadyStateToFile(FS);
    fclose(FS);
    fprintf(stdout, "[II]: Saved steady state to file \"%s\"\n", STAT.SSFile());
  }

  // run unsteady
  print_flag = 0;
  if (OPT.PrintQ() == 1)
    print_flag |= PRT_Q;
  if (OPT.PrintA() == 1)
    print_flag |= PRT_A;
  if (OPT.PrintZ() == 1)
    print_flag |= PRT_Z;
  if (OPT.PrintD() == 1)
    print_flag |= PRT_D;
  sprintf(outname, "%s.output.dat", fname);
  if (print_flag) {
    FS = fopen(outname, "w");
    if (FS == NULL) {
      fprintf(stdout, "[II]: Unable to open output file %s. Request ignored\n",
              outname);
    } else {
      fprintf(stdout, "[II]: Unsteady results will be stored in \"%s\"\n",
              outname);
    }
    if (OPT.PrintXY() == 1)
      print_flag |= PRT_XY; // XY printing is only turned on when
                            // others are enabled
  } else {
    fprintf(stdout, "[II]: No printing request made. Nothing will be stored\n");
    FS = NULL;
  }

  TM.start();
  rc = SUB->UnsteadySolve(OPT.StopTime(), 1, 25, OPT.Tol(), NULL, FS,
                          print_flag, OPT.PrintStart(), NODE_NAMES->Store());
  TM.stop();
  if (FS)
    fclose(FS);

  fprintf(stdout, "[II]: Simulation of the %d-min event took %.3f seconds.\n",
          (int)(OPT.StopTime() / 60.0), TM.read());

  // store the states
  FS = fopen(statefname, "w");
  if (!FS) {
    printf("Bummer: unable to open file \"%s\" to write state. Check "
           "permissions!\n",
           statefname);
  } else {
    SUB->SaveSteadyStateToFile(FS);
    fprintf(stdout, "[II]: Saved final states to file \"%s\"\n", statefname);
  }
  if (FS)
    fclose(FS);

out:
  if (F)
    fclose(F);
  if (SUB)
    delete SUB;
  if (NODE_NAMES)
    delete NODE_NAMES;

  return 0;
}

/* end */
// Local variables:
// mode: c++
// End:
