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

// parser to read the spt netlist from a file and build the interal
// data structure + options + topograph

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "base.h"
#include "subcatch.h"
#include "options.h"
#include "sim_parse.h"
#include "smap.h"
#include "sptchksum.h"

/* use to tokenize the string */
#define P_DELIM ";, \t\n\r="
#define DES_SIZE(A) (sizeof((A))/sizeof(Descriptor))

Descriptor Spt_Descriptor[]  = {
  { "OPTIONS", no, ft_block },
  { "NODE",    yes, ft_block },
  { "SEGMENT", yes, ft_block },
  { "QSOURCE", yes, ft_block },
  { "BOUNDARYCONDITION", yes, ft_block},
  { "LATERALSOURCE", no, ft_block },
  { "JUNCTION",      no, ft_block }
};

Descriptor Node_Descriptor[] = {
  { "ID", yes, ft_ascii },
  { "SR", yes, ft_real  },
  { "N",  yes, ft_real  },
  { "ZR", yes, ft_real  },
  { "HR", yes, ft_real  },
  { "XCOORD", no, ft_real },
  { "YCOORD", no, ft_real },
  { "TRAPEZOIDAL", secondary, ft_block },
  { "RECTANGULAR", secondary, ft_block },
  { "XY"         , secondary, ft_block },
  { "INTRINSIC",   secondary, ft_block }
};

Descriptor Time_Descriptor[] = {
  { "SECOND", no, ft_real },
  { "MINUTE", no, ft_real },
  { "HOUR",   no, ft_real }
};

Descriptor Trap_Descriptor[] = {
  { "BOTTOMWIDTH", yes, ft_real },
  { "SLOPE",       yes, ft_real }
};

Descriptor Rect_Descriptor[] = {
  { "BOTTOMWIDTH", yes, ft_real }
};

Descriptor XY_Descriptor[] = {
  { "X", yes, ft_real },
  { "Y", yes, ft_real },
  { "DUMMY", no, ft_real} // we need this because it has be like TS, with three fields
};

Descriptor TS_Descriptor[] = {
  { "T", yes, ft_real },
  { "V", yes, ft_real },
  { "TIMEUNIT", no, ft_real}
};

Descriptor APYW_Descriptor[] = {
  { "A", yes, ft_real },
  { "P", yes, ft_real },
  { "Y", yes, ft_real },
  { "W", yes, ft_real }  // need quad to store the data
};


Descriptor Segment_Descriptor[] = {
  { "UP",      yes, ft_ascii },
  { "DOWN",    yes, ft_ascii  },
  { "LENGTH",  yes, ft_real  }
};

Descriptor Junc_Descriptor[] = {
  { "UP1",     yes, ft_ascii },
  { "UP2",     yes, ft_ascii },
  { "DOWN",    yes, ft_ascii },
  { "COEFF1",  yes, ft_real},
  { "COEFF2",  yes, ft_real}
};

Descriptor Qsource_Descriptor[] = {
  { "LOCATION", yes,   ft_ascii },
  { "TIMESERIES", yes, ft_block}
};

Descriptor Latsource_Descriptor[] = {
  { "LOCATION",   yes, ft_ascii},
  { "TIMESERIES", yes, ft_block}
};

Descriptor Bdn_Descriptor[] = {
  { "LOCATION",   yes, ft_ascii},
  { "TYPE",       yes, ft_ascii},
  { "TIMESERIES", yes, ft_block}
};

Descriptor BdnType_Descriptor[] = {
  { "AREA",   yes, ft_ascii},
  { "DEPTH",  yes, ft_ascii}
};

Descriptor Op_Descriptor[] = {
  { "METRIC",     no,  ft_ascii }, // 0
  { "TIMESTEP",    no,  ft_real  },
  { "TIMESTEPUNIT",no,  ft_ascii }, // 2
  { "STOPTIME",    yes, ft_real  },
  { "STOPTIMEUNIT",yes, ft_ascii }, // 4
  { "PRTINTERVAL", no,  ft_ascii }, 
  { "PRTINTERVALUNIT", no, ft_ascii }, // 6
  { "PRTSTART", no,  ft_ascii },    
  { "PRTSTARTUNIT", no, ft_ascii }, // 8
  { "CHECKONLY",   no,  ft_ascii }, 
  { "SSFILE",      no,  ft_ascii }, // 10
  { "LMAX",        no,  ft_ascii }, 
  { "LMIN",        no,  ft_ascii }, // 12
  { "PRTDEPTH",    no,  ft_ascii }, 
  { "PRTSURFELEV", no,  ft_ascii }, // 14
  { "PRTQ",        no,  ft_ascii }, 
  { "PRTA",        no,  ft_ascii }, // 16
  { "PRTCOORD",    no,  ft_ascii }, 
  { "VERBOSE",     no,  ft_ascii }, // 18
  { "EPOCH",       no,  ft_ascii }, // 19
  { "SPINUPTIME", no,  ft_real  }  // 20
};


/* return 1 if found anything
   return 0 if found nothing
   return -1 if error 
*/
int cmp_string(char *t1, const char *t2) {
  int rc = 1;
  while (1) {
    if (*t2==0) break;
    rc *= (toupper(*t1++) == *t2++ ? 1 : 0); /* logic AND, need all characters */
  }
  return (rc==1);
}

FirstDef decipher_keyword(char *pt, Descriptor *des, int sz) {
  FirstDef rc = def_not_found;

  for (int jj=0; jj<sz; jj++) {
    if ( cmp_string(pt, des[jj]._key)) {
      rc = (FirstDef)jj;
      break;
    }
  }
  return rc;
}

FirstDef search_block(FILE *fp, int *cur_ln, Descriptor *des, int sz, StrBuffer *p, StrBuffer *s, FILE *out) {
  int state;
  FirstDef found_key=def_not_found;
  int rc_end, rc_def;
  int ln_number=*cur_ln;

  char tmpbuf[1024];
  char *pt;
  const char def[]="DEF";
  const char end[]="END";

  /* simple state 
     state = 0; nothing, or get ready to return
     state = 1; found first def, everything goes into first buffer
     state = 2; found 2nd def, everthing goes into 2nd buffer
  */
  state = 0;
  found_key = def_not_found;

  p->Reset();
  s->Reset();

  while ( (fgets(tmpbuf, MAX_LINE_LENGTH, fp)) != NULL) {  
    ln_number++;
    if ( tmpbuf[0]=='#'||tmpbuf[0]=='*'||tmpbuf[0]=='%') continue; /* comments */

    pt = strtok(tmpbuf, P_DELIM);
    while (pt!=NULL) {

      rc_end = cmp_string(pt, end);
      rc_def = cmp_string(pt, def);

      if (rc_end == 1) {  /* found "end" */
	state--;
	if ( state < 0 ) {
	  if (out) fprintf(out,"Bummer: redundent END statement on line %d\n",
			   ln_number);
	  return def_error;
	} else if ( state == 0 ) {

	  if (found_key < def_not_found) {
	    p->Finish();
	    s->Finish();
	    *cur_ln = ln_number;
	    return(found_key);
	  }
	}

      } else if ( rc_def == 1 ) { /* found "def" */
	state++;
	if ( state == 1 )      p->Ready(ln_number);
	else if ( state == 2 ) s->Ready(ln_number);
	else if ( state > 2 ) {
	  if (out) fprintf(out,"Bummer: multiple defined DEF on line %d, did you forget an END statement?\n",
			   ln_number);
	  return def_error;
	}

      } else {   /* found something else */

	if ( state ==1 && found_key==def_not_found ) {
	  found_key = decipher_keyword(pt, des, sz);
	}

	if ( found_key < def_not_found ) {
	  if (state==1)          p->Copy_add_Space(pt);
	  else if ( state == 2 ) s->Copy_add_Space(pt);
	} else { // we found something we don't recognize
	  if (out) fprintf(out,"Bummer: Don't know how to deciper DEF \"%s\" on line %d\n",
			   pt, ln_number);
	  return (def_error);
	}
      }

      /* get the next token */
      pt=strtok(NULL, P_DELIM);
    }
  }

  // at this point, state should be zero!
  if ( state != 0 ) {
    if (out) fprintf(out,"Bummer: unbalanced DEF and END statements\n");
    found_key = def_error;
  }

  *cur_ln = ln_number;
  return (found_key);
}

int check_statements(StrBuffer *p, StrBuffer *s, Descriptor *des, int sz, FILE *stdout) {
  int rc=0;

  rc = p->Separate(NULL);

  return rc;
}

inline int is_scientific(const char *a) {
  while ( *a) {
    if ( (*a>='0' && *a<='9')||*a=='+'||*a=='-'||*a=='.'||*a=='E'||*a=='e') a++;
    else return 0;
  }
  return 1;
}

/* return 1, 60, or 3600, -1 if error */
double convert_time_scale(char *t, Descriptor *des) {
  double rc = -1.0;

  if ( cmp_string(t, des[0]._key) ) {
    rc = 1.0;
  } else if ( cmp_string(t, des[1]._key) ) {
    rc = 60.0;
  } else if ( cmp_string(t, des[2]._key) ) {
    rc = 3600.0;
  } 
  return rc;
}

// return 0 or -9 if something is not right
// hard coded assumption that the data would be 2D
int read_2d_arrays(StrBuffer *s, Descriptor *des, SimpleDblArray *d, int with_unit, FILE *out) {
  int jj=0;
  int cntr=0;
  char *t;
  double ss;
  if (with_unit) {
    if ( (t=s->Get_ith_Value(des[2]._key, 0 )) == NULL ) {
      if (out) fprintf(out,"Bummer: unable to locate keyword \"%s\" for array on line %d\n", 
		       des[2]._key, s->LineNumber());
      return 0;
    }

    ss = convert_time_scale(t, Time_Descriptor);
    if ( ss < 0 ) {
      if (out) fprintf(out,"Bummer: only accept time scale of SECOND, MINUTE, HOUR for statement on line %d \n", s->LineNumber());
      return 0;
    }

    d->SetScale( ss );
    jj++;  // move to the first one
  }

  for (; jj<s->NumPairs(); jj+=2) {  // jj has been set
    if ( (t=s->Get_ith_Value(des[0]._key, jj))==NULL) {
      if (out) fprintf(out,"Bummer: unable to locate keyword \"%s\" for array on line %d\n", 
		       des[0]._key, s->LineNumber());
      return 0;
    }
    if ( !is_scientific(t) ) {
      if (out) fprintf(out, "Bummer: invalid number found in field \"%s\" on line %d\n",
		       des[0]._key, s->LineNumber());
      return 0;
    } else {
      d->ithVal(0, cntr) = atof(t);
    }

    if ( (t = s->Get_ith_Value(des[1]._key, jj+1))==NULL) {
      if (out) fprintf(out,"Bummer: unable to locate keyword \"%s\" for array on line %d\n", 
		       des[1]._key, s->LineNumber());
      return 0;
    }
    if ( !is_scientific(t) ) {
      if (out) fprintf(out, "Bummer: invalid number found in field \"%s\" on line %d\n",
		       des[1]._key, s->LineNumber());
      return 0;
    } else {
      d->ithVal(1,cntr++) = atof(t);
    }
  }
  return cntr;
}

// hard coded assumption that the data would be 4D, no checking of unit
int read_4d_arrays(StrBuffer *s, Descriptor *des, SimpleDblQuad *d, FILE *out) {
  int jj=0;
  int cntr=0;
  char *t;
  int idx = 0;

  for (; jj<s->NumPairs(); jj++) {  // jj has been set
    if ( (t=s->Get_ith_Value(des[idx]._key, jj))==NULL) {
      if (out) fprintf(out,"Bummer: unable to locate keyword \"%s\" for array on line %d\n", 
		       des[idx]._key, s->LineNumber());
      return 0;
    }
    if ( !is_scientific(t) ) {
      if (out) fprintf(out, "Bummer: invalid number found in field \"%s\" on line %d\n",
			 des[idx]._key, s->LineNumber());
      return 0;
    } else {
      d->ithVal(idx, cntr) = atof(t);
    }

    // everything here is hard coded for quadruple
    idx++;
    if ( idx>3 ) { // make idx cycle from 0 to 3
      idx = 0;
      cntr++;    // increment cntr when we cross
    }
  }

  return cntr;
}

/* check the sanity of the APYW data */
int check_4d_apyw(int num, SimpleDblQuad *d, FILE *out) {
  double *tmp, threshold;

  // check A, it has to be monotonically increasing
  tmp = d->getIthArray(0); 
  threshold = (tmp[num-1]-tmp[0]) < OPT.Epsilon() ? OPT.Epsilon() : (tmp[num-1]-tmp[0])*OPT.Epsilon();
  for (int jj=0; jj< num-1; jj++) {
    if ( (tmp[jj+1]-tmp[jj])<threshold ) {
      if (out) fprintf(out, "Bummer: wetted area specification in INTRINSIC x-section has to be monotonically increasing.\n");
      return (-1);
    }
  }

  // check P it has to be monotonically increasing
  tmp = d->getIthArray(1); 
  threshold = (tmp[num-1]-tmp[0]) < OPT.Epsilon() ? OPT.Epsilon() : (tmp[num-1]-tmp[0])*OPT.Epsilon();
  for (int jj=0; jj< num-1; jj++) {
    if ( (tmp[jj+1]-tmp[jj])<threshold ) {
      if (out) fprintf(out, "Bummer: wetted perimeter specification in INTRINSIC x-section has to be monotonically increasing.\n");
      return (-1);
    }
  }

  // check Y it has to be monotonically increasing
  tmp = d->getIthArray(2); 
  threshold = (tmp[num-1]-tmp[0]) < OPT.Epsilon() ? OPT.Epsilon() : (tmp[num-1]-tmp[0])*OPT.Epsilon();
  for (int jj=0; jj< num-1; jj++) {
    if ( (tmp[jj+1]-tmp[jj])<threshold ) {
      if (out) fprintf(out, "Bummer: depth specification in INTRINSIC x-section has to be monotonically increasing.\n");
      return (-1);
    }
  }

  // check W it has to be monotonically increasing, but warning only
  tmp = d->getIthArray(2); 
  threshold = (tmp[num-1]-tmp[0]) < OPT.Epsilon() ? OPT.Epsilon() : (tmp[num-1]-tmp[0])*OPT.Epsilon();
  for (int jj=0; jj< num-1; jj++) {
    if ( (tmp[jj+1]-tmp[jj])<threshold ) {
      if (out) fprintf(out, "Warning: surface width specification in INTRINSIC x-section is NOT monotonically increasing.\n");
    }
  }

  return 0;

}

// top level function to parse the spt netlist and construct the SUBCATCHMENT object
// as well as modifying the Options and Sgraph objects
int read_spt_from_file(FILE* F, Subcatchment *SUB, NameStore *NODE_NAMES, SMap *HASH, FILE* out) {

  FirstDef rc;
  int r_code;
  int line_num, good, node_index;
  source_id src_id;
  unsigned int n_op,n_node, n_seg, n_junc, n_qs, n_ls, n_bdn;
  StrBuffer *Prime, *Second;
  SimpleDblArray dba;
  SimpleDblQuad dqa4; 
  double *XX=NULL, *YY=NULL, *PP=NULL, *WW=NULL, *AA=NULL; // pointers only

  FlexVec<int, SHORT_ARRY_LENGTH>  qsources; // we store index only
  FlexVec<int, SHORT_ARRY_LENGTH>  lsources;
  FlexVec<int, SHORT_ARRY_LENGTH>  asources;
  
  ChkSum   mychk;
  char     chktmp[10];

  Prime = new StrBuffer(1023);    // both arrayes can automatically grow if needed
  Second = new StrBuffer(4095,4); // by default set to 4 in case we encounter APYW cases
  r_code = -1;      // by default, the return code indicates failure
  
  /* we first read the option blocks */
  rc = def_not_found;
  line_num = 0;
  
  n_op = n_node = n_seg = n_junc = n_qs = n_ls = n_bdn = 0;
  while (1) { // first do the tally, and a rough check
    rc = search_block(F, &line_num, Spt_Descriptor, DES_SIZE(Spt_Descriptor), Prime, Second, out); 
    switch (rc) {
    case p_options: n_op++; break;
    case p_node:  n_node++; break;
    case p_segment: n_seg++; break;
    case p_junction: n_junc++; break;
    case p_qsource: n_qs++; break;
    case p_lateralsource: n_ls++; break;
    case p_boundarycondition: n_bdn++; break;
    case def_not_found: break;
    case def_error: 
      goto back;
      break;
    default: break;
    }
    if ( rc == def_not_found ) break;
  }
  rewind(F);
  
  if (out) {
    fprintf(out,"[II]: Found %d nodes, %d segmetns, %d junctions, %d q_sources, %d lateral sources, %d bdn condtions, %d options in the netlist\n", 
	    n_node, n_seg, n_junc, n_qs, n_ls, n_bdn, n_op);
  }
  if ( n_node < 2 ) {
    if (out) fprintf(out, "Bummer: not enough node statements, unable to procced\n");
    goto back;
  } else if ( n_seg < 1 ) {
    if (out) fprintf(out, "Bummer: not enough segment statements, unable to procced\n");
    goto back;
  } else if ( n_qs < 1 ) {
    if (out) fprintf(out, "Bummer: not enough Q source statements, unable to procced\n");
    goto back;
  } else if ( n_bdn < 1 ) {
    if (out) fprintf(out, "Bummer: not enough boundary condition statements, unable to procced\n");
    goto back;
  }
  
  // Copy the stat into the STAT object
  STAT.N_Nodes() = n_node;
  STAT.N_Segs() = n_seg;
  STAT.N_Juncs() = n_junc;
  STAT.N_Qsrcs() = n_qs;
  STAT.N_Lsrcs() = n_ls;
  STAT.N_Bnds() = n_bdn;

  // do the memory allocation stuff here
  NODE_NAMES->Init(n_node);
  SUB->InitGraph(n_node, n_seg, n_junc); //Testing

  // allocate memory to store the boundary conditions
  // by default the value is -1
  qsources.size(n_node+10);
  lsources.size(n_node+10);
  asources.size(n_node+10);

  for (unsigned int jj=0; jj<n_node+10; jj++) qsources[jj]=-1;
  for (unsigned int jj=0; jj<n_node+10; jj++) lsources[jj]=-1;
  for (unsigned int jj=0; jj<n_node+10; jj++) asources[jj]=-1;

  // we now proess options, otherwise we use default values
  if ( n_op > 0 ) {
    char *tmp;
    double stoptimeunit=1.0, timestepunit=1.0, prtintervalunit=1.0;
    double stoptime=-1.0, timestep=0.0, prtinterval=0.0;
    double prtstartunit = 1.0, prtstart = 0.0;
    double spinuptime = 0.0;

    rewind(F);
    rc = def_not_found;
    line_num = 0;
    while (1) {  
      rc = search_block(F, &line_num, Spt_Descriptor, DES_SIZE(Spt_Descriptor), Prime, Second, NULL); 
      
      if (rc == def_not_found ) break;
      if (rc != p_options ) continue;
      
      // read the options, in particular UseMetric!     
      good = Prime->Separate(out);
      if (good<0) {
	if (out) fprintf(out,"Bummer: Syntax error in netlist specification at line %d\n", line_num);
	goto back;
      }
      tmp = Prime->Find_Value( Op_Descriptor[0]._key); // USEMETRIC
      if (tmp) {
	OPT.UseMetric() = atoi(tmp)>0 ? 1 : 0;
	if (out) fprintf(out,"[II]: Metric set to %d\n", OPT.UseMetric());
      }

      // we process unit first so that we don't have to worry later
      tmp = Prime->Find_Value( Op_Descriptor[2]._key); // TIMESTEPUNIT
      if (tmp) {
	timestepunit = convert_time_scale(tmp, Time_Descriptor);
      }
      if ( timestepunit < 0 ) {
	if (out) fprintf(out,"[WW]: only accept time unit of SECOND, MINUTE, HOUR for statement on line %d. Use default=SECOND. \n", 
			 Prime->LineNumber());
	timestepunit = 1.0;
      }
      
      tmp = Prime->Find_Value( Op_Descriptor[1]._key); // TIMESTEP
      if (tmp) timestep = atof(tmp);

      tmp = Prime->Find_Value( Op_Descriptor[4]._key); // STOPTIMEUNIT
      if (tmp) {
	stoptimeunit = convert_time_scale(tmp, Time_Descriptor);
      }
      if ( stoptimeunit < 0 ) {
	if (out) fprintf(out,"Bummer:only accept time unit of SECOND, MINUTE, HOUR for statement on line %d. Use default=SECOND. \n", 
			 Prime->LineNumber());
	stoptimeunit = 1.0;
      }

      tmp = Prime->Find_Value( Op_Descriptor[3]._key); // STOPTIME
      if (tmp) stoptime = atof(tmp);

      // also deal with print interval first
      tmp = Prime->Find_Value( Op_Descriptor[6]._key); // PRTINTERVALUNIT
      if (tmp) {
	prtintervalunit = convert_time_scale(tmp, Time_Descriptor);
      }
      if ( prtintervalunit < 0 ) {
	if (out) fprintf(out,"[WW]: only accept time unit of SECOND, MINUTE, HOUR for statement on line %d. Use default=SECOND. \n", 
			 Prime->LineNumber());
	prtintervalunit = 1.0;
      }

      
      tmp = Prime->Find_Value( Op_Descriptor[8]._key); // PRTSTARTUNIT
      if (tmp) {
	prtstartunit = convert_time_scale(tmp, Time_Descriptor);
      }
      if ( prtstartunit < 0 ) {
	if (out) fprintf(out,"[WW]: only accept time unit of SECOND, MINUTE, HOUR for statement on line %d. Use default=SECOND. \n", 
			 Prime->LineNumber());
	prtstartunit = 1.0;

      }

      tmp = Prime->Find_Value( Op_Descriptor[7]._key); // PRTSTART
      if (tmp) prtstart = atof(tmp);

      tmp = Prime->Find_Value( Op_Descriptor[5]._key); // PRTINTERVAL
      if (tmp) prtinterval = atof(tmp);
      if ( prtinterval < 0 ) {
	if (out) fprintf(out,"[WW]: negative print interval specified on line %d. Overwritten with zero!\n",
			 Prime->LineNumber());
	prtinterval = 0.0;
      }
      
      tmp = Prime->Find_Value( Op_Descriptor[9]._key); // CHECKONLY
      if (tmp) {
	OPT.CheckOnly() = atoi(tmp) > 0 ? 1 : 0;
	if (out) fprintf(out,"[II]: CheckOnly set to %d\n", OPT.CheckOnly());
      }
      
      tmp = Prime->Find_Value( Op_Descriptor[10]._key); // SSFILE
      if (tmp) {
	STAT.SetSSFile(tmp);
	if (out) fprintf(out, "[II]: SSFile (SteadyState file) set to \"%s\"\n", tmp);
      }

      tmp = Prime->Find_Value( Op_Descriptor[11]._key); // Lmax
      if (tmp) {
	OPT.LMax() = atof(tmp);
	if (out) fprintf(out,"[II]: LMax set to %.4e\n", OPT.LMax() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[12]._key); // Lmin
      if (tmp) {
	OPT.LMin() = atof(tmp);
	if (out) fprintf(out,"[II]: LMin set to %.4e\n", OPT.LMin() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[13]._key); // PrintDepth
      if (tmp) {
	OPT.PrintD() = atoi(tmp) > 0 ? 1 : 0;
	if (out) fprintf(out,"[II]: PrtDepth set to %1d\n", OPT.PrintD() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[14]._key); // PrintSurfElev
      if (tmp) {
	OPT.PrintZ() = atoi(tmp) > 0 ? 1 : 0;
	if (out) fprintf(out,"[II]: PrtSurfElev set to %1d\n", OPT.PrintZ() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[15]._key); // PrintQ
      if (tmp) {
	OPT.PrintQ() = atoi(tmp) > 0 ? 1 : 0;
	if (out) fprintf(out,"[II]: PrtQ set to %1d\n", OPT.PrintQ() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[16]._key); // PrintA
      if (tmp) {
	OPT.PrintA() = atoi(tmp) > 0 ? 1 : 0;
	if (out) fprintf(out,"[II]: PrtA set to %1d\n", OPT.PrintA() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[17]._key); // PrintCoord
      if (tmp) {
	OPT.PrintXY() = atoi(tmp) > 0 ? 1 : 0;
	if (out) fprintf(out,"[II]: PrtCoord set to %1d\n", OPT.PrintXY() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[18]._key); // verbose
      if (tmp) {
	OPT.DebugLevel() = atoi(tmp) > 0 ? atoi(tmp) : 0;
	if (out) fprintf(out,"[II]: Verbose set to %1d\n", OPT.DebugLevel() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[19]._key); // epoch
      if (tmp) {
	STAT.SetEpoch(tmp);
	if (out) fprintf(out, "[II]: Epoch is set to %s\n", STAT.Epoch() );
      }

      tmp = Prime->Find_Value( Op_Descriptor[20]._key); // spinup_time
      if (tmp) spinuptime = atof(tmp);
      if ( spinuptime < 0 ) {
	if (out) fprintf(out,"[WW]: negative spin-up time specified on line %d. Use zero instead (disabled)\n",
			 Prime->LineNumber());
	spinuptime = 0.0;
      }
      
    }

    // check if stoptime is specified (by setting default to 0), 
    // do the scaling,
    // reason: we might have multiple option blocks

    if ( timestep < 0 ) {
      if (out) fprintf(out,"[WW]: negative TimeStep specified. Overwritten with zeor!\n");
    } else if (timestep > 0) {
      timestep *= timestepunit;
      OPT.FixedStep() = timestep;
      if (out) fprintf(out,"[II]: TimeStep set to %.3e second.\n", timestep);
    }
    
    if ( prtinterval < 0 ) {
      if (out) fprintf(out,"[WW]: negative PrInterval specified. Overwritten with zero!\n");
    } else if (prtinterval>0) {
      int tt;
      prtinterval *= prtintervalunit;
      tt = (int) round(prtinterval/60.0);
      tt = tt >=0 ? tt : 0;
      OPT.PrintInterval() = tt;
      if (out) fprintf(out,"[II]: PrtInterval set to %d min\n", tt);
    }

    if (prtstart < 0) {
      if (out) fprintf(out,"[WW]: negative PrtStart specified. Overwritten with zero!\n");
    } else if ( prtstart > 0) {
      int tt;
      prtstart *= prtstartunit;
      tt = (int)round(prtstart/60.0);
      tt = tt>=0 ? tt : 0;
      OPT.PrintStart() = tt;
      if (out) fprintf(out, "[II]: PrtStart set to %d min\n", tt);
    }

    if (stoptime > 0 ) {
      stoptime *= stoptimeunit;
      OPT.StopTime() = stoptime;
      if (out) fprintf(out,"[II]: StopTime set to %.3e second\n", stoptime);
    } else {
      if (out) fprintf(out,"Bummer: a positive StopTime is required in order to run the simulation!\n");
      goto back;
    }

    if ( spinuptime > 0 ) {
      int tt = (int)spinuptime;
      const int min_spin = 300;
      if ( tt > min_spin ) {
	OPT.SpinUpTime() = tt;
	if (out) fprintf(out,"[II]: SpinUpTime set to %d second\n", tt);
      } else {
	OPT.SpinUpTime() = 0;
	if (out) fprintf(out,"[WW]: minimal spinup time is %d seconds. Revert back to zero (spinup disabled)\n", min_spin);
      }
    }
  }
  
  // add nodes,  check duplicate
  rewind(F);
  rc = def_not_found;
  line_num = 0;
  while (1) {  
    char node_id[MAX_LINE_LENGTH];
    char *tmp;
    double s0, n, h0, z0, width, slope;
    XsecType xtp;
    int   num_pairs;
    double q0=0, a0=0;  // we will deal with them later
    double x=-1.0, y=-1.0;  // not used 

    rc = search_block(F, &line_num, Spt_Descriptor, DES_SIZE(Spt_Descriptor), Prime, Second, NULL); 
    
    // we don't check on the error flags since we have done it before
    if (rc == def_not_found ) break;   // we added all nodes
    if (rc == def_error ) goto back;  

    width = 0.0;
    slope = 0.0;
    num_pairs = 0;
    switch (rc) {
    case p_node:
      /* intercept and understand node */
      if ( !Second->HasContent()) {
	if (out) fprintf(out,"Bummer: NODE statement on line %d requires a x-section specification\n", Prime->LineNumber() );
	goto back;
      }
      good = Prime->Separate(out);
      good = Second->Separate(out);
      
      if (good <0) {
	if (out) fprintf(out,"Bummer: syntax error in netlist specification in line %d\n", Prime->LineNumber());
	goto back;
      }

      // parse each field one by one
      tmp = Prime->Find_Value( Node_Descriptor[0]._key); // "id"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the NODE statement on line %d\n", Node_Descriptor[0]._key, 
			 Prime->LineNumber());
	goto back;
      }
      // the clean node id is now in node_id;
      strncpy(node_id, tmp, MAX_WORD_LENGTH);
      node_index = HASH->Check(node_id);
      if ( node_index != -1 ) {
	if (out) fprintf(out,"Bummer: duplicated NODE definition of node \"%s\" on line %d\n", node_id, Prime->LineNumber() );
	goto back;
      }
      node_index = HASH->Create(node_id);
      NODE_NAMES->Insert(node_id, node_index);
      
      tmp = Prime->Find_Value( Node_Descriptor[1]._key);  // "sR"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the NODE statement on line %d\n", Node_Descriptor[1]._key, 
			 Prime->LineNumber());
	goto back;
      }
      s0 = atof(tmp);
      if ( s0 < 1e-6 ) {
	if (out) fprintf(out, "[WW]: reference slope sR at location %s is nonpositive with value of %.4e. \n", node_id, s0);
#if 0
	goto back;
#endif
      }

      tmp = Prime->Find_Value( Node_Descriptor[2]._key);  // "n"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the NODE statement on line %d\n", Node_Descriptor[2]._key, 
			 Prime->LineNumber());
	goto back;
      }
      n = atof(tmp);
      if ( n < OPT.MinN() ) {
	if (out) fprintf(out,"[WW]: Manning's N at location %s is %.4e, less than threshold value %.4e. Supercritical flow might occur\n",
			 node_id, n, OPT.MinN());
      }
      
      tmp = Prime->Find_Value( Node_Descriptor[3]._key);  // "zR"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the NODE statement on line %d\n", Node_Descriptor[3]._key, 
			 Prime->LineNumber());
	goto back;
      }
      z0 = atof(tmp);
      
      
      tmp = Prime->Find_Value( Node_Descriptor[4]._key);  // "hR"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the NODE statement on line %d\n", Node_Descriptor[4]._key, 
			 Prime->LineNumber());
	goto back;
      }
      h0 = atof(tmp);
      
      tmp = Prime->Find_Value( Node_Descriptor[5]._key); // "xcoord"
      if (tmp) {
	x = atof(tmp);    // this is optional
      }

      tmp = Prime->Find_Value( Node_Descriptor[6]._key); // "ycoord"
      if (tmp) {
	y = atof(tmp);
      }

      if ( Second->Cmp_Head(Node_Descriptor[7]._key)) {   // "Trapezoidal
	tmp = Second->Find_Value(Trap_Descriptor[0]._key);  // "bottomwidth"
	if (!tmp) {
	  if (out) fprintf(out, "Bummer: %s x-section on line %d requires field %s\n", Node_Descriptor[5]._key, 
			   Second->LineNumber(), Trap_Descriptor[0]._key);
	  goto back;
	}
	width = atof(tmp);
	
	tmp = Second->Find_Value(Trap_Descriptor[1]._key);  // "slope"
	if (!tmp) {
	  if (out) fprintf(out, "Bummer: %s x-section on line %d requires field %s\n", Node_Descriptor[5]._key, 
			   Second->LineNumber(), Trap_Descriptor[1]._key);
	  goto back;
	} else {
	  slope = atof(tmp);
	}
	xtp = TRAP;
	
      } else if ( Second->Cmp_Head(Node_Descriptor[8]._key)) {  // "Rectangular"
	tmp = Second->Find_Value(Rect_Descriptor[0]._key);  // "bottomwidth"
	if (!tmp) {
	  if (out) fprintf(out, "Bummer: %s x-section on line %d requires field %s\n", Node_Descriptor[6]._key, 
			   Second->LineNumber(), Rect_Descriptor[0]._key);
	  goto back;
	} else {
	  width = atof(tmp);
	}
	xtp = RECT;
	
      } else if ( Second->Cmp_Head(Node_Descriptor[9]._key)) { // "XY"
	if ( Second->NumPairs() % 2 != 0 ) {
	  if (out) fprintf(out,"Bummer: unbalanced XY pairs for %s statement on line %d\n", 
			   Node_Descriptor[7]._key,Second->LineNumber());
	  goto back;
	  
	} else {  // all good, copy the values
	  dba.Reset();
	  num_pairs = read_2d_arrays(Second, XY_Descriptor, &dba, 0, out);
	  
	  if ( num_pairs <= 0 ) {
	    if (out) fprintf(out,"Bummer: error reading 2D data specified on line %d\n", Second->LineNumber() );
	    goto back;
	  }
	}
	xtp = SPLINE;

      } else if ( Second->Cmp_Head(Node_Descriptor[10]._key)) { // "INTRINSIC"
	if ( Second->NumPairs() % 4 != 0 ) {
	  if (out) fprintf(out,"Bummer: unbalanced APYW quad for %s statement on line %d\n", 
			   Node_Descriptor[10]._key,Second->LineNumber());
	  goto back;
	}
	dqa4.Reset();
	num_pairs = read_4d_arrays(Second, APYW_Descriptor, &dqa4, out);
	if ( num_pairs <= 5 ) {
	  if (out) fprintf(out,"Bummer: error reading 4D data specified in line %d, requires at least 5 points \n", Second->LineNumber() );
	  goto back;
	}
	if ( check_4d_apyw(num_pairs, &dqa4, out) < 0 ) {
	  if (out) fprintf(out,"Bummer: INTRINC x-section data specified near line %d are erroneous \n", Second->LineNumber() );
	  goto back;
	}
	xtp = INTRINSIC;
	
      } else {
	if (out) fprintf(out,"Bummer: NODE statement on line %d requires at least one of the x-section specifications: %s, %s, %s or %s\n",
			 Prime->LineNumber(), Node_Descriptor[7]._key, Node_Descriptor[8]._key, Node_Descriptor[9]._key, Node_Descriptor[10]._key);
	goto back;
      }
      
      // at this point, all variables are ready, create the node
      // use node_idx, s0, n, h0, z0, dba (good)
      // first do scaling based UseMetric()
      if ( xtp == SPLINE ) {
	XX = dba.getIthArray(0);
	YY = dba.getIthArray(1);
      } else if ( xtp == INTRINSIC ) {
	AA = dqa4.getIthArray(0);  // sequence is hardcoded in APYW_Descriptor
	PP = dqa4.getIthArray(1);
	YY = dqa4.getIthArray(2);
	WW = dqa4.getIthArray(3);
      }

      if (OPT.UseMetric() == 0 ) { // scaling in case 
	if (xtp == SPLINE) {
	  for (int jj=0; jj<num_pairs; jj++) {
	    XX[jj] *= OPT.FtoM();
	    YY[jj] *= OPT.FtoM();
	  }
	} else if ( xtp == INTRINSIC) {
	  for (int jj=0; jj<num_pairs; jj++) {
	    AA[jj] *= OPT.F2toM2();  // wetted area
	    YY[jj] *= OPT.FtoM();    // everything else is linear
	    PP[jj] *= OPT.FtoM();
	    WW[jj] *= OPT.FtoM();
	  }
	} else {  // RECT or TRAP, only change the bottom width, slope is dimensionless
	  width *= OPT.FtoM();
	}
	z0 *= OPT.FtoM();
	h0 *= OPT.FtoM();
      } 
      
      if ( xtp == SPLINE ) {
#ifdef DBGSPLINE
	printf("node name %s, id %d\n", node_id, (unsigned int)node_index); 
#endif
	SUB->MakeNode((unsigned int)node_index, s0, n, x, y, q0, a0, z0, h0, SPLINE, num_pairs, XX, YY);
      } else if (xtp == INTRINSIC) {
	SUB->MakeNode((unsigned int)node_index, s0, n, x, y, q0, a0, z0, h0, INTRINSIC, num_pairs, AA, PP, YY, WW);
      } else if (xtp == TRAP) {
	SUB->MakeNode((unsigned int)node_index, s0, n, x, y, q0, a0, z0, h0, TRAP, width, slope);
      } else if (xtp == RECT) {
	SUB->MakeNode((unsigned int)node_index, s0, n, x, y, q0, a0, z0, h0, RECT, width);
      }
      
      break;

    default: break;

    } // end switch
  }   // while loop, find all nodes


  // make SUB happy
  SUB->AssignNodes();


  // records Qsrc, Lsrc and Bdn's
  rewind(F);
  rc = def_not_found;
  line_num = 0;
  while (1) {  
    char node_id[MAX_LINE_LENGTH];
    char *tmp;
    int   num_pairs;
    int loc_index, bnd_is_area;
    Node *node=NULL;
    //    double geo_scale=1.0;

    rc = search_block(F, &line_num, Spt_Descriptor, DES_SIZE(Spt_Descriptor), Prime, Second, NULL); 
    
    // we don't check on the error flags since we have done it before
    if (rc == def_not_found ) break;   // we added all nodes
    if (rc == def_error ) goto back;  

    num_pairs = 0;
    switch (rc) {
    case p_qsource:
      if (!Second->HasContent()) {
	if (out) fprintf(out,"Bummer: QSOURCE statement on line %d requires a time series specification\n", Prime->LineNumber() );
	goto back;
      }
      good = Prime->Separate(out);
      good = Second->Separate(out);
      
      if (good<0) {
	if (out) fprintf(out,"Bummer: syntax error in netlist specification in line %d\n", line_num);
	goto back;
      }

      tmp = Prime->Find_Value( Qsource_Descriptor[0]._key); // "location"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the QSOURCE statement on line %d\n", Qsource_Descriptor[0]._key, 
			 Prime->LineNumber());
	goto back;
      }
      if ( !Second->Cmp_Head( Qsource_Descriptor[1]._key) ) { // "timeseries"
	if (out) fprintf(out,"Bummer: QSOURCE statement on line %d requires time series specification: %s\n",
			 Prime->LineNumber(), Qsource_Descriptor[1]._key);
	goto back;
      }

      strncpy(node_id, tmp, MAX_WORD_LENGTH);
      loc_index = HASH->Check(node_id);
      if ( loc_index == -1 ) {
	if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			 node_id, Prime->LineNumber());
	goto back;
      }
      
      dba.Reset();
      num_pairs = read_2d_arrays(Second, TS_Descriptor, &dba, 1, out);
      if ( num_pairs <= 0 ) {
	if (out) fprintf(out,"Bummer: error reading 2D data specified on line %d\n", Second->LineNumber() );
	goto back;
      }

      // ready to insert Q source, with node and time series in dba (good), remember the
      // dba.Scale() field, which is the time scale
      XX = dba.getIthArray(0);
      YY = dba.getIthArray(1);
      src_id = SUB->MakeSource(num_pairs, &XX[0], &YY[0], dba.Scale(), (OPT.UseMetric()==0?OPT.F3toM3():1.0) );
      if (src_id<0) {
	if (out) fprintf(out,"Bummer: problem with QSource at node \"%s\" on line %d\n",
			 node_id, Prime->LineNumber());
	goto back;
      }
      
      if ( qsources[ loc_index] > -1 ) {
	if (out) fprintf(out, "Bummer: duplicated definitino of Qsource at node \"%s\"\n",
			 node_id);
	goto back;
      }

      qsources[ loc_index ] = src_id; // store the index
      break;

    case p_lateralsource:
      if (!Second->HasContent()) {
	if (out) fprintf(out,"Bummer: LATERALSOURCE statement on line %d requires a time series specification\n", Prime->LineNumber() );
	goto back;
      }
      good = Prime->Separate(out);
      good = Second->Separate(out);
      
      if (good<0) {
	if (out) fprintf(out,"Bummer: problem with LateralSource on line %d\n",
			 Prime->LineNumber());
	goto back;
      }

      tmp = Prime->Find_Value( Latsource_Descriptor[0]._key); // "location"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the LATERALSOURCE statement on line %d\n", Latsource_Descriptor[0]._key, 
			 Prime->LineNumber());
	goto back;
      }
      if ( !Second->Cmp_Head( Latsource_Descriptor[1]._key) ) { // "timeseries"
	if (out) fprintf(out,"Bummer: LATERALSOURCE statement on line %d requires time series specification: %s\n",
			 Prime->LineNumber(), Latsource_Descriptor[1]._key);
	goto back;
      }

      strncpy(node_id, tmp, MAX_WORD_LENGTH);
      loc_index = HASH->Check(node_id);
      if ( loc_index == -1 ) {
	if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			 node_id, Prime->LineNumber());
	goto back;
      }
      
      dba.Reset();
      num_pairs = read_2d_arrays(Second, TS_Descriptor, &dba, 1, out);
      if ( num_pairs <= 0 ) {
	if (out) fprintf(out,"Bummer: error reading 2D data specified on line %d\n", Second->LineNumber() );
	goto back;
      }
      // ready to insert lateral source, with loc_index and time series in dba (good),
      // remember dba.Scale()
      XX = dba.getIthArray(0);
      YY = dba.getIthArray(1);
      src_id = SUB->MakeSource(num_pairs, &XX[0], &YY[0], dba.Scale(), (OPT.UseMetric()==0?OPT.F3toM3():1.0) );
      if (src_id<0) {
	if (out) fprintf(out,"Bummer: problem with LATERALSOURCE at node \"%s\" on line %d\n",
			 node_id, Prime->LineNumber());
	goto back;
      }

      if ( lsources[ loc_index ] > -1 ) {
	if (out) fprintf(out, "Bummer: duplicated definitino of LATERALSOURCE at node \"%s\"\n",
			 node_id);
	goto back;
      }
      lsources[ loc_index ] = src_id;

      break;

    case p_boundarycondition:

      if (!Second->HasContent()) {
	if (out) fprintf(out,"Bummer: BOUNDARYCONDITION statement on line %d requires a time series specification\n", Prime->LineNumber() );
	goto back;
      }
      good = Prime->Separate(out);
      good = Second->Separate(out);
      if (good<0) {
	if (out) fprintf(out,"Bummer: syntax error in BoundaryCondition Statement in line %d\n", Prime->LineNumber());
	goto back;
      }

      tmp = Prime->Find_Value( Bdn_Descriptor[1]._key); // "type" 
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the BOUNDARYCONDITION statement on line %d\n", Bdn_Descriptor[1]._key, 
			 Prime->LineNumber());
	goto back;
      }
      if ( cmp_string(tmp, BdnType_Descriptor[0]._key) ) { // "AREA"
	bnd_is_area = 1;
      } else if (cmp_string(tmp, BdnType_Descriptor[1]._key) ) { // "DEPTH"
	bnd_is_area = 0;
      } else {
	if (out) fprintf(out, "Bummer: BoundaryCondition has to be either \"%s\" or \"%s\" for statement on line %d\n",
			 BdnType_Descriptor[0]._key, 
			 BdnType_Descriptor[1]._key,
			 Prime->LineNumber() );
	goto back;
      }
      
      tmp = Prime->Find_Value( Bdn_Descriptor[0]._key); // "location"
      if (!tmp) {
	if (out) fprintf(out, "Bummer: unable to locate field %s for the BOUNDARYCONDITION statement on line %d\n", Bdn_Descriptor[0]._key, 
			 Prime->LineNumber());
	goto back;
      }

      strncpy(node_id, tmp, MAX_WORD_LENGTH);
      loc_index = HASH->Check(node_id);
      if ( loc_index == -1 ) {
	if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			 node_id, Prime->LineNumber());
	goto back;
      }
      
      if ( !Second->Cmp_Head( Bdn_Descriptor[2]._key) ) { // "timeseries"
	if (out) fprintf(out,"Bummer: BOUNDARYCONDITION statement on line %d requires time series specification: %s\n",
			 Prime->LineNumber(), Bdn_Descriptor[2]._key);
	goto back;
      }

      dba.Reset();
      num_pairs = read_2d_arrays(Second, TS_Descriptor, &dba, 1, out);
      if ( num_pairs <= 0 ) {
	if (out) fprintf(out,"Bummer: error reading 2D data specified on line %d\n", Second->LineNumber() );
	goto back;
      }

      // ready to insert Bdn condition, with loc_index and time series in dba (good),
      // remember dba.Scale()
      XX = dba.getIthArray(0);
      YY = dba.getIthArray(1);
      if ( bnd_is_area == 0 ) {
	node = SUB->GetNode( loc_index );
	if ( !node ) {
	  if (out) fprintf(out, "Bummer: node \"%s\" not yet defined!\n", node_id);
	  goto back;
	}
      }
      if ( bnd_is_area) {
	src_id = SUB->MakeSource(num_pairs, &XX[0], &YY[0], dba.Scale(), (OPT.UseMetric()==0?OPT.F2toM2():1.0) );
      } else {
	FlexVec<double> tmp;
	tmp.size(num_pairs);
	node = SUB->GetNode( loc_index);
	for (int jj=0; jj<num_pairs;jj++) { 
	  tmp[jj]= node->GetAbyDepth( YY[jj]*(OPT.UseMetric()==0?OPT.FtoM():1.0) );
	}
	src_id = SUB->MakeSource(num_pairs, &XX[0], &tmp[0], dba.Scale(), 1.0);
      }

      if (src_id < 0) {
	if (out) fprintf(out, "Bummer: problem with BoundaryCondition definition on line %d\n", Prime->LineNumber() );
	goto back;
      }
      if ( asources[ loc_index] >-1 ) {
	if (out) fprintf(out, "Bummer: duplicated definitino of BOUNDARYCONDITION at node \"%s\"\n",
			 node_id);
	goto back;
      }
      asources[ loc_index ] = src_id;
      break;

    default: break;

    } // end switch
  }   // while loop, found all qsource, lsources and asources

  // process the rest, segments, junctions
  rc = def_not_found;
  rewind(F);
  line_num = 0;
  while (1) {  
    int up_idx, down_idx;
    int up1_idx, up2_idx;
    double length, coef1, coef2;
    char *tmp;
    Node *fp, *tp;

    rc = search_block(F, &line_num, Spt_Descriptor, DES_SIZE(Spt_Descriptor), Prime, Second, NULL); 
    
    if (rc == def_error) goto back;
    if (rc == def_not_found) break;

    switch ( rc ) {
    case p_segment:
      good = Prime->Separate(out);
      if (good<0) {
	if (out) fprintf(out,"Bummer: syntax error in netlist specification at line %d\n", Prime->LineNumber());
	goto back;
      }

      tmp = Prime->Find_Value( Segment_Descriptor[0]._key); // "up" 
      if ( !tmp) {
	if (out) fprintf(out, "Bummer: SEGMENT statement on line %d requires the \"%s\" field.\n",
			 Prime->LineNumber(), Segment_Descriptor[0]._key);
	goto back;
      }
      up_idx = HASH->Check(tmp);
      if ( up_idx == -1 ) {
	if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			 tmp, Prime->LineNumber());
	goto back;
      }

      tmp = Prime->Find_Value( Segment_Descriptor[1]._key); // "down:
      if ( !tmp) {
	if (out) fprintf(out, "Bummer: SEGMENT statement on line %d requires the \"%s\" field.\n",
			 Prime->LineNumber(), Segment_Descriptor[1]._key);
	goto back;
      }
      down_idx = HASH->Check(tmp);
      if ( down_idx == -1 ) {
	if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			 tmp, Prime->LineNumber());
	goto back;
      }
      tmp = Prime->Find_Value( Segment_Descriptor[2]._key); // "length"
      if ( !tmp ) {
	if (out) fprintf(out, "Bummer: SEGMENT statement on line %d requires the \"%s\" field.\n",
			 Prime->LineNumber(), Segment_Descriptor[2]._key);
	goto back;
      }
      if ( !is_scientific( tmp ) ) {
	if (out) fprintf(out, "Bummer: invalid real value specification on line %d\n",Prime->LineNumber());
	goto back;
      }
      length = atof(tmp);

      // ready to insert segment, use up_idx, down_idx, length, check length > 0!
      if ( length <= 0.0 ) {
	if (out) fprintf(out,"Bummer: found negative or zero length on line %d\n", Prime->LineNumber() );
	goto back;
      }
      if ( length <= OPT.LMin() ) {
	if (out) fprintf(out,"[WW]: segment length on line %d is %.4e, less than preferred minimal %.4e\n",
			 Prime->LineNumber(), length, OPT.LMin());
      }
      if ( length >= OPT.LMax() ) {
	if (out) fprintf(out,"[WW]: segment length on line %d is %.4e, larger than preferred maxium %.4e\n",
			 Prime->LineNumber(), length, OPT.LMax());
      }

      if (OPT.UseMetric() == 0 ) length *= OPT.FtoM();

      fp = SUB->GetNode( up_idx );
      tp = SUB->GetNode( down_idx );
      if (fp == NULL) {
	if (out) fprintf(out,"Bummer: internal error, unable to locate node at index %d for statement on line %d\n", 
			 up_idx, Prime->LineNumber());
	goto back;
      }      
      if (tp == NULL) {
	if (out) fprintf(out,"Bummer: internal error, unable to locate node at index %d for statement on line %d\n", 
			 down_idx, Prime->LineNumber());
	goto back;
      }

      // in case we have a lateral source
      src_id = lsources[up_idx];
      SUB->MakeStvEquation( up_idx, down_idx, length, src_id );

      break;

    case p_junction:
      good = Prime->Separate(out);
      if (good<0) {
	if (out) fprintf(out,"Bummer: syntax error in netlist specification at line %d\n", Prime->LineNumber());
	goto back;
      }
      
      tmp = Prime->Find_Value( Junc_Descriptor[0]._key); // "up1" 
      if (tmp) {
	up1_idx = HASH->Check(tmp);
	if ( up1_idx == -1 ) {
	  if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			   tmp, Prime->LineNumber());
	  goto back;
	}
      } else up1_idx = -1;


      tmp = Prime->Find_Value( Junc_Descriptor[1]._key); // "up2:
      if (tmp) {
	up2_idx = HASH->Check(tmp);
	if ( up2_idx == -1 ) {
	  if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			   tmp, Prime->LineNumber());
	  goto back;
	}
      } else up2_idx = -1;


      tmp = Prime->Find_Value( Junc_Descriptor[2]._key); // "down"
      if ( !tmp) {
	if (out) fprintf(out, "Bummer: JUNCTION statement on line %d requires the \"%s\" field.\n",
			 Prime->LineNumber(), Junc_Descriptor[2]._key);
	goto back;
      }
      down_idx = HASH->Check(tmp);
      if ( down_idx == -1 ) {
	if (out) fprintf(out, "Bummer: cannot find node \"%s\" specified on line %d\n",
			 tmp, Prime->LineNumber());
	goto back;
      }

      tmp = Prime->Find_Value( Junc_Descriptor[3]._key); // "coeff1"
      if (tmp) {
	if ( !is_scientific( tmp ) ) {
	  if (out) fprintf(out, "Bummer: invalid real value specification on line %d\n",Prime->LineNumber());
	  goto back;
	}
	coef1 = atof(tmp);
	if ( coef1 <= 0.0 ) {
	  if (out) fprintf(out, "Bummer: coefficient should be positive on line %d\n",Prime->LineNumber());
	  goto back;
	}
      } else coef1 = -1.0;

      tmp = Prime->Find_Value( Junc_Descriptor[4]._key); // "coeff2"
      if ( tmp) {
	if ( !is_scientific( tmp ) ) {
	  if (out) fprintf(out, "Bummer: invalid real value specification on line %d\n",Prime->LineNumber());
	  goto back;
	}
	coef2 = atof(tmp);
	if ( coef2 <= 0.0 ) {
	  if (out) fprintf(out, "Bummer: coefficient should be positive on line %d\n",Prime->LineNumber());
	  goto back;
	}
      } else coef2 = -1.0;

      // ready to insert junction, use up1_idx, up2_idx, down_idx, coef1, coef2, check
      if ( up1_idx == -1 && up2_idx == -1 ) {
	if (out) fprintf(out,"Bummer: at least one upstream nodes has to be specified on line %d\n",
			 Prime->LineNumber() );
	goto back;
      }

      if ( coef1 < 0 && coef2 < 0 ) {
	if (out) fprintf(out,"Bummer: at least one upstream coefficients has to be specified on line %d\n",
			 Prime->LineNumber() );
	goto back;
      }

      if ( (up1_idx > 0 && coef1 < 0) || ( up2_idx>0 && coef2<0 ) ) {
	if (out) fprintf(out,"Bummer: mismatch in upstream node and coefficient specifications on line %d\n",
			 Prime->LineNumber());
	goto back;
      }

      /* To add: 
	 distance check among three nodes, use X() and Y() of each node 
       */
      int nn[2];
      double rr[2];

      if ( up1_idx > 0 && up2_idx > 0 ) {  // both
	nn[0] = up1_idx; 
	nn[1] = up2_idx;
	rr[0] = coef1;
	rr[1] = coef2;
	SUB->MakeDepEquation( down_idx, 2, &nn[0], &rr[0]);

      } else if ( up1_idx > 0 ) {    // up1 only
	nn[0] = up1_idx;
	rr[0] = coef1;
	SUB->MakeDepEquation( down_idx, 1, &nn[0], &rr[0]);

      } else {    // up2 only
	nn[0] = up2_idx;
	rr[0] = coef2;
	SUB->MakeDepEquation( down_idx, 1, &nn[0], &rr[0]);

      }

      break;

    default:
      break;
    }
  }
  
  // insert Qsrc, lateral source and Bdn
  for (int jj=0; jj<STAT.N_Nodes(); jj++) {
    if ( qsources[jj] > -1) {
      if ( SUB->GetSource( qsources[jj])->GetParent() >= 0 ) {
	if (out) fprintf(out, "Bummer: trying to re-use qsource at node index %d\n",jj);
	goto back;
      }
      SUB->MakeQrSrcEquation( jj, qsources[jj] );
    }
  }

  
  for (int jj=0; jj<STAT.N_Nodes(); jj++) {
    if ( asources[jj] > -1 ) {
      if ( SUB->GetSource( asources[jj] )->GetParent() >= 0 ) {
	if (out) fprintf(out, "Bummer: trying to re-use asource at node index %d\n",jj);
	goto back;
      }
      SUB->MakeArSrcEquation( jj, asources[jj]);
      break;    // we assume there is only ONE BdnCondition
    }
  }

  r_code = 0; // only at this point we indicate everything is clean

  // calculate the chksum in case we are going to load something from the data
  for (int jj=0; jj<STAT.N_Nodes(); jj++) {
    mychk.add_string( NODE_NAMES->Get(jj) );
  }

  for (int jj=0; jj<STAT.N_Nodes(); jj++) {
    if ( qsources[jj] > -1 ) {
      Source *s = (Source*)SUB->GetSource( qsources[jj] );
      sprintf(chktmp, "%+.1e", s->Evaluate(0.0));
      mychk.add_string(chktmp);
    }
  }
  for (int jj=0; jj<STAT.N_Nodes(); jj++) {
    if ( asources[jj] > -1 ) {
      Source *s = (Source*)SUB->GetSource( asources[jj]);
      sprintf(chktmp, "%+.1e", s->Evaluate(0.0));
      mychk.add_string(chktmp);
    }
  }
  for (int jj=0; jj<STAT.N_Nodes(); jj++) {
    if ( lsources[jj] > -1 ) {
      Source *s = (Source*)SUB->GetSource( lsources[jj] );
      sprintf(chktmp, "%+.1e", s->Evaluate(0.0));
      mychk.add_string(chktmp);
    }
  }
  
  STAT.SetChkSum( mychk.HexCode() );

 back:
  if (Prime) delete Prime;
  if (Second) delete Second;

  return r_code;
}

// equivalent of shell basename() 
char *get_basename(char *path) {
  char *base = strrchr(path,'/');
  return base ? base+1 : path;
}
// end

