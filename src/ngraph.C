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
   methods for ngraph class
 */

#include <stdio.h>

#include "ngraph.h"

// constant
const double ngraph::_tol = 1e-4;

// take care of memory allocation
void ngraph::Init(int max_node_idx, int max_edge_idx, int max_junc_idx) {
  int jj;
  _nodes.size(max_node_idx);
  _edges.size(max_edge_idx);
  _junctions.size(max_junc_idx);
  for (jj = 0; jj < max_node_idx; jj++)
    _nodes[jj] = NULL;
  for (jj = 0; jj < max_edge_idx; jj++)
    _edges[jj] = NULL;
  for (jj = 0; jj < max_junc_idx; jj++)
    _junctions[jj] = NULL;
}

// insert a node, assign id
void ngraph::AddNode(int id) {
  vertex *vv = new vertex(id);
  _nodes[_num_nodes++] = vv;
}

// insert an edge, assign id
void ngraph::AddEdge(int up, int down, double len) {
  edge *ee;
  vertex *v1, *v2;

  assert(_nodes[up]);
  assert(_nodes[down]);
  assert(len > 0.0);

  v1 = _nodes[up];
  v2 = _nodes[down];

  ee = new edge(up, v1, v2, len);

  // set up the linked values in the vertices
  v1->SetDownEdge(ee);
  v2->SetUpEdge(ee);

  _edges[_num_edges++] = ee;
}

// insert a junction, one downstream node, many upstream nodes
void ngraph::AddJunction(int from, int num_inflow, int *to, double *r) {
  int jj;
  junction *junc = new junction(from);

  junc->SetFanin(_nodes[from]);
  for (jj = 0; jj < num_inflow; jj++) {
    junc->AddVertex(_nodes[to[jj]], r[jj]);
    _nodes[to[jj]]->SetOther(junc);
  }
  _nodes[from]->SetOther(junc);

  _junctions[_num_junctions++] = junc;
}

// insert a junction, many downstream nodes, one upstream node
void ngraph::AddJunction(int num_outflow, int *from, double *r, int to) {
  int jj;
  junction *junc = new junction(to);

  junc->SetFanin(_nodes[to]);
  for (jj = 0; jj < num_outflow; jj++) {
    junc->AddVertex(_nodes[from[jj]], r[jj]);
    _nodes[from[jj]]->SetOther(junc);
  }
  _nodes[to]->SetOther(junc);
  _junctions[_num_junctions++] = junc;
}

// compute downstream fan-out nodes
int ngraph::_FindDownStreamNodes(vertex *st, vertex **des, double *vals) {
  int rc = 0;

  if (st->GetDown() != NULL) {
    des[0] = st->GetDown()->GetTo();
    vals[0] = 1.0;
    rc = 1;
  } else {
    if (st->GetOther() == NULL) {
      rc = 0;
    } else {
      if (st->GetOther()->IsFanin(st)) {
        // branching point, we copy all fan-out vertices
        for (int jj = 0; jj < st->GetOther()->GetNum(); jj++) {
          des[jj] = st->GetOther()->GetVertex(jj);
          vals[jj] = st->GetOther()->GetRatio(jj);
        }
        rc = st->GetOther()->GetNum();
      } else {
        // converging, we go down the fan-in vertex
        des[0] = st->GetOther()->GetFanin();
        vals[0] = 1.0;
        rc = 1;
      }
    }
  }
  return rc;
}

// compute upstream fan-out nodes
int ngraph::_FindUpStreamNodes(vertex *st, vertex **des, double *vals) {
  int rc = 0;

  if (st->GetUp() != NULL) {
    des[0] = st->GetUp()->GetFrom();
    vals[0] = 1.0;
    rc = 1;
  } else {
    if (st->GetOther() == NULL) {
      rc = 0;
    } else {
      if (st->GetOther()->IsFanin(st)) {
        // branching point, we copy all fan-out vertices
        for (int jj = 0; jj < st->GetOther()->GetNum(); jj++) {
          des[jj] = st->GetOther()->GetVertex(jj);
          vals[jj] = st->GetOther()->GetRatio(jj);
        }
        rc = st->GetOther()->GetNum();
      } else {
        // converging, we go down the fan-in vertex
        des[0] = st->GetOther()->GetFanin();
        vals[0] = 1.0;
        rc = 1;
      }
    }
  }
  return rc;
}

/*
  T-Check
  Do the following:
    1. Identify the ROOT nodes ( with only one upstream edges )
    2. Identify the LEAF nodes ( with only one downstream edges )
    3. Propagate from all LEAF nodes down to all ROOT nodes;
    4. Propagate from all ROOT nodes up to all LEAF nodes;
    5. All nodes should be visited exactly twice. If not, raise the flag
  returncode:
      0 good
     <0 bad
*/
int ngraph::TChk() {
  int rc;
  int jj;

  edge *eup, *edn;

  rc = 0;
  for (jj = 0; jj < _num_nodes; jj++) {
    eup = _nodes[jj]->GetUp();
    edn = _nodes[jj]->GetDown();
    if (eup != NULL) {
      if (edn != NULL)
        _nodes[jj]->SetType(NORMAL);
      else
        _nodes[jj]->SetType(ROOT);
    } else {
      if (edn != NULL)
        _nodes[jj]->SetType(LEAF);
      else {
        _nodes[jj]->SetType(UNDEFINED);
        _nodes[jj]->SetFlag(HANGING);
        rc = -1;
      }
    }
  }
  if (rc != 0)
    return (rc); // hanging nodes, no need to continue

  // do a quick check to make sure no vertices with both edge and other
  // which should not happen anyway.
  // also set the vertex type to those vertices, they could be mislabelled as
  // leaf or root
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->GetOther() != NULL && _nodes[jj]->GetDown() != NULL &&
        _nodes[jj]->GetUp() != NULL) {
      _nodes[jj]->SetType(UNDEFINED);
      _nodes[jj]->SetFlag(WRONG_J);
      rc = -2; // error flag
    }
    if (_nodes[jj]->GetOther() != NULL)
      _nodes[jj]->SetType(IN_JUNC);
  }
  if (rc != 0)
    return (rc);

  // Each Leaf should have a positive injection
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->GetType() != LEAF)
      continue;
    if (_nodes[jj]->Injection() < this->_tol) {
      _nodes[jj]->SetFlag(NEED_Q);
      rc = -3;
    }
  }
  if (rc != 0)
    return (rc);

  // Each Root should have a positive A
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->GetType() != ROOT)
      continue;
    if (_nodes[jj]->SumUp() < this->_tol) {
      _nodes[jj]->SetFlag(NEED_A);
      rc = -4;
    }
  }
  if (rc != 0)
    return (rc);

  /* propagate Q downstream, from leaves to roots */
  // note injection can present at any nodes, not necessarily only on the leaf
  // starting from the leaves

  // two stacks, always work in sync
  MyStack<vertex *> NS; // node stack
  MyStack<double> VS;   // value stack

  vertex *nn;
  double vv = -1.0;
  TinyVec<vertex *, 64> tmp_nodes; // somehow unsafe, but 64 should be enough
  TinyVec<double, 64> tmp_vals;
  int junk;
  // start from leaves
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->GetType() == LEAF) {
      NS.Push(_nodes[jj]);
      VS.Push(_nodes[jj]->SumDown());
    }
  }

  while (1) {
    if (NS.Pop(nn) == -1)
      break;
    VS.Pop(vv);

    junk = this->_FindDownStreamNodes(nn, &tmp_nodes[0], &tmp_vals[0]);
    vv += nn->Injection();
    // we add the injection at this node when computing the Q, but it has to
    // propagate downstream. This is due to the definition of the Qlateral.
    // Except when the injection is at the LEAF, in which case it's a Qsource
    // and we have to keep it!
    if (LEAF == nn->GetType()) {
      nn->SumDown() += vv;
    } else {
      nn->SumDown() += vv - nn->Injection();
    }
    for (int jj = 0; jj < junk; jj++) {
      NS.Push(tmp_nodes[jj]);
      VS.Push(vv * tmp_vals[jj]);
    }
  }

  // ditto for upstream traversing, but start from roots
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->GetType() == ROOT) {
      NS.Push(_nodes[jj]);
      VS.Push(_nodes[jj]->SumUp());
    }
  }

  while (1) {
    if (NS.Pop(nn) == -1)
      break;
    VS.Pop(vv);

    junk = this->_FindUpStreamNodes(nn, &tmp_nodes[0], &tmp_vals[0]);
    nn->SumUp() = vv; // upstream is slightly different, no accumulation

    for (int jj = 0; jj < junk; jj++) {
      NS.Push(tmp_nodes[jj]);
      VS.Push(vv * tmp_vals[jj]);
    }
  }

  // both sumup and sumdown should be nonzero at all nodes
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->SumDown() < this->_tol) {
      _nodes[jj]->SetFlag(NO_QPATH);
      rc = -5;
    }
    if (_nodes[jj]->SumUp() < _tol) {
      _nodes[jj]->SetFlag(NO_APATH);
      rc = -5;
    }
#if 0
    printf("%d: %.5f %.5f\n", jj, _nodes[jj]->SumUp(), _nodes[jj]->SumDown());
#endif
  }
  if (rc != 0)
    return (rc);

  _chk_ok = 1;
  return (rc);
}

// print an error message to the file handle
void ngraph::PrintErrorMsg(FILE *F, const char *header, char **trtable) {
  int jj;
  for (jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj] == NULL)
      continue; // just in case the id's are not continuous

    if (trtable) {
      if (_nodes[jj]->GetFlag() == HANGING) {
        fprintf(F, "%s node \"%s\" is floating.\n", header,
                trtable[_nodes[jj]->GetId()]);
      }

      if (_nodes[jj]->GetFlag() == WRONG_J) {
        fprintf(F, "%s node \"%s\" has wrong junction.\n", header,
                trtable[_nodes[jj]->GetId()]);
      }

      if (_nodes[jj]->GetFlag() == NEED_A) {
        fprintf(F, "%s node \"%s\" requires outpouring boundary condition.\n",
                header, trtable[_nodes[jj]->GetId()]);
      }

      if (_nodes[jj]->GetFlag() == NEED_Q) {
        fprintf(F, "%s node \"%s\" requires upstreaming Q.\n", header,
                trtable[_nodes[jj]->GetId()]);
      }

      if (_nodes[jj]->GetFlag() == NO_QPATH) {
        fprintf(F, "%s node \"%s\" no flow path to upstream.\n", header,
                trtable[_nodes[jj]->GetId()]);
      }

      if (_nodes[jj]->GetFlag() == NO_APATH) {
        fprintf(F, "%s node \"%s\" no flow path to downsteram.\n", header,
                trtable[_nodes[jj]->GetId()]);
      }

    } else {

      if (_nodes[jj]->GetFlag() == HANGING) {
        fprintf(F, "%s node \"%d\" is floating.\n", header,
                _nodes[jj]->GetId());
      }

      if (_nodes[jj]->GetFlag() == WRONG_J) {
        fprintf(F, "%s node \"%d\" has wrong junction.\n", header,
                _nodes[jj]->GetId());
      }

      if (_nodes[jj]->GetFlag() == NEED_A) {
        fprintf(F, "%s node \"%d\" requires outpouring boundary condition.\n",
                header, _nodes[jj]->GetId());
      }

      if (_nodes[jj]->GetFlag() == NEED_Q) {
        fprintf(F, "%s node \"%d\" requires upstreaming Q.\n", header,
                _nodes[jj]->GetId());
      }

      if (_nodes[jj]->GetFlag() == NO_QPATH) {
        fprintf(F, "%s node \"%d\" no flow path to upstream.\n", header,
                _nodes[jj]->GetId());
      }

      if (_nodes[jj]->GetFlag() == NO_APATH) {
        fprintf(F, "%s node \"%d\" no flow path to downsteram.\n", header,
                _nodes[jj]->GetId());
      }
    }
  }
}

// returns the number of roots, as well as a list of root indices
// called owns memory
int ngraph::GetRoots(int *arr) {
  if (!_chk_ok) {
    arr[0] = 0;
    return (-1.0);
  }
  int rc = 0;
  for (int jj = 0; jj < _num_nodes; jj++) {
    if (_nodes[jj]->GetType() == ROOT)
      arr[rc++] = _nodes[jj]->GetId();
  }
  return (rc);
}

// return upstream length
int ngraph::GetUpstreamLength(int j_idx, double &len) {
  assert(j_idx < _num_nodes);

  len = 0.0;
  vertex *n = _nodes[j_idx];

  if (n->GetUp() != NULL) {
    len = n->GetUp()->Len();
    return (0);
  } else if (n->GetOther() != NULL) {
    return (n->GetOther()->GetNum());
  } else {
    return (-1);
  }
}

// just for junctions, returns the indices of upstream fan-outs
// also returns the number of fanouts
// caller owns memory
int ngraph::GetJunctionFanouts(int j_idx, int *fo, double *r) {
  assert(j_idx < _num_nodes);

  vertex *n = _nodes[j_idx];

  if (n->GetUp() != NULL || n->GetOther() == NULL) {
    fo[0] = -1;
    r[0] = 0.0;
    return (0);
  }
  junction *junc = n->GetOther();

  int num = junc->GetNum();
  for (int jj = 0; jj < num; jj++) {
    fo[jj] = junc->GetVertex(jj)->GetId();
    r[jj] = junc->GetRatio(jj);
  }
  return (num);
}

// debug functions
void ngraph::PrintAll(FILE *F) {
  if (!F)
    return;

  fprintf(F, "**** T-check is %s\n", _chk_ok ? "done" : "not done");

  for (int jj = 0; jj < _num_nodes; jj++) {
    fprintf(F, "%3d: [%4d]: %1d %1d %.2f %.2f %.2f\n", jj, _nodes[jj]->GetId(),
            _nodes[jj]->GetType(), _nodes[jj]->GetFlag(),
            _nodes[jj]->Injection(), _nodes[jj]->SumDown(),
            _nodes[jj]->SumUp());
  }
}
// Local Variables:
// mode: c++
// End:
