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

/***
 a simple graph package
 
 functionality:
    t-check;
    propagates Q downstream and A upstream;
    returns the length to the upsteam node;

 **********************************************/

#ifndef _N_GRAPH_H
#define _N_GRAPH_H

#include <assert.h>
#include <string.h>

#define DEFAULT_STACKSIZE 255
#define DEFAULT_ARRAYSIZE 512
#define DEFAULT_VERTEXSIZE 16
#define DEFAULT_VERTEXPERJUNC 4

enum V_TYPE   {NORMAL=0, LEAF, ROOT, IN_JUNC, UNDEFINED};
enum ERR_TYPE {TCHK_OK=0, WRONG_J, HANGING, NEED_A, NEED_Q, NO_QPATH, NO_APATH };

/* TinyVec: a small flexible array to store vertex array */
template<class T, int FIXSZ=DEFAULT_VERTEXSIZE> 
class TinyVec {
protected:
  T   _fixmem[FIXSZ];
  int _asz;
  T*  _vec;
  
  // if not large enough, allocates more than requested and copies the old data into new memory
  int _do_grow( int sz ) { 
    T* newvec = new T[sz + _asz]; // grow the memory block size
    for( int k=0; k!=_asz; ++k) newvec[k] = _vec[k];
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = newvec;
    _asz += sz;                            
    return _asz;
  }
  
public:
  /* constructors and destructor */
  TinyVec() : _asz(FIXSZ), _vec(_fixmem) {};
  TinyVec(int sz) { 
    if (sz <= FIXSZ) {
      _vec = _fixmem;
      _asz = FIXSZ;
    } else {
      _vec = new T[sz];
      _asz = sz;
    }
  };
  
  TinyVec( const TinyVec<T,FIXSZ>& ) { assert(0); } // disallow copy constructor
  ~TinyVec() {
    if (_vec != (T*)_fixmem) delete[] _vec; 
  };
  
  /* access methods */
  operator T*() { return _vec; }
  T& operator[](int i) { assert(i<_asz); assert(i>=0); return _vec[i]; }; 
  
  int size(int sz) { // if not large enough allocates exactly as requested
    if ( sz <= _asz ) return _asz; 
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = new T[sz];                
    _asz = sz;                            
    return _asz;
  }
  
  int resize( int sz ) { // if not large enough allocates more than requested
    if ( sz <= _asz ) return _asz; 
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = new T[sz + _asz]; // at least doubles the allocated size                
    _asz += sz;                            
    return _asz;
  }
  
  inline int grow (int sz) { return (sz<_asz)?_asz:_do_grow(sz); }
};

/* A simple stack, for DFS search
   two basic methods: push() and pop()
   additionally: 1) it's possible to check the current depth of the stack
   2) it's possible to nuke the content within the stack
*/
template<class T> 
class MyStack {
  
protected:
  T                               *_array;
  int                              _size;
  int                              _current;
  
public:
  MyStack():_size(DEFAULT_STACKSIZE),_current(0) {
    _array = new T [DEFAULT_STACKSIZE];
  }
  ~MyStack() {
    if (_array) delete [] _array;
  }
  
  void Push(T data) {
    int  newsz;
    T   *tmp;
    if ( _current >= _size ) {  // if insufficient, grow and copy
      newsz = 2*_size;
      tmp = new T[newsz];
      memcpy(tmp, _array, _size*sizeof(T));
      delete[] _array;
      _array = tmp;
      _size = newsz;
    }
    _array[_current++] = data;
  }
  
  /* return top of the stack by reference
     return code:
     rc ==  1: good
     rc == -1: stack is empty 
  */
  int Pop(T& data) {  // rc > 1, good
    if ( _current > 0 ) {
      data = _array[--_current];
      return 1;
    } else {
      return -1;
    }
  }
  
  int StackDepth(void) const {
    return _current;
  }
  
  void Nuke(void) {  // flush what ever in the stack 
    _current = 0;
  }
};

/*
  Each edge has two vertices, one upstream, one downstream;
  Each vertex has one or two edges, which defines a LEAF or ROOT;
  Use Junction to link the edges together
*/
class edge;
class junction;

class vertex {
protected:
  int              _id;
  V_TYPE           _type;
  edge*            _up;
  edge*            _down;
  junction*        _other;
  // flag
  ERR_TYPE         _flag;
  
  // we can propagate these quantaties
  double           _injection;
  double           _sum_down;
  double           _sum_up;

public:
  vertex(): _id(-1),_type(UNDEFINED), _up(NULL), _down(NULL),_other(NULL),
	    _flag(TCHK_OK),_injection(0.0), _sum_down(0.0), _sum_up(0.0) {}
  vertex(int id):_type(UNDEFINED),_up(NULL),_down(NULL),_other(NULL),_flag(TCHK_OK),_injection(0.0),_sum_down(0.0),_sum_up(0.0) {
    _id = id;
  }
  vertex(int id, edge* U, edge *D) {
    _id = id;
    _up = U;
    _down = D;
    _other = NULL;
    _flag = TCHK_OK;
    _type = UNDEFINED;
    _injection = 0.0;
    _sum_down = 0.0;
    _sum_up = 0.0;
  }
  vertex( const vertex& ) { assert(0); } // cannot copy 
  ~vertex() {}
  
  // methods
  int GetId() const  { return _id; } 
  void SetId(int id) { _id = id;  }
  ERR_TYPE GetFlag() const { return _flag; }
  void SetFlag(ERR_TYPE f) { _flag = f; }
  
  void SetDownEdge(edge* e) { _down = e; }
  void SetUpEdge(edge* e)   { _up = e; }
  edge* GetDown()    { return _down; }
  edge* GetUp()      { return _up;   }
  
  junction* GetOther() { return _other; }
  void SetOther(junction* junc) { _other = junc; }
  
  V_TYPE  GetType() const { return _type; }
  void SetType(V_TYPE v) { _type = v; }
  double& Injection() { return _injection; }
  double& SumDown()   { return _sum_down; }
  double& SumUp()     { return _sum_up;   }
};

/*
  Each edge *must* have two vertices 
*/
class edge {
protected:
  int              _id;
  unsigned int     _strahler;  // we can compute strahler number, but it's not done yet
  vertex*          _from;
  vertex*          _to;
  
  double           _len; 
  
public:
  edge(): _id(-1), _strahler(0),_from(NULL), _to(NULL), _len(0.0) {}
  edge(int id): _strahler(0), _from(NULL), _to(NULL),_len(0.0) {
    _id = id;
  }
  edge(int id, vertex* from, vertex* to, double len) {
    _id = id;
    _strahler = 0;
    _from = from;
    _to = to;
    _len = len;
  }
  
  edge( const edge& ) { assert(0); } // cannot copy
  ~edge() {}
  
  // methods
  int GetId() const { return _id; }
  void SetId(int id) { _id = id;  }
  vertex* GetFrom()  { return _from; }
  vertex* GetTo()    { return _to;   }
  unsigned int& Strahler() { return _strahler; }
  double&  Len()    { return _len; }
};

/* 
   A junction connects multiple edges together,
   from upstream to downstream, it can be either N-to-1, or 1-to-N
*/

class junction {
protected:
  int        _id;
  int        _num;
  vertex*    _fanin;
  TinyVec< vertex*,DEFAULT_VERTEXPERJUNC>  _fanouts;
  TinyVec< double, DEFAULT_VERTEXPERJUNC>  _ratios;
  
  
public:
  junction():_id(0),_num(0), _fanin(NULL) {}
  junction(int id) { _id=id; _num=0; _fanin=NULL; }
  junction( const junction&) { assert(0); } // cannot copy
  ~junction() {}
  
  // methods
  // specify the fan-in vertex
  void SetFanin(vertex *v) { _fanin = v; }
  int IsFanin(vertex *v) const { return (v==_fanin); }
  vertex* GetFanin() { return (_fanin); }
  
  // add fan-out vertices
  void AddVertex(vertex* v, double r) {
    _fanouts.resize(_num+1);
    _ratios.resize(_num+1);  // just in case there are more than default
    _fanouts[_num] = v;
    _ratios[_num++] = r;
  }
  
  int GetNum() const { return _num; }

  // get the fan-out vertices one-by-one
  vertex* GetVertex(int jj) {
    assert( jj < _num );
    return _fanouts[jj];
  }
  
  double GetRatio(int jj) {
    assert ( jj < _num );
    return _ratios[jj];
  }
  
};

/* 
   ngraph stores the whole graph, also provides methods to do t-check
   as well as method to facilitate the traversal

   The class should work in junction with a string->integer hashtable
   Each node should have a unique id, which corresponds to the order they are inserted
     starting from 0

*/
class ngraph {
  
private:
  int                                   _num_nodes;
  int                                   _num_edges;
  int                                   _num_junctions;
  
  int                                   _chk_ok;
  
  TinyVec<vertex*,  DEFAULT_ARRAYSIZE>  _nodes;
  TinyVec<edge*,    DEFAULT_ARRAYSIZE>  _edges;
  TinyVec<junction*,DEFAULT_ARRAYSIZE>  _junctions;
  
  static const double                   _tol;

  // start from a node, find its all downstream fan-outs
  // return the number of fan-outs
  //    for normal nodes, there should only one downstream fan-out
  //    except for the branchking "junctions"
  int _FindDownStreamNodes(vertex* st, vertex** des, double *vals);
  
  // ditto, this time goes upstream
  int _FindUpStreamNodes(vertex* st, vertex** des, double *vals);
  
public:
  ngraph():_num_nodes(0),_num_edges(0),_num_junctions(0),_chk_ok(0) {}
  ~ngraph() {
    int jj;
    for (jj=0; jj<_num_nodes;jj++) { if (_nodes[jj]) delete _nodes[jj]; } 
    for (jj=0; jj<_num_edges;jj++) { if (_edges[jj]) delete _edges[jj]; } 
    for (jj=0; jj<_num_junctions;jj++) { if (_junctions[jj]) delete _junctions[jj]; } 
  }
  
  // methods
  void Init(int max_node_idx, int max_edge_idx, int max_junc_idx);
  
  // name says all
  void AddNode(int id);
  void AddEdge(int up, int down, double len);
  
  // junction 1-to-many 
  void AddJunction(int from, int num_outflow, int* to, double *r);
  // junction many-to-1
  void AddJunction(int num_inflow, int* from, double *r, int to);
  
  // add quantities for downstream propagation
  void InsertQ(int id, double q) {
    assert ( _nodes[id] );
    _nodes[id]->Injection() = q;
  }
  
  // add quantities for upstream propagation
  void InsertA(int id, double a) {
    assert ( _nodes[id] );
    _nodes[id]->SumUp() = a;
  }
  
  
  // toplevel method for t-checking
  //   rc:  0 if everything is fine
  //       <0 otherwise, check Flag() of each vertex for detailed explanations
  int TChk();
  
  // call this function after Tchk(), 
  // otherwise it's useless
  int GetNumRoots() {
    if (!_chk_ok) return -1;
    int rc = 0;
    for (int jj=0; jj<_num_nodes; jj++) rc += (_nodes[jj]->GetType()==ROOT);
    return rc;
  }
  
  // Get the roots, also returns the number of roots
  // called owns memory
  int GetRoots(int* arr);
  
  // Get propagated Q, useless if TChk() hasn't been run
  // returns negative value if no good
  double GetQ(int j_idx) {
    if (!_chk_ok) return -1.0;  // anything negative is bad
    assert ( j_idx < _num_nodes );
    return (_nodes[j_idx]->SumDown());
  }

  // Retrieve propagated Q, 
  // for compatibility reason only
  double GetA(int j_idx) {
    if (!_chk_ok) return -1.0;  // anything negative is bad
    assert ( j_idx < _num_nodes );
    return (_nodes[j_idx]->SumUp());
  }

  // Get the length of upstream edge
  // return: 0 - everything OK
  //        -1 - end
  //        >0 - actually a junction point, returns number of upstream vertices
  int GetUpstreamLength(int j_idx, double &len);

  // this only applies to a non-leaf, non-junction node
  // otherwise we will get -1
  // should use GetUpstreamLength() first to get type
  int GetUpstreamNode(int j_idx) {
    assert ( j_idx < _num_nodes );
    vertex *n=_nodes[j_idx];
    if ( n->GetUp() == NULL ) return (-1); 
    return ( n->GetUp()->GetId() );
  }

  // In case of junction, get indices and ratios of the upstream points
  int GetJunctionFanouts(int j_idx, int* fanouts, double *ratios);

  // print detailed topological error message, 
  // pass the translation table for the node names in trtable,
  // if set trtable to NULL if only the node indices are needed
  void PrintErrorMsg(FILE *F, const char* header, char** trtable);

  // these methods are for debug only
  void PrintAll(FILE *F);
  int GetNumNodes() const { return _num_nodes; }
  int GetNumEdges() const { return _num_edges; }
  int GetNumJunctions() const { return _num_junctions; }

};

#endif

// Local Variables:
// mode: c++
// End:

