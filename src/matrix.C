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

/************************************************************************
  History
  05/05/11    FYL       Created
*************************************************************************/

#include    <assert.h>

#include    "matrix.h"

/* common facility methods */
void matrix_dump(FILE *fp, char* mnm, int nr, int nc, double* h, int ldh) {
  fprintf( fp, "\n%s = [", mnm );
  for (int j=0; j<nr; ++j) {
    fprintf( fp, "\t" );
    for (int k=0; k<nc; ++k) fprintf( fp, "%le ", h[ldh*k+j] );
    fprintf( fp, "\n" );
  }
  fprintf( fp, "];\n" );
}

void matrix_dump(FILE *fp, char* mnm, int nr, int nc, dComplex* h, int ldh ) {
  fprintf( fp, "\n%s = [", mnm );
  for (int j=0; j<nr; ++j) {
    fprintf( fp, "\t" );
    for (int k=0; k<nc; ++k) {
      double r = real(h[ldh*k+j]);
      double i = imag(h[ldh*k+j]);
      fprintf( fp, "%le + %lei ", r, i );
    }
    fprintf( fp, "\n" );
  }
  fprintf( fp, "];\n" );
}

void vector_dump(FILE *fp, char* mnm, int nr, double* v) {
  fprintf( fp, "\n%s = [", mnm );
  for (int k=0; k<nr; ++k) 
    fprintf( fp, "%le ", v[k] );
  fprintf( fp, "];\n" );
}

void vector_dump(FILE *fp, char* mnm, int nr, dComplex* v) {
  fprintf( fp, "\n%s = [", mnm );
  for (int k=0; k<nr; ++k) {
    double r = real(v[k]);
    double i = imag(v[k]);
    fprintf( fp, "%le + %lei ", r, i );
  }
  fprintf( fp, "];\n" );
}

void DiagMatrix::SimplePrint(FILE *of) {
  assert( _dim[0] >= 0 );
  fprintf(of, "======= row: %d col:%d =============\n",Ndim(0), Ndim(1) );
  fprintf(of, "=========== diagonal entries ===============\n");
  for (int i=0; i!=Nnz(); i++ ) {
    fprintf(of, "%-4d %-4d %lg\n", i, i, Data(i));
  }
}

void DiagMatrix::MatlabDump(FILE *of, const char* name) {
  assert( _dim[0] >= 0 );
  fprintf(of, "%sdiag = [", name); 
  for( int k=0; k!=Nnz(); ++k ) fprintf(of, "%le\n", Data(k));
  fprintf(of, "];\n%s = sparse(1:%d,1:%d,%sdiag,%d,%d);\n\n", name, Nnz(), Nnz(), name, Ndim(0), Ndim(1));
  fprintf(of, "clear %sdiag;\n", name );
}


/* construct identity matrix */
void SparseMatrix::Eye(int rows, int cols, double alpha) {
  _dim[0] = rows;
  _dim[1] = cols;
  _size   = MIN(rows,cols);
  MakeSpace(_size);
  for (int i=0; i<_size; i++) {
    Row(i) = i; 
    Col(i) = i;
    Data(i) = alpha;
  }
  _crow.size(Ndim(0)+1);
  for( int kr=0; kr<=_size; kr++ ) _crow[kr] = kr;
  for( int kr=cols+1; kr<=rows; kr++ ) _crow[kr] = _size;
  _is_sorted = 1;
}

/* 
   printing functions
 */
void SparseMatrix::SimplePrint(FILE *of) {  
  assert( _dim[0] >= 0 );
  fprintf(of, "======= row: %d col:%d nnz: %d =============\n",Ndim(0), Ndim(1), Nnz());
  fprintf(of, "=========== non-zero entries ===============\n");
  for (int i=0; i<Nnz(); i++ ) {
    fprintf(of, "%-4d %-4d %lg\n", Row(i), Col(i), Data(i));
  }
}

void SparseMatrix::MatlabDump(FILE *of,const char* name) {
  assert( _dim[0] >= 0 );
  fprintf(of, "%srow = [", name); 
  for( int k=0; k!=Nnz(); ++k ) fprintf(of, "%d\n", Row(k)+1);
  fprintf(of, "];\n%scol = [", name); 
  for( int k=0; k!=Nnz(); ++k ) fprintf(of, "%d\n", Col(k)+1);
  fprintf(of, "];\n%smat = [", name); 
  for( int k=0; k!=Nnz(); ++k ) fprintf(of, "%le\n", Data(k));
  fprintf(of, "];\n%s = sparse(%srow,%scol,%smat,%d,%d);\n\n", name, name, name, name, Ndim(0), Ndim(1));
  fprintf(of, "clear %srow;\nclear %scol;\nclear %smat;\n", name, name, name );
}


void SelMatrix::MatlabDump( FILE* fn, const char* nm )
{
  fprintf( fn, "%sRows = [\n", nm ); 
  for( int k=0; k!=_sz; ++k ) if (_sm[k]!=-1) fprintf( fn, "%d\n", k );
  fprintf( fn, "];\n%sCols = [\n", nm ); 
  for( int k=0; k!=_sz; ++k ) if (_sm[k]!=-1) fprintf( fn, "%d\n", _sm[k] );
  fprintf( fn, "];\n%s=sparse( %sRows+1, %sCols+1, ones(size(%sRows)), %d, %d);\n"
	   , nm, nm, nm, nm, _dim[0], _dim[1] );
  fprintf( fn, "clear %sRows %sCols;\n", nm, nm );
}

void AdjMatrix::SimplePrint( FILE* fn ) 
{
  fprintf( fn,   "Adjacency matrx %d %d\n", Ndim(0), Ndim(1) );
  for( int k=0; k!=_nr; ++k )
    fprintf( fn, "%d %d %d\n", k, Np(k), Nn(k) );
}

void AdjMatrix::MatlabDump( FILE* fn, const char* nm )
{
  fprintf( fn,   "%sRows = [\n", nm ); 
  for( int k=0; k!=_nr; ++k ) {
    if (Np(k)!=-1) fprintf( fn, "%d\n", k );
    if (Nn(k)!=-1) fprintf( fn, "%d\n", k );
  }
  fprintf( fn, "];\n%sCols = [\n", nm ); 
  for( int k=0; k!=_nr; ++k ) {
    if (Np(k)!=-1) fprintf( fn, "%d\n", Np(k) );
    if (Nn(k)!=-1) fprintf( fn, "%d\n", Nn(k) );
  }
  fprintf( fn, "];\n%sData = [\n", nm ); 
  for( int k=0; k!=_nr; ++k ) {
    if (Np(k)!=-1) fprintf( fn, "%d\n", 1  );
    if (Nn(k)!=-1) fprintf( fn, "%d\n", -1 );
  }
  fprintf( fn, "];\n%s=sparse( %sRows+1, %sCols+1, %sData, %d, %d);\n", 
	   nm, nm, nm, nm, Ndim(0), Ndim(1) );
  fprintf( fn, "clear %sRows %sCols %sData;\n", nm, nm, nm );
}

// End

