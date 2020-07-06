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

/* class to do check sum */

#ifndef _SPT_CHKSUM_H
#define _SPT_CHKSUM_H

#include <stdlib.h>
#include <time.h>

/****************************************
 a simple class to check multiple char strings
   Typical usage:
      CheckSum cs;
      cs.add_string(char*);
      cs.add_string(char*);
      print hexcode() or
      print checksum()
      cs.reset();
      repeat
   can be copied or assigned, in that case the internal checksum will be
inherited change the sequence of adding checksum will not change the check sum
   hence
      cs.add_string(str1); cs.add_string(str2)

   generates identical results as in
      cs.add_string(str2); cs.add_string(str1)

   to compare two different files, the user should also provide the "seed"
   cs.set_seed( seed );

   read large array of strings and calculate the checksum
     while ( fgets(buf, 1023, F) != NULL ) {
      cs.add_string( buf );
     }

     printf("checksum = %x, seed = %x, hex = %s\n", cs.CheckSum(), cs.Seed(),
cs.HexCode() );


   unsigned int overflow could happen.

   rough runtime: 500,000 records, 40 ms, including disk I/O
************************************************/

const char mycheck_hex[] = "0123456789abcdef";
const int mycheck_init = 0x12345678;

class ChkSum {
protected:
  int _counter;
  unsigned int _chk;
  unsigned int _seed;
  char _buf[9];
  void _gen_seed() {
    srand((unsigned)time(0));
    _seed = (unsigned int)rand();
  }

public:
  ChkSum() : _counter(0) {
    _seed = mycheck_init;
    _chk = _seed;
  }
  ChkSum(const ChkSum &other)
      : _counter(0), _chk(other._chk), _seed(other._seed) {}
  ChkSum operator=(const ChkSum &rhs) {
    _counter = 0;
    _chk = rhs._chk;
    _seed = rhs._seed;
    return *this;
  }
  ~ChkSum() {}

  int Counter() const { return _counter; }
  unsigned int Seed() const { return _seed; }
  unsigned int CheckSum() const { return _chk; }
  char *HexCode() {
    int jj;
    unsigned int tmp;
    tmp = _chk;
    for (jj = 0; jj < 4; jj++) {
      _buf[(jj << 1)] = mycheck_hex[((tmp % 256) >> 4)];
      _buf[(jj << 1) + 1] = mycheck_hex[((tmp % 256) % 16)];
      tmp = (tmp >> 8);
    }
    _buf[8] = 0; // hard coded size;
    return _buf;
  }

  void reset() {
    _counter = 0;
    _seed = mycheck_init;
    _chk = _seed;
  }
  void reset(unsigned int s) {
    _counter = 0;
    _seed = s;
    _chk = _seed;
  }

  void set_seed(int s) { _seed = s; }
  void rnd_seed(void) {
    _gen_seed();
    _chk = _seed;
  }

  void add_string(char *str) {
    for (int i = 0; str[i] != '\0'; i++)
      _chk += ((int)str[i] * (i + 1));
    _counter++;
  }
};

#endif

// Local Variables:
// mode: c++
// End:
