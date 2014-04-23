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

// simple wrapper around the sys/time to get execution time

#ifndef _SPTTIMEMEAS_H
#define _SPTTIMEMEAS_H

#include <sys/time.h>

/*
 simple wrapper class to measure runtime
 usage example:
    start(); 
    stop(); 
    read(); 
 repeat
 or 
  restart(); 
  stop(); 
  start(); 
  stop(); 
  read(); repeat

 limitation: the resoltion only to milliseconds. 

 make sure to check the return flag of stop() 
*/
class mytm {
protected:
  double           _accum;
  int              _started;
  struct timeval   _start;
  struct timeval   _end;
  
public:
  mytm():_accum(0.0),_started(0) {}
  ~mytm() {}

  void start() {
    _started = 1;
    gettimeofday(&_start,NULL);
  }
  
  int stop() {
    if (!_started) {
      return 0;
    } else {
      gettimeofday(&_end, NULL);
      long int dus, ds;
      dus = _end.tv_usec - _start.tv_usec;
      ds = _end.tv_sec - _start.tv_sec;
      _accum += (double)(ds*1000000+dus)/1000000.0;
      _started = 0;
      return 1;
    }
  }
  
  void restart() {
    _accum =0.0;
    start();
  }
  
  double read() const { return _accum; }

};

#endif

// Local Variables:
// mode: c++
// End:

