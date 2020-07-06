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

/* implementation of x-section translation */

#include "xstrans.h"

// needed by qsort
int mycompare(const void *a, const void *b) {
  return (((segment *)a)->_y1 - ((segment *)b)->_y1);
}

int li_compare(const void *a, const void *b) {
  return (*(long int *)a - *(long int *)b);
}

void xs_trans::_init_internal() {
  // assuming min and max are available
  long int sz;
  int nsz;

  sz = (_max_y - _min_y) / YGRID;
  nsz = (int)((double)sz * 1.10);

  _my_aa.size(nsz);
  _my_ll.size(nsz);
  _my_pp.size(nsz);
  _my_yy.size(nsz);
  _my_cc.size(nsz);

  _aa.size(nsz);
  _ll.size(nsz);
  _pp.size(nsz);
  _cc.size(nsz);
  _yy.size(nsz);
}

int xs_trans::_record_into_int(int num, double *xin, double *yin, double s) {
  int jj;

  for (jj = 0; jj < num; jj++) {
    _xxin[jj] = (long int)(xin[jj] / s);
    _yyin[jj] = (long int)(yin[jj] / s);
  }
  return jj;
}

int xs_trans::_cmethod(double t) {
  int cntr, jj;
  long int xs, ys, xp, yp;
  double angp, angs;

  _rxx[0] = _xxin[0];
  _ryy[0] = _yyin[0];

  xs = _xxin[0];
  ys = _yyin[0];
  cntr = 1;

  xp = _xxin[1];
  yp = _yyin[1];
  angp = _angle(xp, yp, xs, ys);

  for (jj = 2; jj < _num_pts; jj++) {
    angs = _angle(_xxin[jj], _yyin[jj], xs, ys);

    if (fabs((angs - angp)) < t) {
      /* skip */
      xp = _xxin[jj];
      yp = _yyin[jj];

    } else {
      /* store */
      _rxx[cntr] = xp;
      _ryy[cntr++] = yp;

      xs = xp;
      ys = yp;
      xp = _xxin[jj];
      yp = _yyin[jj];

      angp = _angle(xp, yp, xs, ys);
    }
  }
  /* always keep the last one */
  _rxx[cntr] = _xxin[_num_pts - 1];
  _ryy[cntr++] = _yyin[_num_pts - 1];

  return cntr;
}

double xs_trans::_angle(long int x1, long int y1, long int xs, long int ys) {
  double dx, dy, len;

  dx = (double)x1 - (double)xs;
  dy = (double)y1 - (double)ys;

  len = sqrt(dx * dx + dy * dy);
  dx /= len;
  dy /= len;
  return (acos(dy));
}

int xs_trans::_insert2segs() {
  int cntr = 0;
  int jj;

  for (jj = 1; jj < _num_pts; jj++) { // insert and
    if (_ryy[jj] <= _ryy[jj - 1]) {
      _segs[cntr]._x1 = _rxx[jj];
      _segs[cntr]._y1 = _ryy[jj];
      _segs[cntr]._x2 = _rxx[jj - 1];
      _segs[cntr]._y2 = _ryy[jj - 1];
    } else if (_ryy[jj] > _ryy[jj - 1]) {
      _segs[cntr]._x1 = _rxx[jj - 1];
      _segs[cntr]._y1 = _ryy[jj - 1];
      _segs[cntr]._x2 = _rxx[jj];
      _segs[cntr]._y2 = _ryy[jj];
    }
    cntr++;
  }

  return cntr;
}

int xs_trans::_find_intersects(long int yy, long int *working) {
  int jj, kk, ll;
  int cntr = 0, ncntr;
  double slope;

  for (jj = 0; jj < _num_segs; jj++) {
    if (_segs[jj]._y1 > yy)
      break; /* we are done ! */

    if (_segs[jj]._y1 == yy &&
        _segs[jj]._y2 == yy) { /* this is special case! */
      working[cntr++] = _segs[jj]._x1;
      working[cntr++] = _segs[jj]._x2; /* counted as two! */
    } else if (_segs[jj]._y1 < yy && _segs[jj]._y2 > yy) {
      slope = ((double)_segs[jj]._x2 - (double)_segs[jj]._x1) /
              ((double)_segs[jj]._y2 - (double)_segs[jj]._y1);
      working[cntr++] =
          _segs[jj]._x1 +
          (long int)(((double)yy - (double)_segs[jj]._y1) * slope);
    } else if (_segs[jj]._y1 == yy) {
      working[cntr++] = _segs[jj]._x1;
    } else if (_segs[jj]._y2 == yy) {
      working[cntr++] = _segs[jj]._x2;
    }
  }

  /* do a qsort and remove the redundant entries */
  qsort(&working[0], cntr, sizeof(long int), li_compare);
  for (jj = 0; jj < cntr; jj++) {
    for (kk = jj + 1; kk < cntr; kk++)
      if (working[jj] != working[kk])
        break;
    if (kk == jj + 1)
      continue;
    for (ll = jj + 1; kk < cntr; kk++, ll++)
      working[ll] = working[kk];
  }
  ncntr = 1;
  for (jj = 1; jj < cntr; jj++)
    ncntr = working[jj] == working[jj - 1] ? ncntr : ncntr + 1;

  return ncntr;
}

double xs_trans::_get_len(int num_int, long int *working, double unit) {
  int jj;
  long int sum = 0;
  if (num_int == 1)
    return 0.0;
  for (jj = 0; jj < num_int; jj += 2) {
    sum += working[jj + 1] - working[jj];
  }
  return (double)sum * unit;
}

double xs_trans::_get_pp(int num_now, long int *now, int num_prev,
                         long int *prev, double dy, double unit) {
  int jj, kk, num_tt;
  double pp = 0, left, right, mid;
  double ddy = (double)dy;
  long int ww[256]; // hard coded

  for (jj = 0; jj < num_now; jj += 2) {
    num_tt = 0;
    // find overlapping groups
    for (kk = 0; kk < num_prev; kk++) {
      if (prev[kk] >= now[jj] && prev[kk] <= now[jj + 1]) {
        ww[num_tt++] = prev[kk];
      }
    }

    if (num_tt == 0) { // new segment
      mid = ((double)now[jj] + (double)now[jj + 1]) / 2.0;
      left = (double)now[jj] - mid;
      right = mid - (double)now[jj + 1];
      pp += sqrt(left * left + ddy * ddy) + sqrt(right * right + ddy * ddy);
    } else { //  old segment
      left = (double)now[jj] - (double)ww[0];
      right = (double)now[jj + 1] - (double)ww[num_tt - 1];
      pp += sqrt(left * left + ddy * ddy) + sqrt(right * right + ddy * ddy);
    }
  }

  pp *= unit;
  return pp;
}

void xs_trans::_find_min_max() {
  int jj;
  long int mymin, mymax;

  mymin = _ryy[0];
  for (jj = 1; jj < _num_pts; jj++)
    mymin = mymin < _ryy[jj] ? mymin : _ryy[jj];

  mymax = _ryy[0] > _ryy[_num_pts - 1] ? _ryy[_num_pts - 1] : _ryy[0];

  _min_y = mymin;
  _max_y = mymax;
}

// do a quick check, returns 0 is good, -1 if failed
// situations not allowed:
//    1) there is no channel
//    2) x-coordinates are not increasing
int xs_trans::_quick_check() {
  int rc = -1;
  int jj;

  if (_num_pts < 3) {
    fprintf(stderr, "x-section %d has only %d points, need at least 3\n", _name,
            _num_pts);
    return rc;
  }

  memset(&_tmp[0], 0, _num_pts * sizeof(int));
  // x-coordiates have to be monotonically increasing, this might be relaxed in
  // the future as far as the channel doesn't completely close!
  rc = 0;
  for (jj = 1; jj < _num_pts; jj++) {
    _tmp[jj] =
        _yyin[jj] > _yyin[jj - 1] ? 1 : (_yyin[jj] < _yyin[jj - 1] ? -1 : 0);
    if (_xxin[jj] < _xxin[jj - 1]) {
      rc = -1;
      break;
    } // equal is OK
  }
  if (rc == -1) {
    fprintf(stderr, "x-section %d is not open\n", _name);
    return rc;
  }

  // find local min and local max, tag with 0 or 2
  for (jj = 1; jj < _num_pts - 1; jj++) {
    _tmp[jj] = (_yyin[jj - 1] >= _yyin[jj] && _yyin[jj + 1] >= _yyin[jj])
                   ? 0
                   : _tmp[jj];
    _tmp[jj] = (_yyin[jj - 1] <= _yyin[jj] && _yyin[jj + 1] <= _yyin[jj])
                   ? 2
                   : _tmp[jj];
  }

  // at this point, tmp stores
  //  0 - local min, 2 - local max, 1 - rising ramp, -1 - falling ramp
  // excluding the first and last entry
  rc = -1;
  for (jj = 1; jj < _num_pts - 1; jj++) {
    if (_tmp[jj] == 0) {
      rc = 0;
      break;
    }
  }

  if (rc == -1) {
    fprintf(stderr, "x-section %d does not form a channel\n", _name);
    return rc;
  }

  return 0;
}

int xs_trans::Construct(int num, double *x, double *y, double minb) {
  long int cur, ll, dyl;
  int line, num_int, num_prev_int, jj, rc;
  double len_now, len_prev, dy, c2, ccor;

  _num_pts = num;
  _record_into_int(_num_pts, &x[0], &(y[0]), UNIT);
  rc = _quick_check();
  if (rc == -1)
    return (-1);

#ifdef DBGSPLINE
  printf("original values\n");
  for (int jj = 0; jj < _num_pts; jj++) {
    printf("=1= %.4f %.4f\n", x[jj], y[jj]);
  }
  printf("converted to int\n");
  for (int jj = 0; jj < _num_pts; jj++) {
    printf("=2= %ld %ld\n", _xxin[jj], _yyin[jj]);
  }
#endif

  _num_pts = _cmethod(ANG_THRESH);

#ifdef DBGSPLINE
  printf("after compression\n");
  for (int jj = 0; jj < _num_pts; jj++) {
    printf("=3= %ld %ld\n", _rxx[jj], _ryy[jj]);
  }
#endif

  _find_min_max();
  _num_segs = _insert2segs();

  qsort(_segs, _num_segs, sizeof(segment), mycompare);

  // initialize the rest
  _init_internal();

  // find min y, if odd # of intersection, skip
  cur = _min_y;
  while (1) {
    num_int = _find_intersects(cur, &_working[0]);
    len_now = _get_len(num_int, &_working[0], UNIT);
    if (num_int % 2 == 0 && len_now > minb)
      break; // another magic number
    //    if ( num_int % 2 == 0 || num_int == 1 ) break;

    cur += YGRID;
  }

  // first value
  line = 0;
  _my_yy[line] = (double)cur * UNIT;
  _my_aa[line] = 0.0;
  _my_pp[line] = 0.0;
  _my_cc[line] = 0.0;

  for (jj = 0; jj < num_int; jj += 2) {
    _my_pp[line] += _working[jj + 1] - _working[jj];
  }
  _my_pp[line] *= UNIT;
  _my_ll[line] = _my_pp[line];

  // store the previous
  for (jj = 0; jj < num_int; jj++)
    _prev[jj] = _working[jj];
  num_prev_int = num_int;

  // get ready
  line++;
  ll = cur;
  cur += YGRID;

  // work out the rest
  for (; cur < _max_y; cur += YGRID) {
    num_int = _find_intersects(cur, &_working[0]);
    if (num_int % 2 == 1)
      continue; // skip odd number of intersections

    len_now = _get_len(num_int, &_working[0], UNIT);
    len_prev = _get_len(num_prev_int, &_prev[0], UNIT);
    dyl = cur - ll;
    dy = (double)dyl * UNIT;

    c2 = len_prev * dy * dy / 2.0 + (len_now - len_prev) * dy * dy / 6.0;
    ccor = _my_aa[line - 1] * dy;

    _my_yy[line] = (double)cur * UNIT;

    _my_aa[line] = _my_aa[line - 1] + dy * (len_now + len_prev) / 2.0;
    _my_pp[line] =
        _my_pp[line - 1] +
        _get_pp(num_int, &_working[0], num_prev_int, &_prev[0], dyl, UNIT);
    _my_ll[line] = len_now;
    _my_cc[line] = _my_cc[line - 1] + ccor + c2;

    // get ready for the next one
    line++;
    ll = cur;
    for (jj = 0; jj < num_int; jj++)
      _prev[jj] = _working[jj];
    num_prev_int = num_int;
  }
  _num_samples = line;

#ifdef DBGSPLINE
  printf("sampling dense\n");
  for (int jj = 0; jj < _num_samples; jj++) {
    printf("=4= %.4e %.4e %.4e %.4e %.4e\n", _my_aa[jj], _my_ll[jj], _my_pp[jj],
           _my_cc[jj], _my_yy[jj]);
  }
#endif
  // coarse sampling
  _final_samples = 0;
  for (jj = 0; jj < _num_samples; jj += SGRID) {
    _aa[_final_samples] = _my_aa[jj];
    _ll[_final_samples] = _my_ll[jj];
    _pp[_final_samples] = _my_pp[jj];
    _cc[_final_samples] = _my_cc[jj];
    _yy[_final_samples] = _my_yy[jj];

    _final_samples++;
  }

  // add the last one, if it's not too close
  if (jj - 2 * SGRID / 3 < _num_samples) {
    _aa[_final_samples] = _my_aa[_num_samples - 1];
    _ll[_final_samples] = _my_ll[_num_samples - 1];
    _pp[_final_samples] = _my_pp[_num_samples - 1];
    _cc[_final_samples] = _my_cc[_num_samples - 1];
    _yy[_final_samples] = _my_yy[_num_samples - 1];
    _final_samples++;
  }

#ifdef DBGSPLINE
  printf("sampling coarse\n");
  for (int jj = 0; jj < _final_samples; jj++) {
    printf("=5= %.4e %.4e %.4e %.4e %.4e\n", _aa[jj], _ll[jj], _pp[jj], _cc[jj],
           _yy[jj]);
  }
#endif

  // scale y as depth
  for (jj = 1; jj < _final_samples; jj++)
    _yy[jj] -= _yy[0];
  _yy[0] = 0.0;

  return 0;
}

// reset the container, allocate the memory
void xs_trans::Init(int cname, int sz) {
  assert(sz > 0);

  _xxin.size(sz);
  _yyin.size(sz);
  _rxx.size(sz);
  _ryy.size(sz);
  _tmp.size(sz);

  _segs.size(sz);

  _num_segs = 0;
  _num_pts = 0;
  _num_samples = 0;
  _final_samples = 0;

  _name = cname;
}

// for debuggin only
void xs_trans::PrintDenseAll(FILE *f, const char *head) {
  int jj;
  for (jj = 0; jj < _num_samples; jj++) {
    fprintf(f, "%s %.4e %.4e %.4e %.4e %.4e\n", head, _my_aa[jj], _my_yy[jj],
            _my_ll[jj], _my_cc[jj], _my_pp[jj]);
  }
}

void xs_trans::PrintSparseAll(FILE *f, const char *head) {
  int jj;
  for (jj = 0; jj < _final_samples; jj++) {
    fprintf(f, "%s %.4e %.4e %.4e %.4e %.4e\n", head, _aa[jj], _yy[jj], _ll[jj],
            _cc[jj], _pp[jj]);
  }
}

// End
