c
c   driver for testing cmlib routines
c     adapt
c
c    one input data card is required
c         read(lin,1) kprint,times
c    1    format(i1,e10.0)
c
c     kprint = 0   no printing
c              1   no printing for passed tests, short message
c                  for failed tests
c              2   print short message for passed tests, fuller
c                  information for failed tests
c              3   print complete quick-check results
c
c                ***important note***
c         all quick checks use routines r2mach and d2mach
c         to set the error tolerances.
c     times is a constant multiplier that can be used to scale the
c     values of r1mach and d1mach so that
c               r2mach(i) = r1mach(i) * times   for i=3,4,5
c               d2mach(i) = d1mach(i) * times   for i=3,4,5
c     this makes it easily possible to change the error tolerances
c     used in the quick checks.
c     if times .le. 0.0 then times is defaulted to 1.0
c
c              ***end note***
c
      common/unit/lun
      common/msg/icnt,jtest(38)
      common/xxmult/times
      lun=i1mach(2)
      lin=i1mach(1)
      itest=1
c
c     read kprint,times parameters from data card..
c
      read(lin,1) kprint,times
1     format(i1,e10.0)
      if(times.le.0.) times=1.
      call xsetun(lun)
      call xsetf(1)
      call xermax(1000)
c   test adapt
      call adapqx(kprint,ipass)
      itest=itest*ipass
c
      if(kprint.ge.1.and.itest.ne.1) write(lun,2)
2     format(/' ***** warning -- at least one test for sublibrary adapt
     1  has failed ***** ')
      if(kprint.ge.1.and.itest.eq.1) write(lun,3)
3     format(/' ----- sublibrary adapt passed all tests ----- ')
      end
      double precision function d2mach(i)
      double precision d1mach
      common/xxmult/times
      d2mach=d1mach(i)
      if(i.eq.1.or. i.eq.2) return
      d2mach = d2mach * dble(times)
      end
      real function r2mach(i)
      common/xxmult/times
      r2mach=r1mach(i)
      if(i.eq.1.or. i.eq.2) return
      r2mach = r2mach * times
      return
      end
      subroutine adapqx(kprint,ipass)
      real acc,eps,finval
      integer ifail,k,maxpts,minpts,ndim,ipass,kprint
      real a(20),b(20),wrkstr(3000)
      common/unit/lun
      external functn
      ndim=2
      do 20 k=1,ndim
        a(k)=-1
        b(k)=1
   20 continue
      eps=6.e-4
      maxpts=3000
   40 minpts=100
      call  adapt(ndim,a,b,minpts,maxpts,functn,eps,acc,3000,
     * wrkstr,finval,ifail)
      aerr = abs(finval-6.28267)
      if(ifail.eq.0) then
       if(aerr.lt.amax1(eps,5.e-4) .and. acc.le.2.e-03)then
         ipass=1
         if (kprint .ge. 2) write(lun,99997) eps,finval,acc,aerr,minpts
         return
       endif
      endif
      write(lun,*)' test failing'
      ipass=2
      if(kprint.ge.2) write(lun,99998) ifail
      if(kprint.ge.2) write(lun,99997) eps,finval,acc,aerr,minpts
      return
99998 format(5x,'error ',i5)
99997 format(5x,'requested accuracy = ',e12.3/5x,'estimated ',
     * 'value    = ',e12.6/5x,'estimated accuracy = ',e12.3/
     * 8x,'actual accuracy = ',e12.3/
     * 8x,'integrand calls = ',i6)
      end
      function functn(ndim,x)
      dimension x(ndim)
      one=1
      prod=1
      sum=1
      do 10 j=2,ndim
       sum=sum*(one-x(j-1)**2)
   10  prod=prod*sum
       functn=sqrt(prod/(sum*(one-x(ndim)**2)))
      return
      end
