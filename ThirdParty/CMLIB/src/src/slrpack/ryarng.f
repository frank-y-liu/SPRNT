      subroutine ryarng (b, bprev, angle, prkey)
c
c-----------------------------------------------------------------------
c
c     Package:  SLRPACK
c
c     Version:  October, 1985
c
c-----------------------------------------------------------------------
c
c     PURPOSE
c     -------
c
c     This subroutine is used by user-callable subroutine RYORK to
c     rearrange the slope entries of three lines so that the slope of th
c     third line is associated with the smallest angle with respect to
c     the x-axis of line i and of line 4 for i = 1,2,3.
c
c-----------------------------------------------------------------------
c
      dimension angle(*), b(*)
      logical prkey
      pi = 4.0 * atan(1.0)
c
c      Line i -- slope b(i)  for i = 1,2,3
c
c      Line 4 -- slope bprev
c
c      Compute angles with respect to x-axis of line i for i = 1,2,3
c      and of line 4
c
      do 10 i = 1,3
           angle(i) = atan(b(i))
           if (b(i) .lt. 0.) angle(i) = pi - angle(i)
   10 continue
      angle(4) = atan(bprev)
      if (bprev .lt. 0.) angle(4) = pi - angle(4)
      if (prkey) write (6,901) (angle(i), i = 1,4)
  901 format (' the candidate regression lines form angles with y=0:',/,
     1        ' angle(1) = ', f8.3,/,
     2        ' angle(2) = ', f8.3,/,
     3        ' angle(3) = ', f8.3,/,
     4        ' and for the line selected at the last iteration: ',/,
     5        ' angle(4) = ', f8.3,//)
c
c      Compute differences between angle(i) and angle(4) for i=1,2,3
c
      do 20 i = 1,3
         angle(i) = abs(angle(i) - angle(4))
   20 continue
      if (prkey) write (6,902) (angle(i), i = 1,3)
  902 format (' the differences between angle(i) and angle(4) are:',/,
     1        ' i = 1:', f8.3,/,
     2        ' i = 2:', f8.3,/,
     3        ' i = 3:', f8.3,//)
c
c      Rearrange entries b(i) so b(3) is associated
c      with smallest angle(i)
c
      if (angle(3) .le. angle(1)) go to 30
         btemp = b(3)
         b(3) = b(1)
         b(1) = btemp
         atemp = angle(3)
         angle(3) = angle(1)
         angle(1) = atemp
   30 if (angle(3) .le. angle(2)) go to 40
         btemp = b(3)
         b(3) = b(2)
         b(2) = btemp
         atemp = angle(3)
         angle(3) = angle(2)
         angle(2) = atemp
  40  if (prkey) write (6,903) (i, angle(i), b(i), i = 1,3)
  903 format (' angles and slopes rearranged so angle(3) is smallest:',/
     1        '     i    angle(i)    b(i)',/,
     2        3(i6,f11.3,f13.3/))
      return
      end
