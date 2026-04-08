c     test_bh2dterms.f - differential test for the C translation of
c     subroutine bh2dterms.
c
c     Calls both bh2dterms (Fortran reference, from libfmm2d.a) and
c     bh2dterms_c_ (C translation, from build/bh2dterms.o) on a battery
c     of eps values and verifies that (nterms, ier) match exactly.

      program test_bh2dterms
      implicit none

      external bh2dterms
      external bh2dterms_c

      integer ntests
      parameter (ntests = 8)

      real *8 eps_list(ntests)
      integer i
      integer nterms_f, ier_f, nterms_c, ier_c
      integer nfail

      data eps_list /
     1   1.0d-2,   1.0d-4,   1.0d-6,   5.0d-7,
     2   1.0d-9,   2.5d-10,  1.0d-12,  1.0d-14 /

      nfail = 0

      do i = 1, ntests
         nterms_f = -999
         ier_f    = -999
         nterms_c = -999
         ier_c    = -999

         call bh2dterms(eps_list(i), nterms_f, ier_f)
         call bh2dterms_c(eps_list(i), nterms_c, ier_c)

         if (nterms_f .ne. nterms_c .or. ier_f .ne. ier_c) then
            write(*,1000) i, eps_list(i),
     1           nterms_f, nterms_c, ier_f, ier_c
            nfail = nfail + 1
         else
            write(*,1001) i, eps_list(i), nterms_f, ier_f
         endif
      enddo

 1000 format(' [FAIL] case ',i3,' eps=',1pe10.3,
     1   '  nterms f/c=',i6,'/',i6,'  ier f/c=',i3,'/',i3)
 1001 format(' [ ok ] case ',i3,' eps=',1pe10.3,
     1   '  nterms=',i6,'  ier=',i3)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' of ', ntests, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: all ', ntests, ' cases match'

      end
