c     test_h2dterms.f - differential test for the C translation of
c     src/helmholtz/h2dterms.f (h2dterms only).

      program test_h2dterms
      implicit none

      external h2dterms, h2dterms_c
      real *8 hkrand
      external hkrand

      integer nfail, i
      real *8 dummy
      real *8 bsize, eps
      complex *16 zk
      integer nterms_f, ier_f, nterms_c, ier_c

      nfail = 0
      dummy = hkrand(1234)

c     Test various box sizes, zk values, and eps values
      do i = 1, 20
         bsize = 0.1d0 + 5.0d0*hkrand(0)
         zk = dcmplx(1.0d0 + 10.0d0*hkrand(0), 0.0d0)
         eps = 10.0d0**(-3.0d0 - 10.0d0*hkrand(0))

         nterms_f = -999
         ier_f = -999
         nterms_c = -999
         ier_c = -999

         call h2dterms(bsize, zk, eps, nterms_f, ier_f)
         call h2dterms_c(bsize, zk, eps,
     1        nterms_c, ier_c)

         if (nterms_f .ne. nterms_c .or.
     1       ier_f .ne. ier_c) then
            write(*,1000) i, nterms_f, nterms_c,
     1           ier_f, ier_c
            nfail = nfail + 1
         endif
      enddo

 1000 format(' [FAIL] h2dterms case ',i3,
     1   ' nterms f/c=',i6,'/',i6,
     2   ' ier f/c=',i3,'/',i3)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: h2dterms all 20 cases match'

      end
