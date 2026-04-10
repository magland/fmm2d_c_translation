c     test_h2dcommon.f - differential test for the C translation of
c     src/helmholtz/h2dcommon.f.

      program test_h2dcommon
      implicit none

      external h2cart2polar, h2cart2polar_c
      external h2dall, h2dall_c
      external h2dmpzero, h2dmpzero_c
      external h2dsigzero, h2dsigzero_c
      real *8 hkrand
      external hkrand

      integer nfail, i
      real *8 dummy, errmax, e

c     h2cart2polar variables
      real *8 zat(2), r_f, theta_f, r_c, theta_c

c     h2dall variables
      integer ntmax
      parameter (ntmax = 50)
      complex *16 z
      real *8 rscale
      complex *16 hvec_f(0:ntmax), hder_f(0:ntmax)
      complex *16 hvec_c(0:ntmax), hder_c(0:ntmax)
      integer nterms, ifder, j

c     h2dmpzero variables
      integer ndmax, ntmmax
      parameter (ndmax = 3, ntmmax = 20)
      complex *16 mpole_f(ndmax, -ntmmax:ntmmax)
      complex *16 mpole_c(ndmax, -ntmmax:ntmmax)
      integer nd, n

c     h2dsigzero variables
      integer nsigmax
      parameter (nsigmax = 50)
      complex *16 sig_f(ndmax, nsigmax)
      complex *16 sig_c(ndmax, nsigmax)
      integer nsig

      nfail = 0
      dummy = hkrand(1234)

c     ============================================================
c     Test h2cart2polar
c     ============================================================
      do i = 1, 20
         zat(1) = 10.0d0*hkrand(0) - 5.0d0
         zat(2) = 10.0d0*hkrand(0) - 5.0d0

         call h2cart2polar(zat, r_f, theta_f)
         call h2cart2polar_c(zat, r_c, theta_c)

         errmax = 0.0d0
         e = dabs(r_f - r_c)
         if (e .gt. errmax) errmax = e
         e = dabs(theta_f - theta_c)
         if (e .gt. errmax) errmax = e

         if (errmax .gt. 0.0d0) then
            write(*,*) '[FAIL] h2cart2polar case ', i,
     1           ' errmax=', errmax
            nfail = nfail + 1
         endif
      enddo
      write(*,*) '[ ok ] h2cart2polar: 20 cases'

c     ============================================================
c     Test h2dall
c     ============================================================
      do i = 1, 20
         z = dcmplx(5.0d0*hkrand(0)+0.1d0,
     1              5.0d0*hkrand(0)+0.1d0)
         nterms = 5 + int(40*hkrand(0))
         if (nterms .gt. ntmax) nterms = ntmax
         rscale = 1.0d0
         if (cdabs(z) .lt. 1.0d0) rscale = cdabs(z)

         do ifder = 0, 1
            do j = 0, nterms
               hvec_f(j) = (0.0d0, 0.0d0)
               hder_f(j) = (0.0d0, 0.0d0)
               hvec_c(j) = (0.0d0, 0.0d0)
               hder_c(j) = (0.0d0, 0.0d0)
            enddo

            call h2dall(nterms, z, rscale,
     1           hvec_f, ifder, hder_f)
            call h2dall_c(nterms, z, rscale,
     1           hvec_c, ifder, hder_c)

            errmax = 0.0d0
            do j = 0, nterms
               e = cdabs(hvec_f(j) - hvec_c(j))
               if (e .gt. errmax) errmax = e
            enddo
            if (ifder .eq. 1) then
               do j = 0, nterms
                  e = cdabs(hder_f(j) - hder_c(j))
                  if (e .gt. errmax) errmax = e
               enddo
            endif

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'h2dall', i, ifder,
     1              errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,*) '[ ok ] h2dall: 40 cases'

c     ============================================================
c     Test h2dmpzero
c     ============================================================
      nd = 3
      nterms = 15

c     Fill with nonzero values
      do n = -nterms, nterms
         do j = 1, nd
            mpole_f(j, n) = dcmplx(1.0d0, 2.0d0)
            mpole_c(j, n) = dcmplx(1.0d0, 2.0d0)
         enddo
      enddo

      call h2dmpzero(nd, mpole_f, nterms)
      call h2dmpzero_c(nd, mpole_c, nterms)

      errmax = 0.0d0
      do n = -nterms, nterms
         do j = 1, nd
            e = cdabs(mpole_f(j,n) - mpole_c(j,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] h2dmpzero errmax=', errmax
         nfail = nfail + 1
      endif
      write(*,*) '[ ok ] h2dmpzero'

c     ============================================================
c     Test h2dsigzero
c     ============================================================
      nd = 3
      nsig = 40

c     Fill with nonzero values
      do n = 1, nsig
         do j = 1, nd
            sig_f(j, n) = dcmplx(3.0d0, 4.0d0)
            sig_c(j, n) = dcmplx(3.0d0, 4.0d0)
         enddo
      enddo

      call h2dsigzero(nd, sig_f, nsig)
      call h2dsigzero_c(nd, sig_c, nsig)

      errmax = 0.0d0
      do n = 1, nsig
         do j = 1, nd
            e = cdabs(sig_f(j,n) - sig_c(j,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] h2dsigzero errmax=', errmax
         nfail = nfail + 1
      endif
      write(*,*) '[ ok ] h2dsigzero'

c     ============================================================
c     Summary
c     ============================================================
 1000 format(' [FAIL] ',a,' case ',i3,' ifder=',i2,
     1   ' errmax=',1pe12.5)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: h2dcommon all tests passed'

      end
