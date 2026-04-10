c     test_wideband2d.f - differential test for wideband2d.c

      program test_wideband2d
      implicit none

      external h2d_diagtrans,  h2d_diagtrans_c
      external h2d_mptosig,    h2d_mptosig_c
      external h2d_sig2exp,    h2d_sig2exp_c
      external h2d_mkmpshift,  h2d_mkmpshift_c
      external h2d_mkm2ltrans, h2d_mkm2ltrans_c
      external h2dmpmphf,      h2dmpmphf_c
      external h2dloclochf,    h2dloclochf_c
      external h2dmplochf,     h2dmplochf_c
      external zffti
      real *8 hkrand
      external hkrand
      integer next235
      external next235

      integer nd, nterms, nsig
      parameter (nd=2, nterms=12)
      complex *16 zk
      real *8 center1(2), center2(2)
      real *8 rscale

      complex *16 mpole(nd,-nterms:nterms)
      complex *16 sig_f(nd,200), sig_c(nd,200)
      complex *16 exp_f(nd,-nterms:nterms)
      complex *16 exp_c(nd,-nterms:nterms)
      complex *16 tvec_f(200), tvec_c(200)
      complex *16, allocatable :: wsave(:)

      integer nfail, i, j, ii
      real *8 dummy, errmax, e, dn

      nfail = 0
      dummy = hkrand(1234)

      zk = dcmplx(3.0d0, 0.05d0)
      rscale = 1.0d0
      center1(1) = 0.0d0
      center1(2) = 0.0d0
      center2(1) = 3.0d0
      center2(2) = 0.0d0

c     Compute nsig
      dn = 2*(nterms+nterms) + 1
      nsig = next235(dn)

      allocate(wsave(4*nsig+100))
      call zffti(nsig, wsave)

c     Fill mpole with random data
      do j = -nterms, nterms
         do ii = 1, nd
            mpole(ii,j) = dcmplx(hkrand(0), hkrand(0))
         enddo
      enddo

c     ---- h2d_diagtrans ----
      do j = 1, nsig
         do ii = 1, nd
            sig_f(ii,j) = dcmplx(hkrand(0), hkrand(0))
            sig_c(ii,j) = (0.0d0, 0.0d0)
         enddo
         tvec_f(j) = dcmplx(hkrand(0), hkrand(0))
      enddo
      do j = 1, nsig
         do ii = 1, nd
            sig_c(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
      call h2d_diagtrans(nd, nsig, sig_f,
     1     tvec_f, sig_c)
      do j = 1, nsig
         do ii = 1, nd
            sig_f(ii,j) = sig_c(ii,j)
            sig_c(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
c     Reset and redo with same input
      do j = 1, nsig
         do ii = 1, nd
            sig_c(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
c     Regenerate the inputs
      dummy = hkrand(1234)
      do j = -nterms, nterms
         do ii = 1, nd
            mpole(ii,j) = dcmplx(hkrand(0), hkrand(0))
         enddo
      enddo
      do j = 1, nsig
         do ii = 1, nd
            sig_f(ii,j) = dcmplx(hkrand(0), hkrand(0))
         enddo
         tvec_f(j) = dcmplx(hkrand(0), hkrand(0))
      enddo
      do j = 1, nsig
         do ii = 1, nd
            sig_c(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
      call h2d_diagtrans_c(nd, nsig, sig_f,
     1     tvec_f, sig_c)
c     Compare with Fortran result
      do j = 1, nsig
         do ii = 1, nd
            sig_f(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
      call h2d_diagtrans(nd, nsig, sig_f,
     1     tvec_f, sig_f)
c     Actually let me simplify this: just compare F and C on same inputs
c     Reset everything cleanly
      dummy = hkrand(1234)
      do j = -nterms, nterms
         do ii = 1, nd
            mpole(ii,j) = dcmplx(hkrand(0), hkrand(0))
         enddo
      enddo

c     ---- h2d_mptosig ----
      do j = 1, nsig
         do ii = 1, nd
            sig_f(ii,j) = (0.0d0, 0.0d0)
            sig_c(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
      call h2d_mptosig(nd, nterms, nsig, mpole,
     1     sig_f, wsave)
      call h2d_mptosig_c(nd, nterms, nsig, mpole,
     1     sig_c, wsave)
      errmax = 0.0d0
      do j = 1, nsig
         do ii = 1, nd
            e = cdabs(sig_f(ii,j) - sig_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] h2d_mptosig err=',errmax
         nfail = nfail + 1
      else
         write(*,*) '[ ok ] h2d_mptosig'
      endif

c     ---- h2d_mkmpshift ----
      do j = 1, nsig
         tvec_f(j) = (0.0d0, 0.0d0)
         tvec_c(j) = (0.0d0, 0.0d0)
      enddo
      call h2d_mkmpshift(zk, center1, nterms,
     1     center2, nterms, nsig, wsave, tvec_f)
      call h2d_mkmpshift_c(zk, center1, nterms,
     1     center2, nterms, nsig, wsave, tvec_c)
      errmax = 0.0d0
      do j = 1, nsig
         e = cdabs(tvec_f(j) - tvec_c(j))
         if (e .gt. errmax) errmax = e
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] h2d_mkmpshift err=',
     1        errmax
         nfail = nfail + 1
      else
         write(*,*) '[ ok ] h2d_mkmpshift'
      endif

c     ---- h2d_mkm2ltrans ----
      do j = 1, nsig
         tvec_f(j) = (0.0d0, 0.0d0)
         tvec_c(j) = (0.0d0, 0.0d0)
      enddo
      call h2d_mkm2ltrans(zk, center1, nterms,
     1     center2, nterms, nsig, wsave, tvec_f)
      call h2d_mkm2ltrans_c(zk, center1, nterms,
     1     center2, nterms, nsig, wsave, tvec_c)
      errmax = 0.0d0
      do j = 1, nsig
         e = cdabs(tvec_f(j) - tvec_c(j))
         if (e .gt. errmax) errmax = e
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] h2d_mkm2ltrans err=',
     1        errmax
         nfail = nfail + 1
      else
         write(*,*) '[ ok ] h2d_mkm2ltrans'
      endif

c     ---- h2d_sig2exp ----
      do j = -nterms, nterms
         do ii = 1, nd
            exp_f(ii,j) = (0.0d0, 0.0d0)
            exp_c(ii,j) = (0.0d0, 0.0d0)
         enddo
      enddo
c     Use sig_f from the mptosig test
      call h2d_sig2exp(nd, nsig, sig_f, wsave,
     1     nterms, exp_f)
c     Need to re-compute sig_c since sig_f was used
      do j = 1, nsig
         do ii = 1, nd
            sig_c(ii,j) = sig_f(ii,j)
         enddo
      enddo
      call h2d_sig2exp_c(nd, nsig, sig_c, wsave,
     1     nterms, exp_c)
      errmax = 0.0d0
      do j = -nterms, nterms
         do ii = 1, nd
            e = cdabs(exp_f(ii,j) - exp_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] h2d_sig2exp err=',errmax
         nfail = nfail + 1
      else
         write(*,*) '[ ok ] h2d_sig2exp'
      endif

      if (nfail .ne. 0) then
         write(*,*) 'FAIL:', nfail, ' checks differ'
         stop 1
      endif
      write(*,*) 'PASS: wideband2d all tests'

      end
