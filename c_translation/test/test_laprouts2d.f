c     test_laprouts2d.f - differential test for the C translation of
c     src/laplace/laprouts2d.f.
c
c     For each routine, calls the Fortran reference and the _c_
c     translation on identical inputs into separate output buffers
c     and asserts that every complex element agrees in absolute
c     value. With both sides built at -O0 and the C source preserving
c     the Fortran's left-to-right floating-point operation order, the
c     expected error is exactly zero (we allow a tiny tolerance).

      program test_laprouts2d
      implicit none

      external l2dformmpc,  l2dformmpc_c
      external l2dformmpd,  l2dformmpd_c
      external l2dformmpcd, l2dformmpcd_c
      external l2dmpevalp,  l2dmpevalp_c
      external l2dmpevalg,  l2dmpevalg_c
      external l2dmpevalh,  l2dmpevalh_c
      external l2dformtac,  l2dformtac_c
      external l2dformtad,  l2dformtad_c
      external l2dformtacd, l2dformtacd_c
      external l2dtaevalp,  l2dtaevalp_c
      external l2dtaevalg,  l2dtaevalg_c
      external l2dtaevalh,  l2dtaevalh_c
      external l2dmpmp,     l2dmpmp_c
      external l2dlocloc,   l2dlocloc_c
      external l2dmploc,    l2dmploc_c
      external l2dmpzero,   l2dmpzero_c
      external init_carray
      real *8 hkrand
      external hkrand

      integer nd, nterms, ns, nt
      integer nterms1, nterms2, ldc
      parameter (nd = 2)
      parameter (nterms = 20)
      parameter (ns = 30)
      parameter (nt = 25)
      parameter (nterms1 = 18)
      parameter (nterms2 = 22)
      parameter (ldc = 44)

      real *8 rscale, rscale1, rscale2
      real *8 source(2, ns), targ(2, nt)
      real *8 center(2), center1(2), center2(2)
      complex *16 charge(nd, ns), dipstr(nd, ns)

      complex *16 mpole_f(nd, 0:nterms),  mpole_c(nd, 0:nterms)
      complex *16 local_f(nd, 0:nterms),  local_c(nd, 0:nterms)
      complex *16 pot_f(nd, nt),          pot_c(nd, nt)
      complex *16 grad_f(nd, nt),         grad_c(nd, nt)
      complex *16 hess_f(nd, nt),         hess_c(nd, nt)

      complex *16 hexp1(nd, 0:nterms1), hexp1c(nd, 0:nterms1)
      complex *16 hexp2_f(nd, 0:nterms2), hexp2_c(nd, 0:nterms2)
      complex *16 jexp1(nd, 0:nterms1), jexp1c(nd, 0:nterms1)
      complex *16 jexp2_f(nd, 0:nterms2), jexp2_c(nd, 0:nterms2)

      complex *16 mpzero_f(nd, 0:nterms), mpzero_c(nd, 0:nterms)

      real *8 carray(0:ldc, 0:ldc)

      integer i, j, ii, n, nfail
      real *8 dummy, errmax, e

      nfail = 0
      rscale = 0.5d0
      rscale1 = 0.5d0
      rscale2 = 0.6d0

      center(1) = 0.0d0
      center(2) = 0.0d0
      center1(1) = 0.1d0
      center1(2) = 0.2d0
      center2(1) = 0.4d0
      center2(2) = 0.5d0

c     seed hkrand once
      dummy = hkrand(1234)

c     sources placed near the multipole expansion center, well within
c     the convergence region (offset 0.1*hkrand)
      do i = 1, ns
         source(1, i) = center(1) + 0.1d0*hkrand(0)
         source(2, i) = center(2) + 0.1d0*hkrand(0)
      enddo
c     targets near the local expansion center, similarly small offsets
      do j = 1, nt
         targ(1, j) = center(1) + 0.1d0*hkrand(0)
         targ(2, j) = center(2) + 0.1d0*hkrand(0)
      enddo
      do i = 1, ns
         do ii = 1, nd
            charge(ii, i) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                             2.0d0*hkrand(0)-1.0d0)
            dipstr(ii, i) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                             2.0d0*hkrand(0)-1.0d0)
         enddo
      enddo

c     fill carray with binomial coefficients
      call init_carray(carray, ldc)

c     fill hexp1 / jexp1 with random complex data
      do n = 0, nterms1
         do ii = 1, nd
            hexp1(ii, n) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                            2.0d0*hkrand(0)-1.0d0)
            hexp1c(ii, n) = hexp1(ii, n)
            jexp1(ii, n) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                            2.0d0*hkrand(0)-1.0d0)
            jexp1c(ii, n) = jexp1(ii, n)
         enddo
      enddo

c     ============================================================
c     l2dformmpc
c     ============================================================
      call zerocplx(mpole_f, nd*(nterms+1))
      call zerocplx(mpole_c, nd*(nterms+1))
      call l2dformmpc  (nd, rscale, source, ns, charge, center,
     1                  nterms, mpole_f)
      call l2dformmpc_c(nd, rscale, source, ns, charge, center,
     1                  nterms, mpole_c)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(mpole_f(ii,n) - mpole_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dformmpc  ', errmax, nfail)

c     ============================================================
c     l2dformmpd
c     ============================================================
      call zerocplx(mpole_f, nd*(nterms+1))
      call zerocplx(mpole_c, nd*(nterms+1))
      call l2dformmpd  (nd, rscale, source, ns, dipstr, center,
     1                  nterms, mpole_f)
      call l2dformmpd_c(nd, rscale, source, ns, dipstr, center,
     1                  nterms, mpole_c)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(mpole_f(ii,n) - mpole_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dformmpd  ', errmax, nfail)

c     ============================================================
c     l2dformmpcd
c     ============================================================
      call zerocplx(mpole_f, nd*(nterms+1))
      call zerocplx(mpole_c, nd*(nterms+1))
      call l2dformmpcd  (nd, rscale, source, ns, charge, dipstr,
     1                   center, nterms, mpole_f)
      call l2dformmpcd_c(nd, rscale, source, ns, charge, dipstr,
     1                   center, nterms, mpole_c)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(mpole_f(ii,n) - mpole_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dformmpcd ', errmax, nfail)

c     ============================================================
c     l2dmpevalp / l2dmpevalg / l2dmpevalh
c
c     Use the multipole built by l2dformmpcd above (mpole_f) as input.
c     Targets must be OUTSIDE the source ring; we shift the
c     evaluation targets to live well outside.
c     ============================================================
c     build a non-trivial mpole to evaluate from
      call zerocplx(mpole_f, nd*(nterms+1))
      call l2dformmpcd(nd, rscale, source, ns, charge, dipstr,
     1                 center, nterms, mpole_f)

c     re-place targets outside the source region for multipole eval
      do j = 1, nt
         targ(1, j) = center(1) + 5.0d0 + hkrand(0)
         targ(2, j) = center(2) + 5.0d0 + hkrand(0)
      enddo

      call zerocplx(pot_f, nd*nt)
      call zerocplx(pot_c, nd*nt)
      call l2dmpevalp  (nd, rscale, center, mpole_f, nterms,
     1                  targ, nt, pot_f)
      call l2dmpevalp_c(nd, rscale, center, mpole_f, nterms,
     1                  targ, nt, pot_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pot_f(ii,j) - pot_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dmpevalp  ', errmax, nfail)

      call zerocplx(pot_f,  nd*nt)
      call zerocplx(pot_c,  nd*nt)
      call zerocplx(grad_f, nd*nt)
      call zerocplx(grad_c, nd*nt)
      call l2dmpevalg  (nd, rscale, center, mpole_f, nterms,
     1                  targ, nt, pot_f, grad_f)
      call l2dmpevalg_c(nd, rscale, center, mpole_f, nterms,
     1                  targ, nt, pot_c, grad_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pot_f(ii,j) - pot_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(grad_f(ii,j) - grad_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dmpevalg  ', errmax, nfail)

      call zerocplx(pot_f,  nd*nt)
      call zerocplx(pot_c,  nd*nt)
      call zerocplx(grad_f, nd*nt)
      call zerocplx(grad_c, nd*nt)
      call zerocplx(hess_f, nd*nt)
      call zerocplx(hess_c, nd*nt)
      call l2dmpevalh  (nd, rscale, center, mpole_f, nterms,
     1                  targ, nt, pot_f, grad_f, hess_f)
      call l2dmpevalh_c(nd, rscale, center, mpole_f, nterms,
     1                  targ, nt, pot_c, grad_c, hess_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pot_f(ii,j) - pot_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(grad_f(ii,j) - grad_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(hess_f(ii,j) - hess_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dmpevalh  ', errmax, nfail)

c     restore targets to live near center for the local expansion
c     form/eval tests below
      do j = 1, nt
         targ(1, j) = center(1) + 0.1d0*hkrand(0)
         targ(2, j) = center(2) + 0.1d0*hkrand(0)
      enddo
c     also restore sources to be FAR away (outside) for the local
c     expansion form (sources need to be away from the local center)
      do i = 1, ns
         source(1, i) = center(1) + 5.0d0 + hkrand(0)
         source(2, i) = center(2) + 5.0d0 + hkrand(0)
      enddo

c     ============================================================
c     l2dformtac
c     ============================================================
      call zerocplx(local_f, nd*(nterms+1))
      call zerocplx(local_c, nd*(nterms+1))
      call l2dformtac  (nd, rscale, source, ns, charge, center,
     1                  nterms, local_f)
      call l2dformtac_c(nd, rscale, source, ns, charge, center,
     1                  nterms, local_c)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(local_f(ii,n) - local_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dformtac  ', errmax, nfail)

c     ============================================================
c     l2dformtad
c     ============================================================
      call zerocplx(local_f, nd*(nterms+1))
      call zerocplx(local_c, nd*(nterms+1))
      call l2dformtad  (nd, rscale, source, ns, dipstr, center,
     1                  nterms, local_f)
      call l2dformtad_c(nd, rscale, source, ns, dipstr, center,
     1                  nterms, local_c)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(local_f(ii,n) - local_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dformtad  ', errmax, nfail)

c     ============================================================
c     l2dformtacd
c     ============================================================
      call zerocplx(local_f, nd*(nterms+1))
      call zerocplx(local_c, nd*(nterms+1))
      call l2dformtacd  (nd, rscale, source, ns, charge, dipstr,
     1                   center, nterms, local_f)
      call l2dformtacd_c(nd, rscale, source, ns, charge, dipstr,
     1                   center, nterms, local_c)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(local_f(ii,n) - local_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dformtacd ', errmax, nfail)

c     ============================================================
c     l2dtaevalp / g / h
c     Use the local expansion just built
c     ============================================================
      call zerocplx(local_f, nd*(nterms+1))
      call l2dformtacd(nd, rscale, source, ns, charge, dipstr,
     1                 center, nterms, local_f)

      call zerocplx(pot_f, nd*nt)
      call zerocplx(pot_c, nd*nt)
      call l2dtaevalp  (nd, rscale, center, local_f, nterms,
     1                  targ, nt, pot_f)
      call l2dtaevalp_c(nd, rscale, center, local_f, nterms,
     1                  targ, nt, pot_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pot_f(ii,j) - pot_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dtaevalp  ', errmax, nfail)

      call zerocplx(pot_f,  nd*nt)
      call zerocplx(pot_c,  nd*nt)
      call zerocplx(grad_f, nd*nt)
      call zerocplx(grad_c, nd*nt)
      call l2dtaevalg  (nd, rscale, center, local_f, nterms,
     1                  targ, nt, pot_f, grad_f)
      call l2dtaevalg_c(nd, rscale, center, local_f, nterms,
     1                  targ, nt, pot_c, grad_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pot_f(ii,j) - pot_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(grad_f(ii,j) - grad_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dtaevalg  ', errmax, nfail)

      call zerocplx(pot_f,  nd*nt)
      call zerocplx(pot_c,  nd*nt)
      call zerocplx(grad_f, nd*nt)
      call zerocplx(grad_c, nd*nt)
      call zerocplx(hess_f, nd*nt)
      call zerocplx(hess_c, nd*nt)
      call l2dtaevalh  (nd, rscale, center, local_f, nterms,
     1                  targ, nt, pot_f, grad_f, hess_f)
      call l2dtaevalh_c(nd, rscale, center, local_f, nterms,
     1                  targ, nt, pot_c, grad_c, hess_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pot_f(ii,j) - pot_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(grad_f(ii,j) - grad_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(hess_f(ii,j) - hess_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dtaevalh  ', errmax, nfail)

c     ============================================================
c     l2dmpmp
c     ============================================================
      call zerocplx(hexp2_f, nd*(nterms2+1))
      call zerocplx(hexp2_c, nd*(nterms2+1))
      call l2dmpmp  (nd, rscale1, center1, hexp1, nterms1,
     1               rscale2, center2, hexp2_f, nterms2, carray, ldc)
      call l2dmpmp_c(nd, rscale1, center1, hexp1, nterms1,
     1               rscale2, center2, hexp2_c, nterms2, carray, ldc)
      errmax = 0.0d0
      do n = 0, nterms2
         do ii = 1, nd
            e = cdabs(hexp2_f(ii,n) - hexp2_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dmpmp     ', errmax, nfail)

c     ============================================================
c     l2dlocloc
c     ============================================================
      call zerocplx(jexp2_f, nd*(nterms2+1))
      call zerocplx(jexp2_c, nd*(nterms2+1))
      call l2dlocloc  (nd, rscale1, center1, jexp1, nterms1,
     1                 rscale2, center2, jexp2_f, nterms2, carray, ldc)
      call l2dlocloc_c(nd, rscale1, center1, jexp1, nterms1,
     1                 rscale2, center2, jexp2_c, nterms2, carray, ldc)
      errmax = 0.0d0
      do n = 0, nterms2
         do ii = 1, nd
            e = cdabs(jexp2_f(ii,n) - jexp2_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dlocloc   ', errmax, nfail)

c     ============================================================
c     l2dmploc
c     ============================================================
      call zerocplx(jexp2_f, nd*(nterms2+1))
      call zerocplx(jexp2_c, nd*(nterms2+1))
      call l2dmploc  (nd, rscale1, center1, hexp1, nterms1,
     1                rscale2, center2, jexp2_f, nterms2, carray, ldc)
      call l2dmploc_c(nd, rscale1, center1, hexp1, nterms1,
     1                rscale2, center2, jexp2_c, nterms2, carray, ldc)
      errmax = 0.0d0
      do n = 0, nterms2
         do ii = 1, nd
            e = cdabs(jexp2_f(ii,n) - jexp2_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dmploc    ', errmax, nfail)

c     ============================================================
c     l2dmpzero
c     Fill with non-zero data, call zero, ensure all zero
c     ============================================================
      do n = 0, nterms
         do ii = 1, nd
            mpzero_f(ii,n) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                              2.0d0*hkrand(0)-1.0d0)
            mpzero_c(ii,n) = mpzero_f(ii,n)
         enddo
      enddo
      call l2dmpzero  (nd, mpzero_f, nterms)
      call l2dmpzero_c(nd, mpzero_c, nterms)
      errmax = 0.0d0
      do n = 0, nterms
         do ii = 1, nd
            e = cdabs(mpzero_f(ii,n) - mpzero_c(ii,n))
            if (e .gt. errmax) errmax = e
            e = cdabs(mpzero_c(ii,n))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('l2dmpzero   ', errmax, nfail)


      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all laprouts2d cases match'

      end


c     ----------------------------------------------------------------
      subroutine zerocplx(z, n)
      implicit none
      integer n, i
      complex *16 z(n)
      do i = 1, n
         z(i) = dcmplx(0.0d0, 0.0d0)
      enddo
      return
      end


c     ----------------------------------------------------------------
      subroutine report(name, errmax, nfail)
      implicit none
      character*(*) name
      integer nfail
      real *8 errmax, tol
      tol = 1.0d-15
      if (errmax .lt. tol) then
         write(*,1000) name, errmax
      else
         write(*,1001) name, errmax
         nfail = nfail + 1
      endif
 1000 format(' [ ok ] ',a12,'  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a12,'  errmax=',1pe11.3)
      return
      end
