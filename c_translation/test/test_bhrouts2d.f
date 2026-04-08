c     test_bhrouts2d.f - differential test for the C translation of
c     src/biharmonic/bhrouts2d.f.
c
c     For each of the 14 routines, calls the Fortran reference and the
c     _c_ translation on identical inputs into separate output buffers
c     and asserts that every complex element agrees in absolute value.
c     With both sides built at -O0 and the C source preserving the
c     Fortran's left-to-right floating-point operation order, the
c     expected error is exactly zero (we allow a tiny tolerance).

      program test_bhrouts2d
      implicit none

      external bh2dformmpc,  bh2dformmpc_c
      external bh2dformmpd,  bh2dformmpd_c
      external bh2dformmpcd, bh2dformmpcd_c
      external bh2dformtac,  bh2dformtac_c
      external bh2dformtad,  bh2dformtad_c
      external bh2dformtacd, bh2dformtacd_c
      external bh2dmpevalp,  bh2dmpevalp_c
      external bh2dmpevalg,  bh2dmpevalg_c
      external bh2dtaevalp,  bh2dtaevalp_c
      external bh2dtaevalg,  bh2dtaevalg_c
      external bh2dmpmp,     bh2dmpmp_c
      external bh2dlocloc,   bh2dlocloc_c
      external bh2dmploc,    bh2dmploc_c
      external bh2dmpzero,   bh2dmpzero_c
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
c     biharmonic charge has shape (nd,2,*) and dipole shape (nd,3,*)
      complex *16 charge(nd, 2, ns), dip(nd, 3, ns)

      complex *16 mpole_f(nd, 5, 0:nterms),  mpole_c(nd, 5, 0:nterms)
      complex *16 local_f(nd, 5, 0:nterms),  local_c(nd, 5, 0:nterms)
      complex *16 vel_f(nd, nt),             vel_c(nd, nt)
      complex *16 grad_f(nd, 3, nt),         grad_c(nd, 3, nt)

      complex *16 hexp1(nd, 5, 0:nterms1)
      complex *16 hexp2_f(nd, 5, 0:nterms2)
      complex *16 hexp2_c(nd, 5, 0:nterms2)
      complex *16 jexp1(nd, 5, 0:nterms1)
      complex *16 jexp2_f(nd, 5, 0:nterms2)
      complex *16 jexp2_c(nd, 5, 0:nterms2)

      complex *16 mpzero_f(nd, 5, 0:nterms)
      complex *16 mpzero_c(nd, 5, 0:nterms)

      real *8 carray(0:ldc, 0:ldc)

      integer i, j, ii, kk, n, nfail
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

      dummy = hkrand(1234)

c     Sources placed near the multipole expansion center.
      do i = 1, ns
         source(1, i) = center(1) + 0.1d0*hkrand(0)
         source(2, i) = center(2) + 0.1d0*hkrand(0)
      enddo
      do j = 1, nt
         targ(1, j) = center(1) + 0.1d0*hkrand(0)
         targ(2, j) = center(2) + 0.1d0*hkrand(0)
      enddo
      do i = 1, ns
         do ii = 1, nd
            do kk = 1, 2
               charge(ii, kk, i) =
     1            dcmplx(2.0d0*hkrand(0)-1.0d0,
     2                   2.0d0*hkrand(0)-1.0d0)
            enddo
            do kk = 1, 3
               dip(ii, kk, i) =
     1            dcmplx(2.0d0*hkrand(0)-1.0d0,
     2                   2.0d0*hkrand(0)-1.0d0)
            enddo
         enddo
      enddo

      call init_carray(carray, ldc)

c     Fill hexp1 / jexp1 with random data.
      do n = 0, nterms1
         do kk = 1, 5
            do ii = 1, nd
               hexp1(ii, kk, n) =
     1            dcmplx(2.0d0*hkrand(0)-1.0d0,
     2                   2.0d0*hkrand(0)-1.0d0)
               jexp1(ii, kk, n) =
     1            dcmplx(2.0d0*hkrand(0)-1.0d0,
     2                   2.0d0*hkrand(0)-1.0d0)
            enddo
         enddo
      enddo

c     ============================================================
c     bh2dformmpc
c     ============================================================
      call zerobh(mpole_f, nd, nterms)
      call zerobh(mpole_c, nd, nterms)
      call bh2dformmpc  (nd, rscale, source, ns, charge, center,
     1                   nterms, mpole_f)
      call bh2dformmpc_c(nd, rscale, source, ns, charge, center,
     1                   nterms, mpole_c)
      call bhdiff('bh2dformmpc ', mpole_f, mpole_c,
     1            nd, nterms, nfail)

c     ============================================================
c     bh2dformmpd
c     ============================================================
      call zerobh(mpole_f, nd, nterms)
      call zerobh(mpole_c, nd, nterms)
      call bh2dformmpd  (nd, rscale, source, ns, dip, center,
     1                   nterms, mpole_f)
      call bh2dformmpd_c(nd, rscale, source, ns, dip, center,
     1                   nterms, mpole_c)
      call bhdiff('bh2dformmpd ', mpole_f, mpole_c,
     1            nd, nterms, nfail)

c     ============================================================
c     bh2dformmpcd
c     ============================================================
      call zerobh(mpole_f, nd, nterms)
      call zerobh(mpole_c, nd, nterms)
      call bh2dformmpcd  (nd, rscale, source, ns, charge, dip,
     1                    center, nterms, mpole_f)
      call bh2dformmpcd_c(nd, rscale, source, ns, charge, dip,
     1                    center, nterms, mpole_c)
      call bhdiff('bh2dformmpcd', mpole_f, mpole_c,
     1            nd, nterms, nfail)

c     ============================================================
c     bh2dmpevalp / bh2dmpevalg
c     Build a multipole then evaluate at far targets.
c     ============================================================
      call zerobh(mpole_f, nd, nterms)
      call bh2dformmpcd(nd, rscale, source, ns, charge, dip,
     1                  center, nterms, mpole_f)

c     re-place targets outside the source region for multipole eval
      do j = 1, nt
         targ(1, j) = center(1) + 5.0d0 + hkrand(0)
         targ(2, j) = center(2) + 5.0d0 + hkrand(0)
      enddo

      call zerocplx(vel_f, nd*nt)
      call zerocplx(vel_c, nd*nt)
      call bh2dmpevalp  (nd, rscale, center, mpole_f, nterms,
     1                   targ, nt, vel_f)
      call bh2dmpevalp_c(nd, rscale, center, mpole_f, nterms,
     1                   targ, nt, vel_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(vel_f(ii,j) - vel_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('bh2dmpevalp ', errmax, nfail)

      call zerocplx(vel_f, nd*nt)
      call zerocplx(vel_c, nd*nt)
      call zerocplx(grad_f, nd*3*nt)
      call zerocplx(grad_c, nd*3*nt)
      call bh2dmpevalg  (nd, rscale, center, mpole_f, nterms,
     1                   targ, nt, vel_f, grad_f)
      call bh2dmpevalg_c(nd, rscale, center, mpole_f, nterms,
     1                   targ, nt, vel_c, grad_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(vel_f(ii,j) - vel_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
         do kk = 1, 3
            do ii = 1, nd
               e = cdabs(grad_f(ii,kk,j) - grad_c(ii,kk,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
      enddo
      call report('bh2dmpevalg ', errmax, nfail)

c     restore targets near center for local form/eval
      do j = 1, nt
         targ(1, j) = center(1) + 0.1d0*hkrand(0)
         targ(2, j) = center(2) + 0.1d0*hkrand(0)
      enddo
c     move sources far away for local expansion form
      do i = 1, ns
         source(1, i) = center(1) + 5.0d0 + hkrand(0)
         source(2, i) = center(2) + 5.0d0 + hkrand(0)
      enddo

c     ============================================================
c     bh2dformtac
c     ============================================================
      call zerobh(local_f, nd, nterms)
      call zerobh(local_c, nd, nterms)
      call bh2dformtac  (nd, rscale, source, ns, charge, center,
     1                   nterms, local_f)
      call bh2dformtac_c(nd, rscale, source, ns, charge, center,
     1                   nterms, local_c)
      call bhdiff('bh2dformtac ', local_f, local_c,
     1            nd, nterms, nfail)

c     ============================================================
c     bh2dformtad
c     ============================================================
      call zerobh(local_f, nd, nterms)
      call zerobh(local_c, nd, nterms)
      call bh2dformtad  (nd, rscale, source, ns, dip, center,
     1                   nterms, local_f)
      call bh2dformtad_c(nd, rscale, source, ns, dip, center,
     1                   nterms, local_c)
      call bhdiff('bh2dformtad ', local_f, local_c,
     1            nd, nterms, nfail)

c     ============================================================
c     bh2dformtacd
c     ============================================================
      call zerobh(local_f, nd, nterms)
      call zerobh(local_c, nd, nterms)
      call bh2dformtacd  (nd, rscale, source, ns, charge, dip,
     1                    center, nterms, local_f)
      call bh2dformtacd_c(nd, rscale, source, ns, charge, dip,
     1                    center, nterms, local_c)
      call bhdiff('bh2dformtacd', local_f, local_c,
     1            nd, nterms, nfail)

c     ============================================================
c     bh2dtaevalp / bh2dtaevalg
c     ============================================================
      call zerobh(local_f, nd, nterms)
      call bh2dformtacd(nd, rscale, source, ns, charge, dip,
     1                  center, nterms, local_f)

      call zerocplx(vel_f, nd*nt)
      call zerocplx(vel_c, nd*nt)
      call bh2dtaevalp  (nd, rscale, center, local_f, nterms,
     1                   targ, nt, vel_f)
      call bh2dtaevalp_c(nd, rscale, center, local_f, nterms,
     1                   targ, nt, vel_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(vel_f(ii,j) - vel_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('bh2dtaevalp ', errmax, nfail)

      call zerocplx(vel_f, nd*nt)
      call zerocplx(vel_c, nd*nt)
      call zerocplx(grad_f, nd*3*nt)
      call zerocplx(grad_c, nd*3*nt)
      call bh2dtaevalg  (nd, rscale, center, local_f, nterms,
     1                   targ, nt, vel_f, grad_f)
      call bh2dtaevalg_c(nd, rscale, center, local_f, nterms,
     1                   targ, nt, vel_c, grad_c)
      errmax = 0.0d0
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(vel_f(ii,j) - vel_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
         do kk = 1, 3
            do ii = 1, nd
               e = cdabs(grad_f(ii,kk,j) - grad_c(ii,kk,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
      enddo
      call report('bh2dtaevalg ', errmax, nfail)

c     ============================================================
c     bh2dmpmp
c     ============================================================
      call zerobh(hexp2_f, nd, nterms2)
      call zerobh(hexp2_c, nd, nterms2)
      call bh2dmpmp  (nd, rscale1, center1, hexp1, nterms1,
     1                rscale2, center2, hexp2_f, nterms2, carray, ldc)
      call bh2dmpmp_c(nd, rscale1, center1, hexp1, nterms1,
     1                rscale2, center2, hexp2_c, nterms2, carray, ldc)
      call bhdiff('bh2dmpmp    ', hexp2_f, hexp2_c,
     1            nd, nterms2, nfail)

c     ============================================================
c     bh2dlocloc
c     ============================================================
      call zerobh(jexp2_f, nd, nterms2)
      call zerobh(jexp2_c, nd, nterms2)
      call bh2dlocloc  (nd, rscale1, center1, jexp1, nterms1,
     1                  rscale2, center2, jexp2_f, nterms2,
     2                  carray, ldc)
      call bh2dlocloc_c(nd, rscale1, center1, jexp1, nterms1,
     1                  rscale2, center2, jexp2_c, nterms2,
     2                  carray, ldc)
      call bhdiff('bh2dlocloc  ', jexp2_f, jexp2_c,
     1            nd, nterms2, nfail)

c     ============================================================
c     bh2dmploc
c     ============================================================
      call zerobh(jexp2_f, nd, nterms2)
      call zerobh(jexp2_c, nd, nterms2)
      call bh2dmploc  (nd, rscale1, center1, hexp1, nterms1,
     1                 rscale2, center2, jexp2_f, nterms2,
     2                 carray, ldc)
      call bh2dmploc_c(nd, rscale1, center1, hexp1, nterms1,
     1                 rscale2, center2, jexp2_c, nterms2,
     2                 carray, ldc)
      call bhdiff('bh2dmploc   ', jexp2_f, jexp2_c,
     1            nd, nterms2, nfail)

c     ============================================================
c     bh2dmpzero
c     ============================================================
      do n = 0, nterms
         do kk = 1, 5
            do ii = 1, nd
               mpzero_f(ii,kk,n) =
     1            dcmplx(2.0d0*hkrand(0)-1.0d0,
     2                   2.0d0*hkrand(0)-1.0d0)
               mpzero_c(ii,kk,n) = mpzero_f(ii,kk,n)
            enddo
         enddo
      enddo
      call bh2dmpzero  (nd, mpzero_f, nterms)
      call bh2dmpzero_c(nd, mpzero_c, nterms)
      errmax = 0.0d0
      do n = 0, nterms
         do kk = 1, 5
            do ii = 1, nd
               e = cdabs(mpzero_f(ii,kk,n) - mpzero_c(ii,kk,n))
               if (e .gt. errmax) errmax = e
               e = cdabs(mpzero_c(ii,kk,n))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
      enddo
      call report('bh2dmpzero  ', errmax, nfail)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all bhrouts2d cases match'

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
      subroutine zerobh(mp, nd, nterms)
      implicit none
      integer nd, nterms, i, k, ii
      complex *16 mp(nd, 5, 0:nterms)
      do i = 0, nterms
         do k = 1, 5
            do ii = 1, nd
               mp(ii, k, i) = dcmplx(0.0d0, 0.0d0)
            enddo
         enddo
      enddo
      return
      end

c     ----------------------------------------------------------------
      subroutine bhdiff(name, mp_f, mp_c, nd, nterms, nfail)
      implicit none
      character*(*) name
      integer nd, nterms, nfail
      complex *16 mp_f(nd, 5, 0:nterms), mp_c(nd, 5, 0:nterms)
      integer i, k, ii
      real *8 errmax, e
      external report
      errmax = 0.0d0
      do i = 0, nterms
         do k = 1, 5
            do ii = 1, nd
               e = cdabs(mp_f(ii,k,i) - mp_c(ii,k,i))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
      enddo
      call report(name, errmax, nfail)
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
