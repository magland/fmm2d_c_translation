c     test_mbhrouts2d.f - differential test for the C translation
c     of src/modified-biharmonic/mbhrouts2d.f

      program test_mbhrouts2d
      implicit none

      real *8 hkrand
      external hkrand

      integer nfail
      real *8 beta, rscale, r, dummy
      real *8 center(2), targ(2, 10)
      integer nterms, ntarg, nd, i, j, k

c     outputs for rk/rksc/rmk
      real *8 pow_f(0:20), pow_c(0:20)
      real *8 dpow_f(0:20), dpow_c(0:20)

c     outputs for init_carray
      integer ldc
      parameter (ldc = 10)
      real *8 carray_f(0:ldc, 0:ldc)
      real *8 carray_c(0:ldc, 0:ldc)

c     outputs for convtomp_vec
      integer ns
      parameter (ns = 5)
      real *8 charge(3, ns), dipstr(3, ns)
      real *8 dipvec(3, 2, ns)
      real *8 quadstr(3, ns), quadvec(3, 3, ns)
      real *8 octstr(3, ns), octvec(3, 4, ns)
      complex *16 mbhmp_f(3, 0:3, ns), mbhmp_c(3, 0:3, ns)
      complex *16 ymp_f(3, 0:3, ns), ymp_c(3, 0:3, ns)
      integer ifcharge, ifdipole, ifquad, ifoct

c     for mpeval/taeval
      complex *16 mbhmpole(3, 0:5), ympole(3, 0:5)
      real *8 pot_f(3, 10), pot_c(3, 10)
      real *8 grad_f(3, 2, 10), grad_c(3, 2, 10)

      real *8 errmax, e, tol
      tol = 1.0d-13

      nfail = 0
      dummy = hkrand(1234)
      beta = 1.5d0
      rscale = 0.8d0
      nterms = 8
      r = 2.5d0

c     ---- mbh2d_rk ----
      call mbh2d_rk(pow_f, dpow_f, r, beta,
     1   rscale, nterms)
      call mbh2d_rk_c(pow_c, dpow_c, r, beta,
     1   rscale, nterms)
      errmax = 0.0d0
      do i = 0, nterms
         e = dabs(pow_f(i) - pow_c(i))
         if (e .gt. errmax) errmax = e
         e = dabs(dpow_f(i) - dpow_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('rk           ', errmax, nfail)

c     ---- mbh2d_rksc ----
      call mbh2d_rksc(pow_f, dpow_f, r, beta,
     1   rscale, nterms)
      call mbh2d_rksc_c(pow_c, dpow_c, r, beta,
     1   rscale, nterms)
      errmax = 0.0d0
      do i = 0, nterms
         e = dabs(pow_f(i) - pow_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('rksc         ', errmax, nfail)

c     ---- mbh2d_rmk ----
      call mbh2d_rmk(pow_f, dpow_f, r, beta,
     1   rscale, nterms)
      call mbh2d_rmk_c(pow_c, dpow_c, r, beta,
     1   rscale, nterms)
      errmax = 0.0d0
      do i = 0, nterms
         e = dabs(pow_f(i) - pow_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('rmk          ', errmax, nfail)

c     ---- mbh2d_init_carray ----
      do j = 0, ldc
         do i = 0, ldc
            carray_f(i, j) = 0.0d0
            carray_c(i, j) = 0.0d0
         enddo
      enddo
      call mbh2d_init_carray(carray_f, ldc)
      call mbh2d_init_carray_c(carray_c, ldc)
      errmax = 0.0d0
      do j = 0, ldc
         do i = 0, ldc
            e = dabs(carray_f(i,j)-carray_c(i,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      call report('init_carray  ', errmax, nfail)

c     ---- convtomp_vec ----
      nd = 3
      do i = 1, ns
         do j = 1, nd
            charge(j, i) = hkrand(0)
            dipstr(j, i) = hkrand(0)
            dipvec(j, 1, i) = hkrand(0)
            dipvec(j, 2, i) = hkrand(0)
            quadstr(j, i) = hkrand(0)
            quadvec(j, 1, i) = hkrand(0)
            quadvec(j, 2, i) = hkrand(0)
            quadvec(j, 3, i) = hkrand(0)
            octstr(j, i) = hkrand(0)
            octvec(j, 1, i) = hkrand(0)
            octvec(j, 2, i) = hkrand(0)
            octvec(j, 3, i) = hkrand(0)
            octvec(j, 4, i) = hkrand(0)
         enddo
      enddo
      nterms = 3
      ifcharge = 1
      ifdipole = 1
      ifquad = 1
      ifoct = 1
      do i = 1, ns
         do j = 0, nterms
            do k = 1, nd
               mbhmp_f(k, j, i) = 0
               mbhmp_c(k, j, i) = 0
               ymp_f(k, j, i) = 0
               ymp_c(k, j, i) = 0
            enddo
         enddo
      enddo
      call mbh2dconvtomp_vec(nd, beta, ns,
     1   ifcharge, charge, ifdipole, dipstr,
     2   dipvec, ifquad, quadstr, quadvec,
     3   ifoct, octstr, octvec, nterms,
     4   mbhmp_f, ymp_f)
      call mbh2dconvtomp_vec_c(nd, beta, ns,
     1   ifcharge, charge, ifdipole, dipstr,
     2   dipvec, ifquad, quadstr, quadvec,
     3   ifoct, octstr, octvec, nterms,
     4   mbhmp_c, ymp_c)
      errmax = 0.0d0
      do i = 1, ns
         do j = 0, nterms
            do k = 1, nd
               e = cdabs(mbhmp_f(k,j,i)-
     1            mbhmp_c(k,j,i))
               if (e .gt. errmax) errmax = e
               e = cdabs(ymp_f(k,j,i)-
     1            ymp_c(k,j,i))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
      enddo
      call report('convtomp     ', errmax, nfail)

c     ---- mpeval_p ----
      nterms = 3
      nd = 2
      rscale = 1.0d0
      center(1) = 0.0d0
      center(2) = 0.0d0
      ntarg = 5
      do j = 0, nterms
         do k = 1, nd
            mbhmpole(k, j) = dcmplx(hkrand(0),
     1         hkrand(0))
            ympole(k, j) = dcmplx(hkrand(0),
     1         hkrand(0))
         enddo
      enddo
      do i = 1, ntarg
         targ(1, i) = 2.0d0 + hkrand(0)
         targ(2, i) = 2.0d0 + hkrand(0)
      enddo
      call zeroreal(pot_f, nd*ntarg)
      call zeroreal(pot_c, nd*ntarg)
      call mbh2dmpevalp_vec(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_f)
      call mbh2dmpevalp_vec_c(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_c)
      errmax = 0.0d0
      do i = 1, nd*ntarg
         e = dabs(pot_f(i,1)-pot_c(i,1))
         if (e .gt. errmax) errmax = e
      enddo
      call report('mpevalp      ', errmax, nfail)

c     ---- mpeval_g ----
      call zeroreal(pot_f, nd*ntarg)
      call zeroreal(pot_c, nd*ntarg)
      call zeroreal(grad_f, nd*2*ntarg)
      call zeroreal(grad_c, nd*2*ntarg)
      call mbh2dmpevalg_vec(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_f, grad_f)
      call mbh2dmpevalg_vec_c(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_c, grad_c)
      errmax = 0.0d0
      do i = 1, nd*ntarg
         e = dabs(pot_f(i,1)-pot_c(i,1))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1, nd*2*ntarg
         e = dabs(grad_f(i,1,1)-grad_c(i,1,1))
         if (e .gt. errmax) errmax = e
      enddo
      call report('mpevalg      ', errmax, nfail)

c     ---- taeval_p ----
      call zeroreal(pot_f, nd*ntarg)
      call zeroreal(pot_c, nd*ntarg)
      call mbh2dtaevalp_vec(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_f)
      call mbh2dtaevalp_vec_c(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_c)
      errmax = 0.0d0
      do i = 1, nd*ntarg
         e = dabs(pot_f(i,1)-pot_c(i,1))
         if (e .gt. errmax) errmax = e
      enddo
      call report('taevalp      ', errmax, nfail)

c     ---- taeval_g ----
      call zeroreal(pot_f, nd*ntarg)
      call zeroreal(pot_c, nd*ntarg)
      call zeroreal(grad_f, nd*2*ntarg)
      call zeroreal(grad_c, nd*2*ntarg)
      call mbh2dtaevalg_vec(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_f, grad_f)
      call mbh2dtaevalg_vec_c(nd, beta, rscale,
     1   center, mbhmpole, ympole, nterms,
     2   targ, ntarg, pot_c, grad_c)
      errmax = 0.0d0
      do i = 1, nd*ntarg
         e = dabs(pot_f(i,1)-pot_c(i,1))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1, nd*2*ntarg
         e = dabs(grad_f(i,1,1)-grad_c(i,1,1))
         if (e .gt. errmax) errmax = e
      enddo
      call report('taevalg      ', errmax, nfail)

      if (nfail .gt. 0) then
         write(*,*) 'FAIL: ',nfail,
     1      ' mbhrouts2d tests failed'
         stop 1
      endif
      write(*,*) 'PASS: all mbhrouts2d cases match'

      end


c     ----------------------------------------------------------------
      subroutine zeroreal(a, n)
      implicit none
      integer n, i
      real *8 a(n)
      do i = 1, n
         a(i) = 0.0d0
      enddo
      return
      end


c     ----------------------------------------------------------------
      subroutine report(name, errmax, nfail)
      implicit none
      character*(*) name
      integer nfail
      real *8 errmax, tol
      tol = 1.0d-13
      if (errmax .lt. tol) then
         write(*,1000) name, errmax
      else
         write(*,1001) name, errmax
         nfail = nfail + 1
      endif
 1000 format(' [ ok ] ',a13,
     1   '  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a13,
     1   '  errmax=',1pe11.3)
      return
      end
