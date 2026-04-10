c     test_mbhkernels2d.f - differential test for the C translation
c     of src/modified-biharmonic/mbhkernels2d.f

      program test_mbhkernels2d
      implicit none

      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 20, ntmax = 15, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      real *8 charge(ndmax, nsmax)
      real *8 dipstr(ndmax, nsmax), dipvec(ndmax, 2, nsmax)
      real *8 quadstr(ndmax, nsmax), quadvec(ndmax, 3, nsmax)
      real *8 octstr(ndmax, nsmax), octvec(ndmax, 4, nsmax)
      real *8 beta, thresh, dummy

      real *8, allocatable :: pot_f(:,:), pot_c(:,:)
      real *8, allocatable :: grad_f(:,:,:), grad_c(:,:,:)
      real *8, allocatable :: hess_f(:,:,:), hess_c(:,:,:)

c     for mps tests
      integer ntermsmps
      parameter (ntermsmps = 3)
      complex *16 mbhmps(ndmax, 0:ntermsmps, nsmax)
      complex *16 ymps(ndmax, 0:ntermsmps, nsmax)

      integer ns, nt, nd, nds(2), ind, nfail
      integer i, j, ii, k1, k2
      real *8 errmax, e

      data nds / 1, 3 /

      ns = nsmax
      nt = ntmax
      nfail = 0
      thresh = 1.0d-15
      beta = 1.5d0

c     seed
      dummy = hkrand(1234)

c     fill sources/targets well separated
      do i = 1, ns
         sources(1, i) = hkrand(0) * 5.0d0
         sources(2, i) = hkrand(0) * 5.0d0
      enddo
      do j = 1, nt
         targ(1, j) = 10.0d0 + hkrand(0) * 5.0d0
         targ(2, j) = 10.0d0 + hkrand(0) * 5.0d0
      enddo

c     fill densities
      do i = 1, ns
         do ii = 1, ndmax
            charge(ii, i) = hkrand(0)
            dipstr(ii, i) = hkrand(0)
            dipvec(ii, 1, i) = hkrand(0)
            dipvec(ii, 2, i) = hkrand(0)
            quadstr(ii, i) = hkrand(0)
            quadvec(ii, 1, i) = hkrand(0)
            quadvec(ii, 2, i) = hkrand(0)
            quadvec(ii, 3, i) = hkrand(0)
            octstr(ii, i) = hkrand(0)
            octvec(ii, 1, i) = hkrand(0)
            octvec(ii, 2, i) = hkrand(0)
            octvec(ii, 3, i) = hkrand(0)
            octvec(ii, 4, i) = hkrand(0)
         enddo
      enddo

c     fill mps coefficients
      do i = 1, ns
         do j = 0, ntermsmps
            do ii = 1, ndmax
               mbhmps(ii, j, i) = dcmplx(hkrand(0),
     1            hkrand(0))
               ymps(ii, j, i) = dcmplx(hkrand(0),
     1            hkrand(0))
            enddo
         enddo
      enddo

c     ============================================================
c     Loop over nd values
c     ============================================================
      do ind = 1, 2
         nd = nds(ind)

         allocate(pot_f(nd, nt), pot_c(nd, nt))
         allocate(grad_f(nd, 2, nt), grad_c(nd, 2, nt))
         allocate(hess_f(nd, 3, nt), hess_c(nd, 3, nt))

c        ---- charge potential only ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call mbh2d_directcp_vec(nd, beta, sources,
     1      ns, charge, targ, nt, pot_f, thresh)
         call mbh2d_directcp_vec_c(nd, beta, sources,
     1      ns, charge, targ, nt, pot_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1      'directcp    ', nd, errmax, nfail)

c        ---- charge+dipole gradient ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grad_f, nd*2*nt)
         call zeroreal(grad_c, nd*2*nt)
         call mbh2d_directcdg_vec(nd, beta,
     1      sources, ns, charge, dipstr,
     2      dipvec, targ, nt, pot_f, grad_f,
     3      thresh)
         call mbh2d_directcdg_vec_c(nd, beta,
     1      sources, ns, charge, dipstr,
     2      dipvec, targ, nt, pot_c, grad_c,
     3      thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1      'cdg pot     ', nd, errmax, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1      'cdg grad    ', nd, errmax, nfail)

c        ---- charge+dipole+quad+oct hessian ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grad_f, nd*2*nt)
         call zeroreal(grad_c, nd*2*nt)
         call zeroreal(hess_f, nd*3*nt)
         call zeroreal(hess_c, nd*3*nt)
         call mbh2d_directcdqoh_vec(nd, beta,
     1      sources, ns, charge, dipstr,
     2      dipvec, quadstr, quadvec,
     3      octstr, octvec, targ, nt,
     4      pot_f, grad_f, hess_f, thresh)
         call mbh2d_directcdqoh_vec_c(nd, beta,
     1      sources, ns, charge, dipstr,
     2      dipvec, quadstr, quadvec,
     3      octstr, octvec, targ, nt,
     4      pot_c, grad_c, hess_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1      'cdqoh pot   ', nd, errmax, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1      'cdqoh grad  ', nd, errmax, nfail)
         call checkerr(hess_f, hess_c, nd*3*nt,
     1      'cdqoh hess  ', nd, errmax, nfail)

c        ---- mps potential ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call mbh2d_directmpsp_vec(nd, beta,
     1      sources, ns, mbhmps, ymps,
     2      ntermsmps, targ, nt, pot_f,
     3      thresh)
         call mbh2d_directmpsp_vec_c(nd, beta,
     1      sources, ns, mbhmps, ymps,
     2      ntermsmps, targ, nt, pot_c,
     3      thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1      'mpsp        ', nd, errmax, nfail)

c        ---- mps gradient ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grad_f, nd*2*nt)
         call zeroreal(grad_c, nd*2*nt)
         call mbh2d_directmpsg_vec(nd, beta,
     1      sources, ns, mbhmps, ymps,
     2      ntermsmps, targ, nt, pot_f,
     3      grad_f, thresh)
         call mbh2d_directmpsg_vec_c(nd, beta,
     1      sources, ns, mbhmps, ymps,
     2      ntermsmps, targ, nt, pot_c,
     3      grad_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1      'mpsg pot    ', nd, errmax, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1      'mpsg grad   ', nd, errmax, nfail)

c        ---- mps hessian ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grad_f, nd*2*nt)
         call zeroreal(grad_c, nd*2*nt)
         call zeroreal(hess_f, nd*3*nt)
         call zeroreal(hess_c, nd*3*nt)
         call mbh2d_directmpsh_vec(nd, beta,
     1      sources, ns, mbhmps, ymps,
     2      ntermsmps, targ, nt, pot_f,
     3      grad_f, hess_f, thresh)
         call mbh2d_directmpsh_vec_c(nd, beta,
     1      sources, ns, mbhmps, ymps,
     2      ntermsmps, targ, nt, pot_c,
     3      grad_c, hess_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1      'mpsh pot    ', nd, errmax, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1      'mpsh grad   ', nd, errmax, nfail)
         call checkerr(hess_f, hess_c, nd*3*nt,
     1      'mpsh hess   ', nd, errmax, nfail)

         deallocate(pot_f, pot_c)
         deallocate(grad_f, grad_c)
         deallocate(hess_f, hess_c)
      enddo

      if (nfail .gt. 0) then
         write(*,*) 'FAIL: ',nfail,
     1      ' mbhkernels2d tests failed'
         stop 1
      endif
      write(*,*) 'PASS: all mbhkernels2d cases match'

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
      subroutine checkerr(a_f, a_c, n, name, nd,
     1   errmax, nfail)
      implicit none
      integer n, nd, nfail, i
      real *8 a_f(n), a_c(n), errmax, e, tol
      character*(*) name
      tol = 1.0d-13
      errmax = 0.0d0
      do i = 1, n
         e = dabs(a_f(i) - a_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      if (errmax .lt. tol) then
         write(*,1000) name, nd, errmax
      else
         write(*,1001) name, nd, errmax
         nfail = nfail + 1
      endif
 1000 format(' [ ok ] ',a13,' nd=',i1,
     1   '  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a13,' nd=',i1,
     1   '  errmax=',1pe11.3)
      return
      end
