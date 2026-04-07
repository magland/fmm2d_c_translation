c     test_rfmm2d_e2e.f - end-to-end smoke test for the matlab path
c     (rfmm2d.m -> fmm2d.mex -> rfmm2d_ndiv).
c
c     This test calls rfmm2d_ndiv directly with random sources/targets,
c     compares the FMM output to a brute-force direct summation, and
c     verifies that the relative error is under the requested eps. It
c     does NOT do a Fortran/C diff — it tests the absolute correctness
c     of the C-drop-in library against the underlying physics.

      program test_rfmm2d_e2e
      implicit none
      external rfmm2d_ndiv
      external hndiv2d
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 1)

      real *8 sources(2, ns), targ(2, nt)
      real *8 charges(nd, ns), dipstr(nd, ns)
      real *8 dipvec(nd, 2, ns)
      real *8 pot(nd, ns), grad(nd, 2, ns), hess(nd, 3, ns)
      real *8 pottarg(nd, nt), gradtarg(nd, 2, nt)
      real *8 hesstarg(nd, 3, nt)

      real *8 pot_ref(nd, nt)
      real *8 timeinfo(8)

      real *8 eps, dummy, errmax, ref_norm, dx, dy, r2, contrib
      integer i, j, k, ifcharge, ifdipole, iper, ifpgh, ifpghtarg
      integer ier, ndiv, idivflag, ifnear
      integer nfail

      nfail = 0
      eps = 1.0d-6

      dummy = hkrand(1234)
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
         charges(1, i) = 2.0d0 * hkrand(0) - 1.0d0
         dipstr(1, i)  = 2.0d0 * hkrand(0) - 1.0d0
         dipvec(1, 1, i) = 2.0d0 * hkrand(0) - 1.0d0
         dipvec(1, 2, i) = 2.0d0 * hkrand(0) - 1.0d0
      enddo
      do j = 1, nt
         targ(1, j) = hkrand(0)
         targ(2, j) = hkrand(0)
      enddo

c     Charges-only, potential-only test at targets.
      ifcharge = 1
      ifdipole = 0
      ifpgh = 0
      ifpghtarg = 1
      iper = 0
      ifnear = 1
      ier = 0

      call hndiv2d(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg,
     1     ndiv, idivflag)

      do j = 1, nt
         pottarg(1, j) = 0.0d0
      enddo
      do j = 1, 8
         timeinfo(j) = 0.0d0
      enddo

      call rfmm2d_ndiv(nd, eps, ns, sources, ifcharge, charges,
     1     ifdipole, dipstr, dipvec, iper, ifpgh, pot, grad, hess,
     2     nt, targ, ifpghtarg, pottarg, gradtarg, hesstarg,
     3     ndiv, idivflag, ifnear, timeinfo, ier)

c     Brute-force reference: pot(j) = sum_i charge_i * log(|targ_j - src_i|)
      do j = 1, nt
         pot_ref(1, j) = 0.0d0
         do i = 1, ns
            dx = targ(1, j) - sources(1, i)
            dy = targ(2, j) - sources(2, i)
            r2 = dx * dx + dy * dy
            if (r2 .gt. 0.0d0) then
               contrib = 0.5d0 * dlog(r2) * charges(1, i)
               pot_ref(1, j) = pot_ref(1, j) + contrib
            endif
         enddo
      enddo

c     Compare FMM result vs direct sum
      errmax = 0.0d0
      ref_norm = 0.0d0
      do j = 1, nt
         if (dabs(pottarg(1, j) - pot_ref(1, j)) .gt. errmax) then
            errmax = dabs(pottarg(1, j) - pot_ref(1, j))
         endif
         if (dabs(pot_ref(1, j)) .gt. ref_norm) then
            ref_norm = dabs(pot_ref(1, j))
         endif
      enddo

      write(*,1000) eps, errmax, errmax / ref_norm
 1000 format(' charge potential at targets: eps=',1pe10.3,
     1   ' max abs err=',1pe10.3, ' rel err=',1pe10.3)

      if (errmax / ref_norm .gt. eps * 100.0d0) then
         write(*,*) ' FAIL: relative error exceeds 100*eps'
         nfail = nfail + 1
      endif

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: rfmm2d_ndiv end-to-end smoke test'

      end
