c     test_stfmm2d_e2e.f - end-to-end smoke test for the matlab path
c     (stfmm2d.m -> fmm2d.mex -> stfmm2d) through the C drop-in.
c
c     Calls stfmm2d directly with random sources/targets and a single
c     Stokeslet, then compares the FMM velocity at targets to a
c     brute-force direct sum.
c
c     The Stokeslet kernel (fmm2d's convention; G is 2*pi larger than
c     the standard one):
c       u_i(x) = sum_m G_{ij}(x, y^(m)) sigma^(m)_j
c       G_{ij}(x,y) = (-delta_{ij} log(r) + r_i r_j / r^2) / 2
c     with r = x - y, r = |r|.

      program test_stfmm2d_e2e
      implicit none
      external stfmm2d
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 1)

      real *8 source(2, ns), targ(2, nt)
      real *8 stoklet(nd, 2, ns)
      real *8 strslet(nd, 2, ns), strsvec(nd, 2, ns)

      real *8 pot(nd, 2, ns), pre(nd, ns)
      real *8 grad(nd, 2, 2, ns)
      real *8 pottarg(nd, 2, nt), pretarg(nd, nt)
      real *8 gradtarg(nd, 2, 2, nt)

      real *8 ref(2, nt)
      real *8 eps, dummy, errmax, ref_norm
      real *8 dx, dy, r2, lr, sx, sy, dot, inv_r2
      integer i, j
      integer ifstoklet, ifstrslet, ifppreg, ifppregtarg, ier
      integer nfail

      nfail = 0
      eps = 1.0d-6

      dummy = hkrand(1234)
      do i = 1, ns
         source(1, i) = hkrand(0)
         source(2, i) = hkrand(0)
         stoklet(1, 1, i) = 2.0d0 * hkrand(0) - 1.0d0
         stoklet(1, 2, i) = 2.0d0 * hkrand(0) - 1.0d0
         strslet(1, 1, i) = 0.0d0
         strslet(1, 2, i) = 0.0d0
         strsvec(1, 1, i) = 0.0d0
         strsvec(1, 2, i) = 0.0d0
      enddo
      do j = 1, nt
         targ(1, j) = 2.0d0 + hkrand(0)
         targ(2, j) = 2.0d0 + hkrand(0)
      enddo

      ifstoklet = 1
      ifstrslet = 0
      ifppreg = 0
      ifppregtarg = 1
      ier = 0

      call stfmm2d(nd, eps, ns, source, ifstoklet, stoklet,
     1     ifstrslet, strslet, strsvec, ifppreg, pot, pre, grad,
     2     nt, targ, ifppregtarg, pottarg, pretarg, gradtarg, ier)

c     Brute-force reference:
c       u_1 = sum_m 0.5 * [ -log(r)*sx + dx*(dx*sx + dy*sy)/r^2 ]
c       u_2 = sum_m 0.5 * [ -log(r)*sy + dy*(dx*sx + dy*sy)/r^2 ]
      do j = 1, nt
         ref(1, j) = 0.0d0
         ref(2, j) = 0.0d0
         do i = 1, ns
            dx = targ(1, j) - source(1, i)
            dy = targ(2, j) - source(2, i)
            r2 = dx * dx + dy * dy
            if (r2 .gt. 0.0d0) then
               sx = stoklet(1, 1, i)
               sy = stoklet(1, 2, i)
               lr = 0.5d0 * dlog(r2)
               inv_r2 = 1.0d0 / r2
               dot = dx * sx + dy * sy
               ref(1, j) = ref(1, j)
     1              + 0.5d0 * (-lr * sx + dx * dot * inv_r2)
               ref(2, j) = ref(2, j)
     1              + 0.5d0 * (-lr * sy + dy * dot * inv_r2)
            endif
         enddo
      enddo

      errmax = 0.0d0
      ref_norm = 0.0d0
      do j = 1, nt
         if (dabs(pottarg(1, 1, j) - ref(1, j)) .gt. errmax) then
            errmax = dabs(pottarg(1, 1, j) - ref(1, j))
         endif
         if (dabs(pottarg(1, 2, j) - ref(2, j)) .gt. errmax) then
            errmax = dabs(pottarg(1, 2, j) - ref(2, j))
         endif
         if (dabs(ref(1, j)) .gt. ref_norm) ref_norm = dabs(ref(1, j))
         if (dabs(ref(2, j)) .gt. ref_norm) ref_norm = dabs(ref(2, j))
      enddo

      write(*,1000) eps, errmax, errmax / ref_norm
 1000 format(' stfmm2d Stokeslet: eps=',1pe10.3,
     1   ' max abs err=',1pe10.3, ' rel err=',1pe10.3)

      if (errmax / ref_norm .gt. eps * 100.0d0) then
         write(*,*) ' FAIL: relative error exceeds 100*eps'
         nfail = nfail + 1
      endif

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: stfmm2d end-to-end smoke test'

      end
