c     test_cfmm2d_e2e.f - end-to-end smoke test for the matlab path
c     (cfmm2d.m -> fmm2d.mex -> cfmm2d_ndiv) through the C drop-in.
c
c     Calls cfmm2d_ndiv directly with random sources/targets and complex
c     charges and verifies the FMM output against a brute-force direct
c     sum.
c
c     The cfmm2d charge kernel:
c       Re(pot(z))  = Re(sum_i charge_i * log(|z - z_i|))
c       grad(z)     = sum_i charge_i * 1/(z - z_i)        (d/dz)
c
c     The IMAGINARY part of pot depends on the harmonic-conjugate gauge
c     (cfmm2d uses a Cauchy expansion which integrates 1/z to the
c     complex log(z), not log(|z|)). The upstream test_cfmm2d.f uses
c     REAL charges and compares only Re(pot); we do the same for the
c     potential test. The gradient is unambiguous (1/(z-z_i)) and we
c     test it with COMPLEX charges to exercise the full complex
c     output code path.

      program test_cfmm2d_e2e
      implicit none
      external cfmm2d_ndiv
      external hndiv2d
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 1)

      real *8 sources(2, ns), targ(2, nt)
      complex *16 charges(nd, ns), dipstr(nd, ns)
      complex *16 pot(nd, ns), grad(nd, ns), hess(nd, ns)
      complex *16 pottarg(nd, nt), gradtarg(nd, nt), hesstarg(nd, nt)

      complex *16 pot_ref(nd, nt), grad_ref(nd, nt)
      real *8 timeinfo(8)

      real *8 eps, dummy, errp, errg, refp, refg, dx, dy, r2, lr
      complex *16 contrib, dz, gcontrib
      integer i, j, ifcharge, ifdipole, iper, ifpgh, ifpghtarg
      integer ier, ndiv, idivflag, ifnear
      integer nfail

      complex *16 eye
      data eye /(0.0d0, 1.0d0)/

      nfail = 0
      eps = 1.0d-6

c     Test 1 input data: REAL charges (for the Re(pot) check)
      dummy = hkrand(1234)
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
         charges(1, i) = (2.0d0 * hkrand(0) - 1.0d0)
         dipstr(1, i) = (0.0d0, 0.0d0)
      enddo
      do j = 1, nt
         targ(1, j) = hkrand(0)
         targ(2, j) = hkrand(0)
      enddo

c     Charges-only, pot+grad at targets.
      ifcharge = 1
      ifdipole = 0
      ifpgh = 0
      ifpghtarg = 2
      iper = 0
      ifnear = 1
      ier = 0

      call hndiv2d(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg,
     1     ndiv, idivflag)

      do j = 1, nt
         pottarg(1, j) = (0.0d0, 0.0d0)
         gradtarg(1, j) = (0.0d0, 0.0d0)
         hesstarg(1, j) = (0.0d0, 0.0d0)
      enddo
      do j = 1, 8
         timeinfo(j) = 0.0d0
      enddo

      call cfmm2d_ndiv(nd, eps, ns, sources, ifcharge, charges,
     1     ifdipole, dipstr, iper, ifpgh, pot, grad, hess,
     2     nt, targ, ifpghtarg, pottarg, gradtarg, hesstarg,
     3     ndiv, idivflag, ifnear, timeinfo, ier)

c     Brute-force reference (test 1, real charges):
c       pot_ref(j)  = sum_i charge_i * log(|z_j - z_i|)
      do j = 1, nt
         pot_ref(1, j) = (0.0d0, 0.0d0)
         do i = 1, ns
            dx = targ(1, j) - sources(1, i)
            dy = targ(2, j) - sources(2, i)
            r2 = dx * dx + dy * dy
            if (r2 .gt. 0.0d0) then
               lr = 0.5d0 * dlog(r2)
               contrib = lr * charges(1, i)
               pot_ref(1, j) = pot_ref(1, j) + contrib
            endif
         enddo
      enddo

c     Compare FMM Re(pot) vs direct sum Re(pot). Charges are real here
c     so the imaginary part of pot is the harmonic conjugate (which we
c     do not check — its gauge is set by the multipole expansion).
      errp = 0.0d0
      refp = 0.0d0
      do j = 1, nt
         if (dabs(dble(pottarg(1, j)) - dble(pot_ref(1, j)))
     1       .gt. errp) then
            errp = dabs(dble(pottarg(1, j)) - dble(pot_ref(1, j)))
         endif
         if (dabs(dble(pot_ref(1, j))) .gt. refp) then
            refp = dabs(dble(pot_ref(1, j)))
         endif
      enddo

c     Test 2 input data: COMPLEX charges (for the gradient check)
      do i = 1, ns
         charges(1, i) = (2.0d0 * hkrand(0) - 1.0d0)
     1        + eye * (2.0d0 * hkrand(0) - 1.0d0)
      enddo
      do j = 1, nt
         pottarg(1, j) = (0.0d0, 0.0d0)
         gradtarg(1, j) = (0.0d0, 0.0d0)
         hesstarg(1, j) = (0.0d0, 0.0d0)
      enddo

      call cfmm2d_ndiv(nd, eps, ns, sources, ifcharge, charges,
     1     ifdipole, dipstr, iper, ifpgh, pot, grad, hess,
     2     nt, targ, ifpghtarg, pottarg, gradtarg, hesstarg,
     3     ndiv, idivflag, ifnear, timeinfo, ier)

c     Brute-force reference (test 2, complex charges):
c       grad_ref(j) = sum_i charge_i / (z_j - z_i)
      do j = 1, nt
         grad_ref(1, j) = (0.0d0, 0.0d0)
         do i = 1, ns
            dx = targ(1, j) - sources(1, i)
            dy = targ(2, j) - sources(2, i)
            r2 = dx * dx + dy * dy
            if (r2 .gt. 0.0d0) then
               dz = dx + eye * dy
               gcontrib = charges(1, i) / dz
               grad_ref(1, j) = grad_ref(1, j) + gcontrib
            endif
         enddo
      enddo

c     Compare FMM grad vs direct sum (full complex)
      errg = 0.0d0
      refg = 0.0d0
      do j = 1, nt
         if (cdabs(gradtarg(1, j) - grad_ref(1, j)) .gt. errg) then
            errg = cdabs(gradtarg(1, j) - grad_ref(1, j))
         endif
         if (cdabs(grad_ref(1, j)) .gt. refg) then
            refg = cdabs(grad_ref(1, j))
         endif
      enddo

      write(*,1000) eps, errp, errp / refp
      write(*,1001) eps, errg, errg / refg
 1000 format(' cfmm2d charge Re(pot): eps=',1pe10.3,
     1   ' max abs err=',1pe10.3, ' rel err=',1pe10.3)
 1001 format(' cfmm2d charge grad   : eps=',1pe10.3,
     1   ' max abs err=',1pe10.3, ' rel err=',1pe10.3)

      if (errp / refp .gt. eps * 100.0d0) then
         write(*,*) ' FAIL: pot relative error exceeds 100*eps'
         nfail = nfail + 1
      endif
      if (errg / refg .gt. eps * 100.0d0) then
         write(*,*) ' FAIL: grad relative error exceeds 100*eps'
         nfail = nfail + 1
      endif

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: cfmm2d_ndiv end-to-end smoke test'

      end
