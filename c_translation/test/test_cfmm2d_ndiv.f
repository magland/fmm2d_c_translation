c     test_cfmm2d_ndiv.f - differential test for the C translation of
c     src/laplace/cfmm2d_ndiv.f.
c
c     Compares the Fortran reference cfmm2d_ndiv to the C port
c     cfmm2d_ndiv_c over a sweep of the 9 combinations of
c     (ifcharge, ifdipole, ifpgh, ifpghtarg). Both sides use the same
c     random inputs; the outputs are compared element-by-element to
c     tolerance 1e-15 (the -O0 build makes the computation
c     bit-deterministic so any nonzero error is a bug).
c
c     timeinfo is not compared (wall-clock noise). ier is compared
c     (should be zero on both sides).
c
c     ndiv = 8, idivflag = 0, ifnear = 1, iper = 1.

      program test_cfmm2d_ndiv
      implicit none

      external cfmm2d_ndiv, cfmm2d_ndiv_c
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 2)
      real *8 sources(2, ns), targ(2, nt)
      complex *16 charge(nd, ns), dipstr(nd, ns)
      complex *16 pot_f(nd, ns),  grad_f(nd, ns),  hess_f(nd, ns)
      complex *16 pot_c(nd, ns),  grad_c(nd, ns),  hess_c(nd, ns)
      complex *16 pottarg_f(nd, nt), gradtarg_f(nd, nt)
      complex *16 hesstarg_f(nd, nt)
      complex *16 pottarg_c(nd, nt), gradtarg_c(nd, nt)
      complex *16 hesstarg_c(nd, nt)
      real *8 timeinfo_f(8), timeinfo_c(8)

      real *8 eps, dummy, errmax, e
      integer i, j, ic
      integer ifcharge, ifdipole, iper_f, iper_c, ifpgh, ifpghtarg
      integer ier_f, ier_c, nfail
      integer ndiv, idivflag, ifnear

c     Test configurations: 9 combos of (ifcharge, ifdipole, ifpgh,
c     ifpghtarg). Each must pass independently.
      integer ncases
      parameter (ncases = 9)
      integer cases(4, ncases)

      data cases /
     1  1, 0, 1, 1,
     1  1, 0, 2, 2,
     1  1, 0, 3, 3,
     1  0, 1, 1, 1,
     1  0, 1, 2, 2,
     1  0, 1, 3, 3,
     1  1, 1, 1, 1,
     1  1, 1, 2, 2,
     1  1, 1, 3, 3 /

      nfail = 0
      eps = 1.0d-6

      ndiv = 8
      idivflag = 0
      ifnear = 1

c     Initialize RNG and input data
      dummy = hkrand(1234)
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
         do j = 1, nd
            charge(j, i) = dcmplx(hkrand(0) - 0.5d0,
     1                            hkrand(0) - 0.5d0)
            dipstr(j, i) = dcmplx(hkrand(0) - 0.5d0,
     1                            hkrand(0) - 0.5d0)
         enddo
      enddo
      do i = 1, nt
         targ(1, i) = hkrand(0)
         targ(2, i) = hkrand(0)
      enddo

c     ================================================================
c     Sweep test cases for cfmm2d_ndiv.
c     ================================================================
      do ic = 1, ncases
         ifcharge  = cases(1, ic)
         ifdipole  = cases(2, ic)
         ifpgh     = cases(3, ic)
         ifpghtarg = cases(4, ic)

c        zero output buffers
         do i = 1, ns
            do j = 1, nd
               pot_f(j, i)  = (0.0d0, 0.0d0)
               grad_f(j, i) = (0.0d0, 0.0d0)
               hess_f(j, i) = (0.0d0, 0.0d0)
               pot_c(j, i)  = (0.0d0, 0.0d0)
               grad_c(j, i) = (0.0d0, 0.0d0)
               hess_c(j, i) = (0.0d0, 0.0d0)
            enddo
         enddo
         do i = 1, nt
            do j = 1, nd
               pottarg_f(j, i)  = (0.0d0, 0.0d0)
               gradtarg_f(j, i) = (0.0d0, 0.0d0)
               hesstarg_f(j, i) = (0.0d0, 0.0d0)
               pottarg_c(j, i)  = (0.0d0, 0.0d0)
               gradtarg_c(j, i) = (0.0d0, 0.0d0)
               hesstarg_c(j, i) = (0.0d0, 0.0d0)
            enddo
         enddo

         do i = 1, 8
            timeinfo_f(i) = 0.0d0
            timeinfo_c(i) = 0.0d0
         enddo

         iper_f = 1
         iper_c = 1
         ier_f = 0
         ier_c = 0

         call cfmm2d_ndiv(nd, eps, ns, sources, ifcharge, charge,
     1        ifdipole, dipstr, iper_f, ifpgh, pot_f, grad_f,
     2        hess_f, nt, targ, ifpghtarg, pottarg_f, gradtarg_f,
     3        hesstarg_f, ndiv, idivflag, ifnear, timeinfo_f, ier_f)

         call cfmm2d_ndiv_c(nd, eps, ns, sources, ifcharge, charge,
     1        ifdipole, dipstr, iper_c, ifpgh, pot_c, grad_c,
     2        hess_c, nt, targ, ifpghtarg, pottarg_c, gradtarg_c,
     3        hesstarg_c, ndiv, idivflag, ifnear, timeinfo_c, ier_c)

         errmax = 0.0d0
         if (ifpgh .ge. 1) then
            do i = 1, ns
               do j = 1, nd
                  e = abs(pot_f(j, i) - pot_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif
         if (ifpgh .ge. 2) then
            do i = 1, ns
               do j = 1, nd
                  e = abs(grad_f(j, i) - grad_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif
         if (ifpgh .ge. 3) then
            do i = 1, ns
               do j = 1, nd
                  e = abs(hess_f(j, i) - hess_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif
         if (ifpghtarg .ge. 1) then
            do i = 1, nt
               do j = 1, nd
                  e = abs(pottarg_f(j, i) - pottarg_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif
         if (ifpghtarg .ge. 2) then
            do i = 1, nt
               do j = 1, nd
                  e = abs(gradtarg_f(j, i) - gradtarg_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif
         if (ifpghtarg .ge. 3) then
            do i = 1, nt
               do j = 1, nd
                  e = abs(hesstarg_f(j, i) - hesstarg_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif

         if (ier_f .ne. ier_c) then
            write(*,1002) ifcharge, ifdipole, ifpgh, ifpghtarg,
     1           ier_f, ier_c
            nfail = nfail + 1
         else if (errmax .gt. 1.0d-15) then
            write(*,1001) ifcharge, ifdipole, ifpgh, ifpghtarg, errmax
            nfail = nfail + 1
         else
            write(*,1000) ifcharge, ifdipole, ifpgh, ifpghtarg, errmax
         endif
      enddo

 1000 format(' [ ok ] cfmm2d_ndiv ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  errmax=',1pe10.3)
 1001 format(' [FAIL] cfmm2d_ndiv ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  errmax=',1pe10.3)
 1002 format(' [FAIL] cfmm2d_ndiv ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  ier mismatch ',i6,' vs ',i6)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all cfmm2d_ndiv cases match'

      end
