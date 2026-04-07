c     test_rfmm2d_ndiv.f - differential test for the C translation of
c     src/laplace/rfmm2d_ndiv.f.
c
c     Compares the Fortran reference rfmm2d_ndiv to the C port
c     rfmm2d_ndiv_c over a sweep of the 9 combinations of
c     (ifcharge, ifdipole, ifpgh, ifpghtarg). Both sides use the same
c     random inputs (real-valued charges, dipole strengths and dipole
c     orientation vectors); the outputs are compared element-by-element
c     to tolerance 1e-15 (the -O0 build makes the computation
c     bit-deterministic so any nonzero error is a bug).
c
c     timeinfo is not compared (wall-clock noise). ier is compared
c     (should be zero on both sides).
c
c     ndiv = 8, idivflag = 0, ifnear = 1, iper = 1.

      program test_rfmm2d_ndiv
      implicit none

      external rfmm2d_ndiv, rfmm2d_ndiv_c
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 2)
      real *8 sources(2, ns), targ(2, nt)
      real *8 charge(nd, ns), dipstr(nd, ns)
      real *8 dipvec(nd, 2, ns)
      real *8 pot_f(nd, ns),  pot_c(nd, ns)
      real *8 grad_f(nd, 2, ns), grad_c(nd, 2, ns)
      real *8 hess_f(nd, 3, ns), hess_c(nd, 3, ns)
      real *8 pottarg_f(nd, nt), pottarg_c(nd, nt)
      real *8 gradtarg_f(nd, 2, nt), gradtarg_c(nd, 2, nt)
      real *8 hesstarg_f(nd, 3, nt), hesstarg_c(nd, 3, nt)
      real *8 timeinfo_f(8), timeinfo_c(8)

      real *8 eps, dummy, errmax, e
      integer i, j, k, ic
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
            charge(j, i) = hkrand(0) - 0.5d0
            dipstr(j, i) = hkrand(0) - 0.5d0
            dipvec(j, 1, i) = 2.0d0 * hkrand(0) - 1.0d0
            dipvec(j, 2, i) = 2.0d0 * hkrand(0) - 1.0d0
         enddo
      enddo
      do i = 1, nt
         targ(1, i) = hkrand(0)
         targ(2, i) = hkrand(0)
      enddo

c     ================================================================
c     Sweep test cases for rfmm2d_ndiv.
c     ================================================================
      do ic = 1, ncases
         ifcharge  = cases(1, ic)
         ifdipole  = cases(2, ic)
         ifpgh     = cases(3, ic)
         ifpghtarg = cases(4, ic)

c        zero output buffers
         do i = 1, ns
            do j = 1, nd
               pot_f(j, i) = 0.0d0
               pot_c(j, i) = 0.0d0
            enddo
            do k = 1, 2
               do j = 1, nd
                  grad_f(j, k, i) = 0.0d0
                  grad_c(j, k, i) = 0.0d0
               enddo
            enddo
            do k = 1, 3
               do j = 1, nd
                  hess_f(j, k, i) = 0.0d0
                  hess_c(j, k, i) = 0.0d0
               enddo
            enddo
         enddo
         do i = 1, nt
            do j = 1, nd
               pottarg_f(j, i) = 0.0d0
               pottarg_c(j, i) = 0.0d0
            enddo
            do k = 1, 2
               do j = 1, nd
                  gradtarg_f(j, k, i) = 0.0d0
                  gradtarg_c(j, k, i) = 0.0d0
               enddo
            enddo
            do k = 1, 3
               do j = 1, nd
                  hesstarg_f(j, k, i) = 0.0d0
                  hesstarg_c(j, k, i) = 0.0d0
               enddo
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

         call rfmm2d_ndiv(nd, eps, ns, sources, ifcharge, charge,
     1        ifdipole, dipstr, dipvec, iper_f, ifpgh, pot_f,
     2        grad_f, hess_f, nt, targ, ifpghtarg, pottarg_f,
     3        gradtarg_f, hesstarg_f, ndiv, idivflag, ifnear,
     4        timeinfo_f, ier_f)

         call rfmm2d_ndiv_c(nd, eps, ns, sources, ifcharge, charge,
     1        ifdipole, dipstr, dipvec, iper_c, ifpgh, pot_c,
     2        grad_c, hess_c, nt, targ, ifpghtarg, pottarg_c,
     3        gradtarg_c, hesstarg_c, ndiv, idivflag, ifnear,
     4        timeinfo_c, ier_c)

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
               do k = 1, 2
                  do j = 1, nd
                     e = abs(grad_f(j, k, i) - grad_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
               enddo
            enddo
         endif
         if (ifpgh .ge. 3) then
            do i = 1, ns
               do k = 1, 3
                  do j = 1, nd
                     e = abs(hess_f(j, k, i) - hess_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
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
               do k = 1, 2
                  do j = 1, nd
                     e = abs(gradtarg_f(j, k, i)
     1                    - gradtarg_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
               enddo
            enddo
         endif
         if (ifpghtarg .ge. 3) then
            do i = 1, nt
               do k = 1, 3
                  do j = 1, nd
                     e = abs(hesstarg_f(j, k, i)
     1                    - hesstarg_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
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

 1000 format(' [ ok ] rfmm2d_ndiv ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  errmax=',1pe10.3)
 1001 format(' [FAIL] rfmm2d_ndiv ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  errmax=',1pe10.3)
 1002 format(' [FAIL] rfmm2d_ndiv ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  ier mismatch ',i6,' vs ',i6)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all rfmm2d_ndiv cases match'

      end
