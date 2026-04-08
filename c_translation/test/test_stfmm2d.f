c     test_stfmm2d.f - differential test for the C translation of
c     src/stokes/stfmm2d.f.
c
c     Compares the Fortran reference stfmm2d to the C port stfmm2d_c
c     over a sweep of (ifstoklet, ifstrslet, ifppreg, ifppregtarg)
c     configurations. Both sides use the same random inputs; the
c     outputs are compared element-by-element to tolerance 1e-15.

      program test_stfmm2d
      implicit none

      external stfmm2d, stfmm2d_c
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 2)
      real *8 source(2, ns), targ(2, nt)
      real *8 stoklet(nd, 2, ns)
      real *8 strslet(nd, 2, ns), strsvec(nd, 2, ns)

      real *8 pot_f(nd, 2, ns), pre_f(nd, ns)
      real *8 grad_f(nd, 2, 2, ns)
      real *8 pot_c(nd, 2, ns), pre_c(nd, ns)
      real *8 grad_c(nd, 2, 2, ns)
      real *8 pottarg_f(nd, 2, nt), pretarg_f(nd, nt)
      real *8 gradtarg_f(nd, 2, 2, nt)
      real *8 pottarg_c(nd, 2, nt), pretarg_c(nd, nt)
      real *8 gradtarg_c(nd, 2, 2, nt)

      real *8 eps, dummy, errmax, e
      integer i, j, k, ic
      integer ifstoklet, ifstrslet, ifppreg, ifppregtarg
      integer ier_f, ier_c, nfail

c     Test configurations: stoklet only / strslet only / both,
c     each with ifppreg and ifppregtarg up to 2.
      integer ncases
      parameter (ncases = 6)
      integer cases(4, ncases)

      data cases /
     1  1, 0, 1, 1,
     1  1, 0, 2, 2,
     1  0, 1, 1, 1,
     1  0, 1, 2, 2,
     1  1, 1, 1, 1,
     1  1, 1, 2, 2 /

      nfail = 0
      eps = 1.0d-6

c     Initialize RNG and input data.
      dummy = hkrand(1234)
      do i = 1, ns
         source(1, i) = hkrand(0)
         source(2, i) = hkrand(0)
         do j = 1, nd
            do k = 1, 2
               stoklet(j, k, i) = hkrand(0) - 0.5d0
               strslet(j, k, i) = hkrand(0) - 0.5d0
               strsvec(j, k, i) = hkrand(0) - 0.5d0
            enddo
         enddo
      enddo
      do i = 1, nt
         targ(1, i) = hkrand(0)
         targ(2, i) = hkrand(0)
      enddo

      do ic = 1, ncases
         ifstoklet  = cases(1, ic)
         ifstrslet  = cases(2, ic)
         ifppreg    = cases(3, ic)
         ifppregtarg = cases(4, ic)

c        zero output buffers
         do i = 1, ns
            do j = 1, nd
               pre_f(j, i) = 0.0d0
               pre_c(j, i) = 0.0d0
            enddo
            do k = 1, 2
               do j = 1, nd
                  pot_f(j, k, i) = 0.0d0
                  pot_c(j, k, i) = 0.0d0
               enddo
            enddo
         enddo
         do i = 1, nt
            do j = 1, nd
               pretarg_f(j, i) = 0.0d0
               pretarg_c(j, i) = 0.0d0
            enddo
            do k = 1, 2
               do j = 1, nd
                  pottarg_f(j, k, i) = 0.0d0
                  pottarg_c(j, k, i) = 0.0d0
               enddo
            enddo
         enddo

         ier_f = 0
         ier_c = 0

         call stfmm2d(nd, eps, ns, source, ifstoklet, stoklet,
     1        ifstrslet, strslet, strsvec, ifppreg, pot_f, pre_f,
     2        grad_f, nt, targ, ifppregtarg, pottarg_f, pretarg_f,
     3        gradtarg_f, ier_f)

         call stfmm2d_c(nd, eps, ns, source, ifstoklet, stoklet,
     1        ifstrslet, strslet, strsvec, ifppreg, pot_c, pre_c,
     2        grad_c, nt, targ, ifppregtarg, pottarg_c, pretarg_c,
     3        gradtarg_c, ier_c)

         errmax = 0.0d0
         if (ifppreg .ge. 1) then
            do i = 1, ns
               do k = 1, 2
                  do j = 1, nd
                     e = abs(pot_f(j, k, i) - pot_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
               enddo
            enddo
         endif
         if (ifppreg .ge. 2) then
            do i = 1, ns
               do j = 1, nd
                  e = abs(pre_f(j, i) - pre_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif
         if (ifppregtarg .ge. 1) then
            do i = 1, nt
               do k = 1, 2
                  do j = 1, nd
                     e = abs(pottarg_f(j, k, i) - pottarg_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
               enddo
            enddo
         endif
         if (ifppregtarg .ge. 2) then
            do i = 1, nt
               do j = 1, nd
                  e = abs(pretarg_f(j, i) - pretarg_c(j, i))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         endif

         if (ier_f .ne. ier_c) then
            write(*,1002) ifstoklet, ifstrslet, ifppreg, ifppregtarg,
     1           ier_f, ier_c
            nfail = nfail + 1
         else if (errmax .gt. 1.0d-15) then
            write(*,1001) ifstoklet, ifstrslet, ifppreg, ifppregtarg,
     1           errmax
            nfail = nfail + 1
         else
            write(*,1000) ifstoklet, ifstrslet, ifppreg, ifppregtarg,
     1           errmax
         endif
      enddo

 1000 format(' [ ok ] stfmm2d ifst=',i1,' ifstr=',i1,' ifppr=',i1,
     1     ' ifpprt=',i1,'  errmax=',1pe10.3)
 1001 format(' [FAIL] stfmm2d ifst=',i1,' ifstr=',i1,' ifppr=',i1,
     1     ' ifpprt=',i1,'  errmax=',1pe10.3)
 1002 format(' [FAIL] stfmm2d ifst=',i1,' ifstr=',i1,' ifppr=',i1,
     1     ' ifpprt=',i1,'  ier mismatch ',i6,' vs ',i6)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all stfmm2d cases match'

      end
