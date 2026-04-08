c     test_bhfmm2d.f - differential test for the C translation of
c     src/biharmonic/bhfmm2d.f.
c
c     Compares the Fortran reference bhfmm2d to the C port bhfmm2d_c
c     (and bhfmm2dmain, bhfmm2dpart_direct, bh2dmpalloc indirectly via
c     the full end-to-end driver) over a sweep of (ifcharge, ifdipole,
c     ifpgh, ifpghtarg) configurations. Both sides use the same random
c     inputs; the outputs are compared element-by-element to tolerance
c     1e-15 (the -O0 build makes the computation bit-deterministic so
c     any nonzero error is a bug).
c
c     ier is compared (should be zero on both sides).
c
c     Note: bhfmm2d only supports ifpgh up to 2 (no hessian) so we
c     test ifpgh in {1, 2} only.
c
c     A separate direct test of bh2dmpalloc is also included.

      program test_bhfmm2d
      implicit none

      external bhfmm2d, bhfmm2d_c
      external bh2dmpalloc, bh2dmpalloc_c
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 2)
      real *8 sources(2, ns), targ(2, nt)
c     biharmonic charge has shape (nd,2,*) and dipole shape (nd,3,*)
      complex *16 charge(nd, 2, ns), dip(nd, 3, ns)
      complex *16 pot_f(nd, ns)
      complex *16 grad_f(nd, 3, ns), hess_f(nd, 3, ns)
      complex *16 pot_c(nd, ns)
      complex *16 grad_c(nd, 3, ns), hess_c(nd, 3, ns)
      complex *16 pottarg_f(nd, nt)
      complex *16 gradtarg_f(nd, 3, nt), hesstarg_f(nd, 3, nt)
      complex *16 pottarg_c(nd, nt)
      complex *16 gradtarg_c(nd, 3, nt), hesstarg_c(nd, 3, nt)

      real *8 eps, dummy, errmax, e
      integer i, j, k, ic
      integer ifcharge, ifdipole, iper_f, iper_c, ifpgh, ifpghtarg
      integer ier_f, ier_c, nfail

c     Test configurations: 6 combos of (ifcharge, ifdipole, ifpgh,
c     ifpghtarg), ifpgh in {1,2}.
      integer ncases
      parameter (ncases = 6)
      integer cases(4, ncases)

c     bh2dmpalloc direct test
      integer nlevels_mp
      parameter (nlevels_mp = 3)
      integer laddr_mp(2, 0:nlevels_mp), nterms_mp(0:nlevels_mp)
      integer, allocatable :: iaddr_f(:), iaddr_c(:)
      integer lmptot_f, lmptot_c
      integer nboxes_mp, nd_mp

      data cases /
     1  1, 0, 1, 1,
     1  1, 0, 2, 2,
     1  0, 1, 1, 1,
     1  0, 1, 2, 2,
     1  1, 1, 1, 1,
     1  1, 1, 2, 2 /

      nfail = 0
      eps = 1.0d-6

c     Initialize RNG and input data
      dummy = hkrand(1234)
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
         do j = 1, nd
            do k = 1, 2
               charge(j, k, i) = dcmplx(hkrand(0) - 0.5d0,
     1                                  hkrand(0) - 0.5d0)
            enddo
            do k = 1, 3
               dip(j, k, i) = dcmplx(hkrand(0) - 0.5d0,
     1                               hkrand(0) - 0.5d0)
            enddo
         enddo
      enddo
      do i = 1, nt
         targ(1, i) = hkrand(0)
         targ(2, i) = hkrand(0)
      enddo

c     ================================================================
c     Sweep test cases for bhfmm2d.
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
               pot_c(j, i)  = (0.0d0, 0.0d0)
            enddo
            do k = 1, 3
               do j = 1, nd
                  grad_f(j, k, i) = (0.0d0, 0.0d0)
                  grad_c(j, k, i) = (0.0d0, 0.0d0)
                  hess_f(j, k, i) = (0.0d0, 0.0d0)
                  hess_c(j, k, i) = (0.0d0, 0.0d0)
               enddo
            enddo
         enddo
         do i = 1, nt
            do j = 1, nd
               pottarg_f(j, i)  = (0.0d0, 0.0d0)
               pottarg_c(j, i)  = (0.0d0, 0.0d0)
            enddo
            do k = 1, 3
               do j = 1, nd
                  gradtarg_f(j, k, i) = (0.0d0, 0.0d0)
                  gradtarg_c(j, k, i) = (0.0d0, 0.0d0)
                  hesstarg_f(j, k, i) = (0.0d0, 0.0d0)
                  hesstarg_c(j, k, i) = (0.0d0, 0.0d0)
               enddo
            enddo
         enddo

         iper_f = 0
         iper_c = 0
         ier_f = 0
         ier_c = 0

         call bhfmm2d(nd, eps, ns, sources, ifcharge, charge,
     1        ifdipole, dip, iper_f, ifpgh, pot_f, grad_f,
     2        hess_f, nt, targ, ifpghtarg, pottarg_f, gradtarg_f,
     3        hesstarg_f, ier_f)

         call bhfmm2d_c(nd, eps, ns, sources, ifcharge, charge,
     1        ifdipole, dip, iper_c, ifpgh, pot_c, grad_c,
     2        hess_c, nt, targ, ifpghtarg, pottarg_c, gradtarg_c,
     3        hesstarg_c, ier_c)

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
               do k = 1, 3
                  do j = 1, nd
                     e = abs(grad_f(j, k, i) - grad_c(j, k, i))
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
               do k = 1, 3
                  do j = 1, nd
                     e = abs(gradtarg_f(j, k, i) - gradtarg_c(j, k, i))
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

 1000 format(' [ ok ] bhfmm2d ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  errmax=',1pe10.3)
 1001 format(' [FAIL] bhfmm2d ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  errmax=',1pe10.3)
 1002 format(' [FAIL] bhfmm2d ifc=',i1,' ifd=',i1,' ifpgh=',i1,
     1     ' ifpght=',i1,'  ier mismatch ',i6,' vs ',i6)

c     ================================================================
c     Direct test of bh2dmpalloc with a hand-built laddr/nterms.
c     ================================================================
      nd_mp = 3
      nterms_mp(0) = 4
      nterms_mp(1) = 6
      nterms_mp(2) = 8
      nterms_mp(3) = 10
      laddr_mp(1, 0) = 1
      laddr_mp(2, 0) = 1
      laddr_mp(1, 1) = 2
      laddr_mp(2, 1) = 5
      laddr_mp(1, 2) = 6
      laddr_mp(2, 2) = 13
      laddr_mp(1, 3) = 14
      laddr_mp(2, 3) = 18
      nboxes_mp = 18
      allocate(iaddr_f(2*nboxes_mp), iaddr_c(2*nboxes_mp))
      do i = 1, 2*nboxes_mp
         iaddr_f(i) = -1
         iaddr_c(i) = -1
      enddo
      lmptot_f = -1
      lmptot_c = -1
      call bh2dmpalloc(nd_mp, laddr_mp, iaddr_f, nlevels_mp, lmptot_f,
     1     nterms_mp)
      call bh2dmpalloc_c(nd_mp, laddr_mp, iaddr_c, nlevels_mp, lmptot_c,
     1     nterms_mp)
      if (lmptot_f .ne. lmptot_c) then
         write(*,*) ' [FAIL] bh2dmpalloc lmptot mismatch: ',
     1        lmptot_f, lmptot_c
         nfail = nfail + 1
      else
         j = 0
         do i = 1, 2*nboxes_mp
            if (iaddr_f(i) .ne. iaddr_c(i)) j = j + 1
         enddo
         if (j .ne. 0) then
            write(*,*) ' [FAIL] bh2dmpalloc iaddr mismatches: ', j
            nfail = nfail + 1
         else
            write(*,*) ' [ ok ] bh2dmpalloc  lmptot=', lmptot_f
         endif
      endif
      deallocate(iaddr_f, iaddr_c)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all bhfmm2d cases match'

      end
