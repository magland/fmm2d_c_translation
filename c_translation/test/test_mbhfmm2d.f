c     test_mbhfmm2d.f - differential test for the C translation of
c     src/modified-biharmonic/mbhfmm2d.f.
c
c     Compares the Fortran reference mbhfmm2d to the C port mbhfmm2d_c
c     (and mbhfmm2dmain, mbhfmm2dmps_direct_vec, mbh2dmpalloc indirectly
c     via the full end-to-end driver) over a sweep of source-type and
c     ifpgh/ifpghtarg configurations. Both sides use the same random
c     inputs; outputs are compared element-by-element to tolerance
c     1e-15 (the -O0 build makes the computation bit-deterministic so
c     any nonzero error is a bug).
c
c     A separate direct test of mbh2dmpalloc is also included.

      program test_mbhfmm2d
      implicit none

      external mbhfmm2d, mbhfmm2d_c
      external mbh2dmpalloc, mbh2dmpalloc_c
      external hkrand
      real *8 hkrand

      integer ns, nt, nd
      parameter (ns = 200, nt = 100, nd = 2)
      real *8 sources(2, ns), targ(2, nt)
      real *8 charge(nd, ns)
      real *8 dipstr(nd, ns), dipvec(nd, 2, ns)
      real *8 quadstr(nd, ns), quadvec(nd, 3, ns)
      real *8 octstr(nd, ns), octvec(nd, 4, ns)

      real *8 pot_f(nd, ns), pot_c(nd, ns)
      real *8 grad_f(nd, 2, ns), grad_c(nd, 2, ns)
      real *8 hess_f(nd, 3, ns), hess_c(nd, 3, ns)
      real *8 pottarg_f(nd, nt), pottarg_c(nd, nt)
      real *8 gradtarg_f(nd, 2, nt), gradtarg_c(nd, 2, nt)
      real *8 hesstarg_f(nd, 3, nt), hesstarg_c(nd, 3, nt)

      real *8 eps, beta, dummy, errmax, e
      integer i, j, k, ic
      integer ifcharge, ifdipole, ifquad, ifoct
      integer iper_f, iper_c, ifpgh, ifpghtarg
      integer ier_f, ier_c, nfail

c     Test configurations: (ifcharge,ifdipole,ifquad,ifoct,ifpgh,ifpghtarg)
      integer ncases
      parameter (ncases = 6)
      integer cases(6, ncases)

c     mbh2dmpalloc direct test
      integer nlevels_mp
      parameter (nlevels_mp = 3)
      integer laddr_mp(2, 0:nlevels_mp), nterms_mp(0:nlevels_mp)
      integer, allocatable :: iaddr_f(:), iaddr_c(:)
      integer lmptot_f, lmptot_c
      integer nboxes_mp, nd_mp

      data cases /
     1  1, 0, 0, 0, 1, 1,
     1  1, 0, 0, 0, 2, 2,
     1  1, 0, 0, 0, 3, 3,
     1  1, 1, 0, 0, 1, 1,
     1  1, 1, 1, 1, 2, 2,
     1  1, 1, 1, 1, 3, 3 /

      nfail = 0
      eps = 1.0d-6
      beta = 1.5d0

c     Initialize RNG and input data
      dummy = hkrand(1234)
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
         do j = 1, nd
            charge(j, i) = hkrand(0) - 0.5d0
            dipstr(j, i) = hkrand(0) - 0.5d0
            dipvec(j, 1, i) = hkrand(0) - 0.5d0
            dipvec(j, 2, i) = hkrand(0) - 0.5d0
            quadstr(j, i) = hkrand(0) - 0.5d0
            quadvec(j, 1, i) = hkrand(0) - 0.5d0
            quadvec(j, 2, i) = hkrand(0) - 0.5d0
            quadvec(j, 3, i) = hkrand(0) - 0.5d0
            octstr(j, i) = hkrand(0) - 0.5d0
            octvec(j, 1, i) = hkrand(0) - 0.5d0
            octvec(j, 2, i) = hkrand(0) - 0.5d0
            octvec(j, 3, i) = hkrand(0) - 0.5d0
            octvec(j, 4, i) = hkrand(0) - 0.5d0
         enddo
      enddo
      do i = 1, nt
         targ(1, i) = hkrand(0)
         targ(2, i) = hkrand(0)
      enddo

c     ================================================================
c     Sweep test cases for mbhfmm2d.
c     ================================================================
      do ic = 1, ncases
         ifcharge  = cases(1, ic)
         ifdipole  = cases(2, ic)
         ifquad    = cases(3, ic)
         ifoct     = cases(4, ic)
         ifpgh     = cases(5, ic)
         ifpghtarg = cases(6, ic)

c        zero output buffers
         do i = 1, ns
            do j = 1, nd
               pot_f(j, i)  = 0.0d0
               pot_c(j, i)  = 0.0d0
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

         iper_f = 0
         iper_c = 0
         ier_f = 0
         ier_c = 0

         call mbhfmm2d(nd, eps, beta, ns, sources, ifcharge, charge,
     1        ifdipole, dipstr, dipvec, ifquad, quadstr, quadvec,
     2        ifoct, octstr, octvec, iper_f, ifpgh,
     3        pot_f, grad_f, hess_f, nt, targ, ifpghtarg,
     4        pottarg_f, gradtarg_f, hesstarg_f, ier_f)

         call mbhfmm2d_c(nd, eps, beta, ns, sources, ifcharge, charge,
     1        ifdipole, dipstr, dipvec, ifquad, quadstr, quadvec,
     2        ifoct, octstr, octvec, iper_c, ifpgh,
     3        pot_c, grad_c, hess_c, nt, targ, ifpghtarg,
     4        pottarg_c, gradtarg_c, hesstarg_c, ier_c)

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
                     e = abs(gradtarg_f(j, k, i) - gradtarg_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
               enddo
            enddo
         endif
         if (ifpghtarg .ge. 3) then
            do i = 1, nt
               do k = 1, 3
                  do j = 1, nd
                     e = abs(hesstarg_f(j, k, i) - hesstarg_c(j, k, i))
                     if (e .gt. errmax) errmax = e
                  enddo
               enddo
            enddo
         endif

         if (ier_f .ne. ier_c) then
            write(*,1002) ifcharge, ifdipole, ifquad, ifoct,
     1           ifpgh, ifpghtarg, ier_f, ier_c
            nfail = nfail + 1
         else if (errmax .gt. 1.0d-15) then
            write(*,1001) ifcharge, ifdipole, ifquad, ifoct,
     1           ifpgh, ifpghtarg, errmax
            nfail = nfail + 1
         else
            write(*,1000) ifcharge, ifdipole, ifquad, ifoct,
     1           ifpgh, ifpghtarg, errmax
         endif
      enddo

 1000 format(' [ ok ] mbhfmm2d ifc=',i1,' ifd=',i1,' ifq=',i1,
     1     ' ifo=',i1,' ifpgh=',i1,' ifpght=',i1,'  errmax=',1pe10.3)
 1001 format(' [FAIL] mbhfmm2d ifc=',i1,' ifd=',i1,' ifq=',i1,
     1     ' ifo=',i1,' ifpgh=',i1,' ifpght=',i1,'  errmax=',1pe10.3)
 1002 format(' [FAIL] mbhfmm2d ifc=',i1,' ifd=',i1,' ifq=',i1,
     1     ' ifo=',i1,' ifpgh=',i1,' ifpght=',i1,
     2     '  ier mismatch ',i6,' vs ',i6)

c     ================================================================
c     Direct test of mbh2dmpalloc with a hand-built laddr/nterms.
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
      allocate(iaddr_f(4*nboxes_mp), iaddr_c(4*nboxes_mp))
      do i = 1, 4*nboxes_mp
         iaddr_f(i) = -1
         iaddr_c(i) = -1
      enddo
      lmptot_f = -1
      lmptot_c = -1
      call mbh2dmpalloc(nd_mp, laddr_mp, iaddr_f, nlevels_mp,
     1     lmptot_f, nterms_mp)
      call mbh2dmpalloc_c(nd_mp, laddr_mp, iaddr_c, nlevels_mp,
     1     lmptot_c, nterms_mp)
      if (lmptot_f .ne. lmptot_c) then
         write(*,*) ' [FAIL] mbh2dmpalloc lmptot mismatch: ',
     1        lmptot_f, lmptot_c
         nfail = nfail + 1
      else
         j = 0
         do i = 1, 4*nboxes_mp
            if (iaddr_f(i) .ne. iaddr_c(i)) j = j + 1
         enddo
         if (j .ne. 0) then
            write(*,*) ' [FAIL] mbh2dmpalloc iaddr mismatches: ', j
            nfail = nfail + 1
         else
            write(*,*) ' [ ok ] mbh2dmpalloc  lmptot=', lmptot_f
         endif
      endif
      deallocate(iaddr_f, iaddr_c)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all mbhfmm2d cases match'

      end
