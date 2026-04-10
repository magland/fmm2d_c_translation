c     test_hfmm2d.f - differential test for the C translation of
c     src/helmholtz/hfmm2d.f.
c
c     Tests hfmm2d and h2dmpalloc by calling both Fortran and C
c     versions on identical inputs and comparing outputs.

      program test_hfmm2d
      implicit none

      external hfmm2d, hfmm2d_c
      external h2dmpalloc, h2dmpalloc_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax=200, ntmax=100, ndmax=2)

      integer nd, ns, nt, i, j, ii
      real *8 eps
      complex *16 zk
      real *8 sources(2,nsmax), targ(2,ntmax)
      complex *16 charge(ndmax,nsmax), dipstr(ndmax,nsmax)
      real *8 dipvec(ndmax,2,nsmax)
      integer ifcharge, ifdipole, ifpgh, ifpghtarg, iper
      integer ier_f, ier_c

      complex *16 pot_f(ndmax,nsmax), pot_c(ndmax,nsmax)
      complex *16 grad_f(ndmax,2,nsmax), grad_c(ndmax,2,nsmax)
      complex *16 hess_f(ndmax,3,nsmax), hess_c(ndmax,3,nsmax)
      complex *16 pott_f(ndmax,ntmax), pott_c(ndmax,ntmax)
      complex *16 gradt_f(ndmax,2,ntmax)
      complex *16 gradt_c(ndmax,2,ntmax)
      complex *16 hesst_f(ndmax,3,ntmax)
      complex *16 hesst_c(ndmax,3,ntmax)

      integer nfail
      real *8 dummy, errmax, e

      nfail = 0
      dummy = hkrand(1234)

      nd = ndmax
      ns = nsmax
      nt = ntmax
      eps = 1.0d-6
      zk = dcmplx(3.0d0, 0.05d0)

      do i = 1, ns
         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)
      enddo
      do j = 1, nt
         targ(1,j) = hkrand(0) + 2.0d0
         targ(2,j) = hkrand(0) + 2.0d0
      enddo
      do i = 1, ns
         do ii = 1, nd
            charge(ii,i) = dcmplx(hkrand(0), hkrand(0))
            dipstr(ii,i) = dcmplx(hkrand(0), hkrand(0))
            dipvec(ii,1,i) = hkrand(0)
            dipvec(ii,2,i) = hkrand(0)
         enddo
      enddo

c     Test charge-only, pot only
      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 1
      call zero16(pot_f, nd*ns)
      call zero16(pot_c, nd*ns)
      call zero16(pott_f, nd*nt)
      call zero16(pott_c, nd*nt)
      ier_f = 0
      ier_c = 0
      iper = 0
      call hfmm2d(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_f,
     1     grad_f, hess_f, nt, targ,
     1     ifpghtarg, pott_f, gradt_f, hesst_f,
     1     ier_f)
      iper = 0
      call hfmm2d_c(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_c,
     1     grad_c, hess_c, nt, targ,
     1     ifpghtarg, pott_c, gradt_c, hesst_c,
     1     ier_c)

      errmax = 0.0d0
      do i = 1, ns
         do ii = 1, nd
            e = cdabs(pot_f(ii,i) - pot_c(ii,i))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pott_f(ii,j) - pott_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,1000) 'hfmm2d ifc=1 ifd=0',errmax
         nfail = nfail + 1
      else
         write(*,1001) 'hfmm2d ifc=1 ifd=0'
      endif

c     Test dipole-only, pot+grad
      ifcharge = 0
      ifdipole = 1
      ifpgh = 2
      ifpghtarg = 2
      call zero16(pot_f, nd*ns)
      call zero16(pot_c, nd*ns)
      call zero16(grad_f, nd*2*ns)
      call zero16(grad_c, nd*2*ns)
      call zero16(pott_f, nd*nt)
      call zero16(pott_c, nd*nt)
      call zero16(gradt_f, nd*2*nt)
      call zero16(gradt_c, nd*2*nt)
      iper = 0
      call hfmm2d(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_f,
     1     grad_f, hess_f, nt, targ,
     1     ifpghtarg, pott_f, gradt_f, hesst_f,
     1     ier_f)
      iper = 0
      call hfmm2d_c(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_c,
     1     grad_c, hess_c, nt, targ,
     1     ifpghtarg, pott_c, gradt_c, hesst_c,
     1     ier_c)

      errmax = 0.0d0
      do i = 1, ns
         do ii = 1, nd
            e = cdabs(pot_f(ii,i) - pot_c(ii,i))
            if (e .gt. errmax) errmax = e
            e = cdabs(grad_f(ii,1,i)-grad_c(ii,1,i))
            if (e .gt. errmax) errmax = e
            e = cdabs(grad_f(ii,2,i)-grad_c(ii,2,i))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pott_f(ii,j)-pott_c(ii,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(gradt_f(ii,1,j)-gradt_c(ii,1,j))
            if (e .gt. errmax) errmax = e
            e = cdabs(gradt_f(ii,2,j)-gradt_c(ii,2,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,1000) 'hfmm2d ifc=0 ifd=1',errmax
         nfail = nfail + 1
      else
         write(*,1001) 'hfmm2d ifc=0 ifd=1'
      endif

c     Test charge+dipole, pot+grad+hess
      ifcharge = 1
      ifdipole = 1
      ifpgh = 3
      ifpghtarg = 3
      call zero16(pot_f, nd*ns)
      call zero16(pot_c, nd*ns)
      call zero16(grad_f, nd*2*ns)
      call zero16(grad_c, nd*2*ns)
      call zero16(hess_f, nd*3*ns)
      call zero16(hess_c, nd*3*ns)
      call zero16(pott_f, nd*nt)
      call zero16(pott_c, nd*nt)
      call zero16(gradt_f, nd*2*nt)
      call zero16(gradt_c, nd*2*nt)
      call zero16(hesst_f, nd*3*nt)
      call zero16(hesst_c, nd*3*nt)
      iper = 0
      call hfmm2d(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_f,
     1     grad_f, hess_f, nt, targ,
     1     ifpghtarg, pott_f, gradt_f, hesst_f,
     1     ier_f)
      iper = 0
      call hfmm2d_c(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_c,
     1     grad_c, hess_c, nt, targ,
     1     ifpghtarg, pott_c, gradt_c, hesst_c,
     1     ier_c)

      errmax = 0.0d0
      do i = 1, ns
         do ii = 1, nd
            e = cdabs(pot_f(ii,i) - pot_c(ii,i))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      do j = 1, nt
         do ii = 1, nd
            e = cdabs(pott_f(ii,j)-pott_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,1000)'hfmm2d ifc=1 ifd=1',errmax
         nfail = nfail + 1
      else
         write(*,1001) 'hfmm2d ifc=1 ifd=1'
      endif

c     ---- h2dmpalloc test ----
c     Test that h2dmpalloc produces the same iaddr and lmptot
      call test_mpalloc(nfail)

 1000 format(' [FAIL] ',a,' errmax=',1pe12.5)
 1001 format(' [ ok ] ',a)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' checks differ'
         stop 1
      endif
      write(*,*) 'PASS: hfmm2d all tests'

      end

      subroutine zero16(buf, n)
      implicit none
      integer n, i
      complex *16 buf(n)
      do i = 1, n
         buf(i) = (0.0d0, 0.0d0)
      enddo
      return
      end

      subroutine test_mpalloc(nfail)
      implicit none
      external h2dmpalloc, h2dmpalloc_c
      integer nfail
      integer nd, nlevels, nboxes
      parameter (nd=2, nlevels=3, nboxes=20)
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer iaddr_f(4,nboxes), iaddr_c(4,nboxes)
      integer lmptot_f, lmptot_c
      integer i, j, errcount

      laddr(1,0) = 1
      laddr(2,0) = 1
      laddr(1,1) = 2
      laddr(2,1) = 5
      laddr(1,2) = 6
      laddr(2,2) = 13
      laddr(1,3) = 14
      laddr(2,3) = 20
      nterms(0) = 5
      nterms(1) = 8
      nterms(2) = 12
      nterms(3) = 15

      call h2dmpalloc(nd, laddr, iaddr_f,
     1     nlevels, lmptot_f, nterms)
      call h2dmpalloc_c(nd, laddr, iaddr_c,
     1     nlevels, lmptot_c, nterms)

      errcount = 0
      do j = 1, nboxes
         do i = 1, 4
            if (iaddr_f(i,j) .ne. iaddr_c(i,j)) then
               errcount = errcount + 1
            endif
         enddo
      enddo
      if (lmptot_f .ne. lmptot_c) errcount = errcount+1

      if (errcount .gt. 0) then
         write(*,*) ' [FAIL] h2dmpalloc ',
     1        errcount, ' diffs'
         write(*,*) '  lmptot f/c=',lmptot_f,lmptot_c
         nfail = nfail + 1
      else
         write(*,*) ' [ ok ] h2dmpalloc  lmptot=',
     1        lmptot_f
      endif
      return
      end
