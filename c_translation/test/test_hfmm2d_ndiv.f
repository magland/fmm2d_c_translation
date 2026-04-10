c     test_hfmm2d_ndiv.f - differential test for hfmm2d_ndiv.c

      program test_hfmm2d_ndiv
      implicit none

      external hfmm2d_ndiv, hfmm2d_ndiv_c
      real *8 hkrand
      external hkrand

      integer ns, nt, nd
      parameter (ns=200, nt=100, nd=2)

      real *8 sources(2,ns), targ(2,nt)
      complex *16 charge(nd,ns), dipstr(nd,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 zk
      real *8 eps
      integer ifcharge, ifdipole, ifpgh, ifpghtarg
      integer iper, ndiv, idivflag, ifnear, ier_f, ier_c

      complex *16 pot_f(nd,ns), pot_c(nd,ns)
      complex *16 grad_f(nd,2,ns), grad_c(nd,2,ns)
      complex *16 hess_f(nd,3,ns), hess_c(nd,3,ns)
      complex *16 pott_f(nd,nt), pott_c(nd,nt)
      complex *16 gradt_f(nd,2,nt), gradt_c(nd,2,nt)
      complex *16 hesst_f(nd,3,nt), hesst_c(nd,3,nt)
      real *8 timeinfo_f(8), timeinfo_c(8)

      integer nfail, i, j, ii
      real *8 dummy, errmax, e

      nfail = 0
      dummy = hkrand(1234)

      eps = 1.0d-6
      zk = dcmplx(3.0d0, 0.05d0)
      ndiv = 15
      idivflag = 0
      ifnear = 1

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

c     Test charge+dipole, pot+grad
      ifcharge = 1
      ifdipole = 1
      ifpgh = 2
      ifpghtarg = 2
      call zz16(pot_f, nd*ns)
      call zz16(pot_c, nd*ns)
      call zz16(grad_f, nd*2*ns)
      call zz16(grad_c, nd*2*ns)
      call zz16(pott_f, nd*nt)
      call zz16(pott_c, nd*nt)
      call zz16(gradt_f, nd*2*nt)
      call zz16(gradt_c, nd*2*nt)
      iper = 0
      call hfmm2d_ndiv(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_f,
     1     grad_f, hess_f, nt, targ,
     1     ifpghtarg, pott_f, gradt_f, hesst_f,
     1     ndiv, idivflag, ifnear, timeinfo_f,
     1     ier_f)
      iper = 0
      call hfmm2d_ndiv_c(nd, eps, zk, ns, sources,
     1     ifcharge, charge, ifdipole, dipstr,
     1     dipvec, iper, ifpgh, pot_c,
     1     grad_c, hess_c, nt, targ,
     1     ifpghtarg, pott_c, gradt_c, hesst_c,
     1     ndiv, idivflag, ifnear, timeinfo_c,
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
            e = cdabs(pott_f(ii,j) - pott_c(ii,j))
            if (e .gt. errmax) errmax = e
         enddo
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,*) '[FAIL] hfmm2d_ndiv err=',errmax
         nfail = nfail + 1
      else
         write(*,*) '[ ok ] hfmm2d_ndiv'
      endif

      if (nfail .ne. 0) then
         write(*,*) 'FAIL:', nfail, ' checks differ'
         stop 1
      endif
      write(*,*) 'PASS: hfmm2d_ndiv all tests'

      end

      subroutine zz16(buf, n)
      implicit none
      integer n, i
      complex *16 buf(n)
      do i = 1, n
         buf(i) = (0.0d0, 0.0d0)
      enddo
      return
      end
