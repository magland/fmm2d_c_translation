c     test_helmkernels2d.f - differential test for the C translation of
c     src/helmholtz/helmkernels2d.f.

      program test_helmkernels2d
      implicit none

      external h2d_directcp,  h2d_directcp_c
      external h2d_directcg,  h2d_directcg_c
      external h2d_directch,  h2d_directch_c
      external h2d_directdp,  h2d_directdp_c
      external h2d_directdg,  h2d_directdg_c
      external h2d_directdh,  h2d_directdh_c
      external h2d_directcdp, h2d_directcdp_c
      external h2d_directcdg, h2d_directcdg_c
      external h2d_directcdh, h2d_directcdh_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 30, ntmax = 20, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      complex *16 wavek
      complex *16 charge(ndmax,nsmax), dipstr(ndmax,nsmax)
      real *8 dipvec(ndmax, 2, nsmax)
      real *8 thresh

      complex *16, allocatable :: pot_f(:,:), pot_c(:,:)
      complex *16, allocatable :: grad_f(:,:,:)
      complex *16, allocatable :: grad_c(:,:,:)
      complex *16, allocatable :: hess_f(:,:,:)
      complex *16, allocatable :: hess_c(:,:,:)
      complex *16, allocatable :: charge_l(:,:)
      complex *16, allocatable :: dipstr_l(:,:)
      real *8, allocatable :: dipvec_l(:,:,:)

      integer ns, nt, nds(2)
      data nds / 1, 3 /

      integer ind, nd, i, j, ii, k
      integer nfail
      real *8 dummy, errmax, e

      ns = nsmax
      nt = ntmax
      thresh = 1.0d-15
      nfail = 0
      wavek = dcmplx(3.5d0, 0.1d0)

      dummy = hkrand(1234)

      do i = 1, ns
         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)
      enddo
      do j = 1, nt
         targ(1,j) = 2.0d0 + hkrand(0)
         targ(2,j) = 2.0d0 + hkrand(0)
      enddo
      do i = 1, ns
         do ii = 1, ndmax
            charge(ii,i) = dcmplx(2*hkrand(0)-1,
     1                            2*hkrand(0)-1)
            dipstr(ii,i) = dcmplx(2*hkrand(0)-1,
     1                            2*hkrand(0)-1)
            dipvec(ii,1,i) = 2*hkrand(0)-1
            dipvec(ii,2,i) = 2*hkrand(0)-1
         enddo
      enddo

      do ind = 1, 2
         nd = nds(ind)

         allocate(charge_l(nd,ns))
         allocate(dipstr_l(nd,ns))
         allocate(dipvec_l(nd,2,ns))
         do i = 1, ns
            do ii = 1, nd
               charge_l(ii,i) = charge(ii,i)
               dipstr_l(ii,i) = dipstr(ii,i)
               dipvec_l(ii,1,i) = dipvec(ii,1,i)
               dipvec_l(ii,2,i) = dipvec(ii,2,i)
            enddo
         enddo

         allocate(pot_f(nd,nt), pot_c(nd,nt))
         allocate(grad_f(nd,2,nt), grad_c(nd,2,nt))
         allocate(hess_f(nd,3,nt), hess_c(nd,3,nt))

c        ---- h2d_directcp ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call h2d_directcp(nd, wavek, sources,
     1        ns, charge_l, targ, nt, pot_f,
     1        thresh)
         call h2d_directcp_c(nd, wavek, sources,
     1        ns, charge_l, targ, nt, pot_c,
     1        thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directcp', nd, nfail)

c        ---- h2d_directcg ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call zerobuf(grad_f, nd*2*nt)
         call zerobuf(grad_c, nd*2*nt)
         call h2d_directcg(nd, wavek, sources,
     1        ns, charge_l, targ, nt, pot_f,
     1        grad_f, thresh)
         call h2d_directcg_c(nd, wavek, sources,
     1        ns, charge_l, targ, nt, pot_c,
     1        grad_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directcg:p', nd, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1        'h2d_directcg:g', nd, nfail)

c        ---- h2d_directch ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call zerobuf(grad_f, nd*2*nt)
         call zerobuf(grad_c, nd*2*nt)
         call zerobuf(hess_f, nd*3*nt)
         call zerobuf(hess_c, nd*3*nt)
         call h2d_directch(nd, wavek, sources,
     1        ns, charge_l, targ, nt, pot_f,
     1        grad_f, hess_f, thresh)
         call h2d_directch_c(nd, wavek, sources,
     1        ns, charge_l, targ, nt, pot_c,
     1        grad_c, hess_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directch:p', nd, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1        'h2d_directch:g', nd, nfail)
         call checkerr(hess_f, hess_c, nd*3*nt,
     1        'h2d_directch:h', nd, nfail)

c        ---- h2d_directdp ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call h2d_directdp(nd, wavek, sources,
     1        ns, dipstr_l, dipvec_l, targ,
     1        nt, pot_f, thresh)
         call h2d_directdp_c(nd, wavek, sources,
     1        ns, dipstr_l, dipvec_l, targ,
     1        nt, pot_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directdp', nd, nfail)

c        ---- h2d_directdg ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call zerobuf(grad_f, nd*2*nt)
         call zerobuf(grad_c, nd*2*nt)
         call h2d_directdg(nd, wavek, sources,
     1        ns, dipstr_l, dipvec_l, targ,
     1        nt, pot_f, grad_f, thresh)
         call h2d_directdg_c(nd, wavek, sources,
     1        ns, dipstr_l, dipvec_l, targ,
     1        nt, pot_c, grad_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directdg:p', nd, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1        'h2d_directdg:g', nd, nfail)

c        ---- h2d_directdh ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call zerobuf(grad_f, nd*2*nt)
         call zerobuf(grad_c, nd*2*nt)
         call zerobuf(hess_f, nd*3*nt)
         call zerobuf(hess_c, nd*3*nt)
         call h2d_directdh(nd, wavek, sources,
     1        ns, dipstr_l, dipvec_l, targ,
     1        nt, pot_f, grad_f, hess_f, thresh)
         call h2d_directdh_c(nd, wavek, sources,
     1        ns, dipstr_l, dipvec_l, targ,
     1        nt, pot_c, grad_c, hess_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directdh:p', nd, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1        'h2d_directdh:g', nd, nfail)
         call checkerr(hess_f, hess_c, nd*3*nt,
     1        'h2d_directdh:h', nd, nfail)

c        ---- h2d_directcdp ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call h2d_directcdp(nd, wavek, sources,
     1        ns, charge_l, dipstr_l, dipvec_l,
     1        targ, nt, pot_f, thresh)
         call h2d_directcdp_c(nd, wavek, sources,
     1        ns, charge_l, dipstr_l, dipvec_l,
     1        targ, nt, pot_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directcdp', nd, nfail)

c        ---- h2d_directcdg ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call zerobuf(grad_f, nd*2*nt)
         call zerobuf(grad_c, nd*2*nt)
         call h2d_directcdg(nd, wavek, sources,
     1        ns, charge_l, dipstr_l, dipvec_l,
     1        targ, nt, pot_f, grad_f, thresh)
         call h2d_directcdg_c(nd, wavek, sources,
     1        ns, charge_l, dipstr_l, dipvec_l,
     1        targ, nt, pot_c, grad_c, thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directcdg:p', nd, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1        'h2d_directcdg:g', nd, nfail)

c        ---- h2d_directcdh ----
         call zerobuf(pot_f, nd*nt)
         call zerobuf(pot_c, nd*nt)
         call zerobuf(grad_f, nd*2*nt)
         call zerobuf(grad_c, nd*2*nt)
         call zerobuf(hess_f, nd*3*nt)
         call zerobuf(hess_c, nd*3*nt)
         call h2d_directcdh(nd, wavek, sources,
     1        ns, charge_l, dipstr_l, dipvec_l,
     1        targ, nt, pot_f, grad_f, hess_f,
     1        thresh)
         call h2d_directcdh_c(nd, wavek, sources,
     1        ns, charge_l, dipstr_l, dipvec_l,
     1        targ, nt, pot_c, grad_c, hess_c,
     1        thresh)
         call checkerr(pot_f, pot_c, nd*nt,
     1        'h2d_directcdh:p', nd, nfail)
         call checkerr(grad_f, grad_c, nd*2*nt,
     1        'h2d_directcdh:g', nd, nfail)
         call checkerr(hess_f, hess_c, nd*3*nt,
     1        'h2d_directcdh:h', nd, nfail)

         deallocate(pot_f, pot_c)
         deallocate(grad_f, grad_c)
         deallocate(hess_f, hess_c)
         deallocate(charge_l, dipstr_l, dipvec_l)
      enddo

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' checks differ'
         stop 1
      endif
      write(*,*) 'PASS: helmkernels2d all tests'

      end

c     Helper: zero a complex buffer
      subroutine zerobuf(buf, n)
      implicit none
      integer n, i
      complex *16 buf(n)
      do i = 1, n
         buf(i) = (0.0d0, 0.0d0)
      enddo
      return
      end

c     Helper: check max abs diff
      subroutine checkerr(f, c, n, name, nd, nfail)
      implicit none
      integer n, nd, nfail, i
      complex *16 f(n), c(n)
      character*(*) name
      real *8 errmax, e
      errmax = 0.0d0
      do i = 1, n
         e = cdabs(f(i) - c(i))
         if (e .gt. errmax) errmax = e
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,1000) name, nd, errmax
         nfail = nfail + 1
      else
         write(*,1001) name, nd
      endif
 1000 format(' [FAIL] ',a,' nd=',i2,
     1   ' errmax=',1pe12.5)
 1001 format(' [ ok ] ',a,' nd=',i2)
      return
      end
