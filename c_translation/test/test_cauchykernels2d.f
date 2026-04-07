c     test_cauchykernels2d.f - differential test for the C translation of
c     src/laplace/cauchykernels2d.f.
c
c     For each kernel, allocates output buffers with leading dimension
c     equal to the actual nd (NOT ndmax) so that the routines see the
c     storage stride they expect. Calls the Fortran reference and the
c     _c_ translation on identical inputs and asserts that every
c     complex element agrees in absolute value.

      program test_cauchykernels2d
      implicit none

      external c2d_directcp,  c2d_directcp_c
      external c2d_directcg,  c2d_directcg_c
      external c2d_directch,  c2d_directch_c
      external c2d_directdp,  c2d_directdp_c
      external c2d_directdg,  c2d_directdg_c
      external c2d_directdh,  c2d_directdh_c
      external c2d_directcdp, c2d_directcdp_c
      external c2d_directcdg, c2d_directcdg_c
      external c2d_directcdh, c2d_directcdh_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 50, ntmax = 40, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      complex *16 charge(ndmax, nsmax), dipstr(ndmax, nsmax)

c     Per-test allocatable buffers, sized to (nd, nt) per iteration
      complex *16, allocatable :: pot_f(:,:),  pot_c(:,:)
      complex *16, allocatable :: grad_f(:,:), grad_c(:,:)
      complex *16, allocatable :: hess_f(:,:), hess_c(:,:)
c     Local copies of charge / dipstr resized to (nd, ns) for the call
      complex *16, allocatable :: charge_loc(:,:), dipstr_loc(:,:)

      integer ns, nt
      integer nds(2)
      data nds / 1, 3 /

      integer ind, nd, i, j, ii
      integer nfail
      real *8 thresh, dummy, errmax, e

      ns = nsmax
      nt = ntmax
      thresh = 1.0d-15
      nfail = 0

c     seed hkrand once
      dummy = hkrand(1234)

c     fill sources, targets, charge, dipstr at maximum sizes
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
      enddo
      do j = 1, nt
         targ(1, j) = hkrand(0)
         targ(2, j) = hkrand(0)
      enddo
      do i = 1, ns
         do ii = 1, ndmax
            charge(ii, i) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                             2.0d0*hkrand(0)-1.0d0)
            dipstr(ii, i) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                             2.0d0*hkrand(0)-1.0d0)
         enddo
      enddo

c     ============================================================
c     Loop over the two nd values and run all 9 kernels for each
c     ============================================================
      do ind = 1, 2
         nd = nds(ind)

c        Resize charge / dipstr to (nd, ns) so the routine's view of
c        their storage matches their declared shape exactly.
         allocate(charge_loc(nd, ns), dipstr_loc(nd, ns))
         do i = 1, ns
            do ii = 1, nd
               charge_loc(ii, i) = charge(ii, i)
               dipstr_loc(ii, i) = dipstr(ii, i)
            enddo
         enddo

c        Allocate output buffers at the correct (nd, nt) shape.
         allocate(pot_f(nd, nt),  pot_c(nd, nt))
         allocate(grad_f(nd, nt), grad_c(nd, nt))
         allocate(hess_f(nd, nt), hess_c(nd, nt))

c        ---- c2d_directcp ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call c2d_directcp  (nd, sources, ns, charge_loc, targ, nt,
     1                       pot_f, thresh)
         call c2d_directcp_c(nd, sources, ns, charge_loc, targ, nt,
     1                       pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directcp ', nd, ns, nt, errmax, nfail)

c        ---- c2d_directcg ----
         call zerocplx(pot_f,  nd*nt)
         call zerocplx(pot_c,  nd*nt)
         call zerocplx(grad_f, nd*nt)
         call zerocplx(grad_c, nd*nt)
         call c2d_directcg  (nd, sources, ns, charge_loc, targ, nt,
     1                       pot_f, grad_f, thresh)
         call c2d_directcg_c(nd, sources, ns, charge_loc, targ, nt,
     1                       pot_c, grad_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(grad_f(ii,j) - grad_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directcg ', nd, ns, nt, errmax, nfail)

c        ---- c2d_directch ----
         call zerocplx(pot_f,  nd*nt)
         call zerocplx(pot_c,  nd*nt)
         call zerocplx(grad_f, nd*nt)
         call zerocplx(grad_c, nd*nt)
         call zerocplx(hess_f, nd*nt)
         call zerocplx(hess_c, nd*nt)
         call c2d_directch  (nd, sources, ns, charge_loc, targ, nt,
     1                       pot_f, grad_f, hess_f, thresh)
         call c2d_directch_c(nd, sources, ns, charge_loc, targ, nt,
     1                       pot_c, grad_c, hess_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(grad_f(ii,j) - grad_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(hess_f(ii,j) - hess_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directch ', nd, ns, nt, errmax, nfail)

c        ---- c2d_directdp ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call c2d_directdp  (nd, sources, ns, dipstr_loc, targ, nt,
     1                       pot_f, thresh)
         call c2d_directdp_c(nd, sources, ns, dipstr_loc, targ, nt,
     1                       pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directdp ', nd, ns, nt, errmax, nfail)

c        ---- c2d_directdg ----
         call zerocplx(pot_f,  nd*nt)
         call zerocplx(pot_c,  nd*nt)
         call zerocplx(grad_f, nd*nt)
         call zerocplx(grad_c, nd*nt)
         call c2d_directdg  (nd, sources, ns, dipstr_loc, targ, nt,
     1                       pot_f, grad_f, thresh)
         call c2d_directdg_c(nd, sources, ns, dipstr_loc, targ, nt,
     1                       pot_c, grad_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(grad_f(ii,j) - grad_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directdg ', nd, ns, nt, errmax, nfail)

c        ---- c2d_directdh ----
         call zerocplx(pot_f,  nd*nt)
         call zerocplx(pot_c,  nd*nt)
         call zerocplx(grad_f, nd*nt)
         call zerocplx(grad_c, nd*nt)
         call zerocplx(hess_f, nd*nt)
         call zerocplx(hess_c, nd*nt)
         call c2d_directdh  (nd, sources, ns, dipstr_loc, targ, nt,
     1                       pot_f, grad_f, hess_f, thresh)
         call c2d_directdh_c(nd, sources, ns, dipstr_loc, targ, nt,
     1                       pot_c, grad_c, hess_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(grad_f(ii,j) - grad_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(hess_f(ii,j) - hess_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directdh ', nd, ns, nt, errmax, nfail)

c        ---- c2d_directcdp ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call c2d_directcdp  (nd, sources, ns, charge_loc, dipstr_loc,
     1                        targ, nt, pot_f, thresh)
         call c2d_directcdp_c(nd, sources, ns, charge_loc, dipstr_loc,
     1                        targ, nt, pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directcdp', nd, ns, nt, errmax, nfail)

c        ---- c2d_directcdg ----
         call zerocplx(pot_f,  nd*nt)
         call zerocplx(pot_c,  nd*nt)
         call zerocplx(grad_f, nd*nt)
         call zerocplx(grad_c, nd*nt)
         call c2d_directcdg  (nd, sources, ns, charge_loc, dipstr_loc,
     1                        targ, nt, pot_f, grad_f, thresh)
         call c2d_directcdg_c(nd, sources, ns, charge_loc, dipstr_loc,
     1                        targ, nt, pot_c, grad_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(grad_f(ii,j) - grad_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directcdg', nd, ns, nt, errmax, nfail)

c        ---- c2d_directcdh ----
         call zerocplx(pot_f,  nd*nt)
         call zerocplx(pot_c,  nd*nt)
         call zerocplx(grad_f, nd*nt)
         call zerocplx(grad_c, nd*nt)
         call zerocplx(hess_f, nd*nt)
         call zerocplx(hess_c, nd*nt)
         call c2d_directcdh  (nd, sources, ns, charge_loc, dipstr_loc,
     1                        targ, nt, pot_f, grad_f, hess_f, thresh)
         call c2d_directcdh_c(nd, sources, ns, charge_loc, dipstr_loc,
     1                        targ, nt, pot_c, grad_c, hess_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(grad_f(ii,j) - grad_c(ii,j))
               if (e .gt. errmax) errmax = e
               e = cdabs(hess_f(ii,j) - hess_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('c2d_directcdh', nd, ns, nt, errmax, nfail)

         deallocate(pot_f, pot_c, grad_f, grad_c, hess_f, hess_c)
         deallocate(charge_loc, dipstr_loc)

      enddo

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all cauchykernels2d cases match'

      end


c     ----------------------------------------------------------------
      subroutine zerocplx(z, n)
      implicit none
      integer n, i
      complex *16 z(n)
      do i = 1, n
         z(i) = dcmplx(0.0d0, 0.0d0)
      enddo
      return
      end


c     ----------------------------------------------------------------
      subroutine report(name, nd, ns, nt, errmax, nfail)
      implicit none
      character*(*) name
      integer nd, ns, nt, nfail
      real *8 errmax, tol
      tol = 1.0d-13
      if (errmax .lt. tol) then
         write(*,1000) name, nd, ns, nt, errmax
      else
         write(*,1001) name, nd, ns, nt, errmax
         nfail = nfail + 1
      endif
 1000 format(' [ ok ] ',a13,'  nd=',i1,' ns=',i3,' nt=',i3,
     1       '  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a13,'  nd=',i1,' ns=',i3,' nt=',i3,
     1       '  errmax=',1pe11.3)
      return
      end
