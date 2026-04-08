c     test_bhkernels2d.f - differential test for the C translation of
c     src/biharmonic/bhkernels2d.f.
c
c     For each kernel, allocates output buffers with leading dimension
c     equal to the actual nd (NOT ndmax) so that the routines see the
c     storage stride they expect. Calls the Fortran reference and the
c     _c_ translation on identical inputs and asserts that every
c     complex element agrees in absolute value.

      program test_bhkernels2d
      implicit none

      external bh2d_directcp,  bh2d_directcp_c
      external bh2d_directcg,  bh2d_directcg_c
      external bh2d_directdp,  bh2d_directdp_c
      external bh2d_directdg,  bh2d_directdg_c
      external bh2d_directcdp, bh2d_directcdp_c
      external bh2d_directcdg, bh2d_directcdg_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 50, ntmax = 40, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      complex *16 charge(ndmax, 2, nsmax)
      complex *16 dippar(ndmax, 3, nsmax)

c     Per-test allocatable buffers, sized per iteration
      complex *16, allocatable :: vel_f(:,:),  vel_c(:,:)
      complex *16, allocatable :: grad_f(:,:,:), grad_c(:,:,:)
c     Local copies of charges / dippar resized for the call
      complex *16, allocatable :: charge_loc(:,:,:)
      complex *16, allocatable :: dippar_loc(:,:,:)

      integer ns, nt
      integer nds(2)
      data nds / 1, 3 /

      integer ind, nd, i, j, ii, kk
      integer nfail
      real *8 thresh, dummy, errmax, e

      ns = nsmax
      nt = ntmax
      thresh = 1.0d-14
      nfail = 0

c     seed hkrand once
      dummy = hkrand(1234)

c     fill sources in the unit square
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
      enddo
c     targets scaled by 2 so they lie away from sources
      do j = 1, nt
         targ(1, j) = 2.0d0 * hkrand(0)
         targ(2, j) = 2.0d0 * hkrand(0)
      enddo
c     fill charge (nd, 2, ns) and dippar (nd, 3, ns) at maximum sizes
      do i = 1, ns
         do kk = 1, 2
            do ii = 1, ndmax
               charge(ii, kk, i) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                                    2.0d0*hkrand(0)-1.0d0)
            enddo
         enddo
         do kk = 1, 3
            do ii = 1, ndmax
               dippar(ii, kk, i) = dcmplx(2.0d0*hkrand(0)-1.0d0,
     1                                    2.0d0*hkrand(0)-1.0d0)
            enddo
         enddo
      enddo

c     ============================================================
c     Loop over the two nd values and run all 6 kernels for each
c     ============================================================
      do ind = 1, 2
         nd = nds(ind)

c        Resize charge / dippar to the actual nd so the routine's
c        view of their storage matches their declared shape exactly.
         allocate(charge_loc(nd, 2, ns))
         allocate(dippar_loc(nd, 3, ns))
         do i = 1, ns
            do kk = 1, 2
               do ii = 1, nd
                  charge_loc(ii, kk, i) = charge(ii, kk, i)
               enddo
            enddo
            do kk = 1, 3
               do ii = 1, nd
                  dippar_loc(ii, kk, i) = dippar(ii, kk, i)
               enddo
            enddo
         enddo

c        Allocate output buffers at the correct shape.
         allocate(vel_f(nd, nt),  vel_c(nd, nt))
         allocate(grad_f(nd, 3, nt), grad_c(nd, 3, nt))

c        ---- bh2d_directcp ----
         call zerocplx(vel_f, nd*nt)
         call zerocplx(vel_c, nd*nt)
         call bh2d_directcp  (nd, sources, ns, charge_loc, targ, nt,
     1                        vel_f, thresh)
         call bh2d_directcp_c(nd, sources, ns, charge_loc, targ, nt,
     1                        vel_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(vel_f(ii,j) - vel_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('bh2d_directcp ', nd, ns, nt, errmax, nfail)

c        ---- bh2d_directcg ----
         call zerocplx(vel_f,  nd*nt)
         call zerocplx(vel_c,  nd*nt)
         call zerocplx(grad_f, nd*3*nt)
         call zerocplx(grad_c, nd*3*nt)
         call bh2d_directcg  (nd, sources, ns, charge_loc, targ, nt,
     1                        vel_f, grad_f, thresh)
         call bh2d_directcg_c(nd, sources, ns, charge_loc, targ, nt,
     1                        vel_c, grad_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(vel_f(ii,j) - vel_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
            do kk = 1, 3
               do ii = 1, nd
                  e = cdabs(grad_f(ii,kk,j) - grad_c(ii,kk,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('bh2d_directcg ', nd, ns, nt, errmax, nfail)

c        ---- bh2d_directdp ----
         call zerocplx(vel_f, nd*nt)
         call zerocplx(vel_c, nd*nt)
         call bh2d_directdp  (nd, sources, ns, dippar_loc, targ, nt,
     1                        vel_f, thresh)
         call bh2d_directdp_c(nd, sources, ns, dippar_loc, targ, nt,
     1                        vel_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(vel_f(ii,j) - vel_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('bh2d_directdp ', nd, ns, nt, errmax, nfail)

c        ---- bh2d_directdg ----
         call zerocplx(vel_f,  nd*nt)
         call zerocplx(vel_c,  nd*nt)
         call zerocplx(grad_f, nd*3*nt)
         call zerocplx(grad_c, nd*3*nt)
         call bh2d_directdg  (nd, sources, ns, dippar_loc, targ, nt,
     1                        vel_f, grad_f, thresh)
         call bh2d_directdg_c(nd, sources, ns, dippar_loc, targ, nt,
     1                        vel_c, grad_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(vel_f(ii,j) - vel_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
            do kk = 1, 3
               do ii = 1, nd
                  e = cdabs(grad_f(ii,kk,j) - grad_c(ii,kk,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('bh2d_directdg ', nd, ns, nt, errmax, nfail)

c        ---- bh2d_directcdp ----
         call zerocplx(vel_f, nd*nt)
         call zerocplx(vel_c, nd*nt)
         call bh2d_directcdp  (nd, sources, ns, charge_loc, dippar_loc,
     1                         targ, nt, vel_f, thresh)
         call bh2d_directcdp_c(nd, sources, ns, charge_loc, dippar_loc,
     1                         targ, nt, vel_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(vel_f(ii,j) - vel_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('bh2d_directcdp', nd, ns, nt, errmax, nfail)

c        ---- bh2d_directcdg ----
         call zerocplx(vel_f,  nd*nt)
         call zerocplx(vel_c,  nd*nt)
         call zerocplx(grad_f, nd*3*nt)
         call zerocplx(grad_c, nd*3*nt)
         call bh2d_directcdg  (nd, sources, ns, charge_loc, dippar_loc,
     1                         targ, nt, vel_f, grad_f, thresh)
         call bh2d_directcdg_c(nd, sources, ns, charge_loc, dippar_loc,
     1                         targ, nt, vel_c, grad_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(vel_f(ii,j) - vel_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
            do kk = 1, 3
               do ii = 1, nd
                  e = cdabs(grad_f(ii,kk,j) - grad_c(ii,kk,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('bh2d_directcdg', nd, ns, nt, errmax, nfail)

         deallocate(vel_f, vel_c, grad_f, grad_c)
         deallocate(charge_loc, dippar_loc)

      enddo

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all bhkernels2d cases match'

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
 1000 format(' [ ok ] ',a14,'  nd=',i1,' ns=',i3,' nt=',i3,
     1       '  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a14,'  nd=',i1,' ns=',i3,' nt=',i3,
     1       '  errmax=',1pe11.3)
      return
      end
