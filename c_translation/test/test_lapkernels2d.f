c     test_lapkernels2d.f - differential test for the C translation of
c     src/laplace/lapkernels2d.f.
c
c     For each kernel, allocates output buffers with leading dimension
c     equal to the actual nd. Calls the Fortran reference and the _c_
c     translation on identical inputs and asserts element-by-element
c     agreement.

      program test_lapkernels2d
      implicit none

      external l2d_directcp,  l2d_directcp_c
      external l2d_directcg,  l2d_directcg_c
      external l2d_directch,  l2d_directch_c
      external l2d_directdp,  l2d_directdp_c
      external l2d_directdg,  l2d_directdg_c
      external l2d_directdh,  l2d_directdh_c
      external l2d_directcdp, l2d_directcdp_c
      external l2d_directcdg, l2d_directcdg_c
      external l2d_directcdh, l2d_directcdh_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 50, ntmax = 40, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      complex *16 charge(ndmax, nsmax)
      complex *16 dipstr(ndmax, nsmax)
      real *8 dipvec(ndmax, 2, nsmax)

c     Per-test allocatable buffers
      complex *16, allocatable :: pot_f(:,:)
      complex *16, allocatable :: pot_c(:,:)
      complex *16, allocatable :: grd_f(:,:,:)
      complex *16, allocatable :: grd_c(:,:,:)
      complex *16, allocatable :: hes_f(:,:,:)
      complex *16, allocatable :: hes_c(:,:,:)
      complex *16, allocatable :: chloc(:,:)
      complex *16, allocatable :: dsloc(:,:)
      real *8, allocatable :: dvloc(:,:,:)

      integer ns, nt
      integer nds(2)
      data nds / 1, 3 /

      integer ind, nd, i, j, ii, k
      integer nfail
      real *8 thresh, dummy, errmax, e

      ns = nsmax
      nt = ntmax
      thresh = 1.0d-15
      nfail = 0

c     seed hkrand once
      dummy = hkrand(1234)

c     fill sources, targets
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
      enddo
      do j = 1, nt
         targ(1, j) = hkrand(0)
         targ(2, j) = hkrand(0)
      enddo
c     fill charge, dipstr, dipvec at maximum sizes
      do i = 1, ns
         do ii = 1, ndmax
            charge(ii, i) = dcmplx(hkrand(0), hkrand(0))
            dipstr(ii, i) = dcmplx(hkrand(0), hkrand(0))
            dipvec(ii, 1, i) = hkrand(0)
            dipvec(ii, 2, i) = hkrand(0)
         enddo
      enddo

c     ============================================================
c     Loop over nd values
c     ============================================================
      do ind = 1, 2
         nd = nds(ind)

c        Copy inputs to (nd, ns) sized arrays
         allocate(chloc(nd, ns))
         allocate(dsloc(nd, ns))
         allocate(dvloc(nd, 2, ns))
         do i = 1, ns
            do ii = 1, nd
               chloc(ii, i) = charge(ii, i)
               dsloc(ii, i) = dipstr(ii, i)
               dvloc(ii, 1, i) = dipvec(ii, 1, i)
               dvloc(ii, 2, i) = dipvec(ii, 2, i)
            enddo
         enddo

         allocate(pot_f(nd, nt), pot_c(nd, nt))
         allocate(grd_f(nd, 2, nt))
         allocate(grd_c(nd, 2, nt))
         allocate(hes_f(nd, 3, nt))
         allocate(hes_c(nd, 3, nt))

c        ---- l2d_directcp ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call l2d_directcp(nd, sources, ns,
     1      chloc, targ, nt, pot_f, thresh)
         call l2d_directcp_c(nd, sources, ns,
     1      chloc, targ, nt, pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('l2d_directcp ', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directcg ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call zerocplx(grd_f, nd*2*nt)
         call zerocplx(grd_c, nd*2*nt)
         call l2d_directcg(nd, sources, ns,
     1      chloc, targ, nt, pot_f, grd_f, thresh)
         call l2d_directcg_c(nd, sources, ns,
     1      chloc, targ, nt, pot_c, grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = cdabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('l2d_directcg ', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directch ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call zerocplx(grd_f, nd*2*nt)
         call zerocplx(grd_c, nd*2*nt)
         call zerocplx(hes_f, nd*3*nt)
         call zerocplx(hes_c, nd*3*nt)
         call l2d_directch(nd, sources, ns,
     1      chloc, targ, nt, pot_f, grd_f,
     2      hes_f, thresh)
         call l2d_directch_c(nd, sources, ns,
     1      chloc, targ, nt, pot_c, grd_c,
     2      hes_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = cdabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         do j = 1, nt
            do k = 1, 3
               do ii = 1, nd
                  e = cdabs(hes_f(ii,k,j)
     1                - hes_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('l2d_directch ', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directdp ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call l2d_directdp(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_f,
     2      thresh)
         call l2d_directdp_c(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_c,
     2      thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('l2d_directdp ', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directdg ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call zerocplx(grd_f, nd*2*nt)
         call zerocplx(grd_c, nd*2*nt)
         call l2d_directdg(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_f,
     2      grd_f, thresh)
         call l2d_directdg_c(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_c,
     2      grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = cdabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('l2d_directdg ', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directdh ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call zerocplx(grd_f, nd*2*nt)
         call zerocplx(grd_c, nd*2*nt)
         call zerocplx(hes_f, nd*3*nt)
         call zerocplx(hes_c, nd*3*nt)
         call l2d_directdh(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_f,
     2      grd_f, hes_f, thresh)
         call l2d_directdh_c(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_c,
     2      grd_c, hes_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = cdabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         do j = 1, nt
            do k = 1, 3
               do ii = 1, nd
                  e = cdabs(hes_f(ii,k,j)
     1                - hes_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('l2d_directdh ', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directcdp ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call l2d_directcdp(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_f, thresh)
         call l2d_directcdp_c(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('l2d_directcdp', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directcdg ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call zerocplx(grd_f, nd*2*nt)
         call zerocplx(grd_c, nd*2*nt)
         call l2d_directcdg(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_f, grd_f, thresh)
         call l2d_directcdg_c(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_c, grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = cdabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('l2d_directcdg', nd, ns,
     1      nt, errmax, nfail)

c        ---- l2d_directcdh ----
         call zerocplx(pot_f, nd*nt)
         call zerocplx(pot_c, nd*nt)
         call zerocplx(grd_f, nd*2*nt)
         call zerocplx(grd_c, nd*2*nt)
         call zerocplx(hes_f, nd*3*nt)
         call zerocplx(hes_c, nd*3*nt)
         call l2d_directcdh(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_f, grd_f, hes_f, thresh)
         call l2d_directcdh_c(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_c, grd_c, hes_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = cdabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = cdabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         do j = 1, nt
            do k = 1, 3
               do ii = 1, nd
                  e = cdabs(hes_f(ii,k,j)
     1                - hes_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('l2d_directcdh', nd, ns,
     1      nt, errmax, nfail)

         deallocate(pot_f, pot_c)
         deallocate(grd_f, grd_c)
         deallocate(hes_f, hes_c)
         deallocate(chloc, dsloc, dvloc)
      enddo

      if (nfail .gt. 0) then
         write(*,*) 'FAIL: ',nfail,
     1      ' lapkernels2d tests failed'
         stop 1
      endif
      write(*,*) 'PASS: all lapkernels2d cases match'

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
      subroutine report(name, nd, ns, nt,
     1   errmax, nfail)
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
 1000 format(' [ ok ] ',a13,'  nd=',i1,
     1   ' ns=',i3,' nt=',i3,
     2   '  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a13,'  nd=',i1,
     1   ' ns=',i3,' nt=',i3,
     2   '  errmax=',1pe11.3)
      return
      end
