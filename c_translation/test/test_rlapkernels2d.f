c     test_rlapkernels2d.f - differential test for the C translation of
c     src/laplace/rlapkernels2d.f.

      program test_rlapkernels2d
      implicit none

      external r2d_directcp,  r2d_directcp_c
      external r2d_directcg,  r2d_directcg_c
      external r2d_directch,  r2d_directch_c
      external r2d_directdp,  r2d_directdp_c
      external r2d_directdg,  r2d_directdg_c
      external r2d_directdh,  r2d_directdh_c
      external r2d_directcdp, r2d_directcdp_c
      external r2d_directcdg, r2d_directcdg_c
      external r2d_directcdh, r2d_directcdh_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 50, ntmax = 40, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      real *8 charge(ndmax, nsmax)
      real *8 dipstr(ndmax, nsmax)
      real *8 dipvec(ndmax, 2, nsmax)

c     Per-test allocatable buffers
      real *8, allocatable :: pot_f(:,:)
      real *8, allocatable :: pot_c(:,:)
      real *8, allocatable :: grd_f(:,:,:)
      real *8, allocatable :: grd_c(:,:,:)
      real *8, allocatable :: hes_f(:,:,:)
      real *8, allocatable :: hes_c(:,:,:)
      real *8, allocatable :: chloc(:,:)
      real *8, allocatable :: dsloc(:,:)
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
            charge(ii, i) = hkrand(0)
            dipstr(ii, i) = hkrand(0)
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

c        ---- r2d_directcp ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call r2d_directcp(nd, sources, ns,
     1      chloc, targ, nt, pot_f, thresh)
         call r2d_directcp_c(nd, sources, ns,
     1      chloc, targ, nt, pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('r2d_directcp ', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directcg ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grd_f, nd*2*nt)
         call zeroreal(grd_c, nd*2*nt)
         call r2d_directcg(nd, sources, ns,
     1      chloc, targ, nt, pot_f, grd_f,
     2      thresh)
         call r2d_directcg_c(nd, sources, ns,
     1      chloc, targ, nt, pot_c, grd_c,
     2      thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = dabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('r2d_directcg ', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directch ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grd_f, nd*2*nt)
         call zeroreal(grd_c, nd*2*nt)
         call zeroreal(hes_f, nd*3*nt)
         call zeroreal(hes_c, nd*3*nt)
         call r2d_directch(nd, sources, ns,
     1      chloc, targ, nt, pot_f, grd_f,
     2      hes_f, thresh)
         call r2d_directch_c(nd, sources, ns,
     1      chloc, targ, nt, pot_c, grd_c,
     2      hes_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = dabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         do j = 1, nt
            do k = 1, 3
               do ii = 1, nd
                  e = dabs(hes_f(ii,k,j)
     1                - hes_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('r2d_directch ', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directdp ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call r2d_directdp(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_f,
     2      thresh)
         call r2d_directdp_c(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_c,
     2      thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('r2d_directdp ', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directdg ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grd_f, nd*2*nt)
         call zeroreal(grd_c, nd*2*nt)
         call r2d_directdg(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_f,
     2      grd_f, thresh)
         call r2d_directdg_c(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_c,
     2      grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = dabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('r2d_directdg ', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directdh ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grd_f, nd*2*nt)
         call zeroreal(grd_c, nd*2*nt)
         call zeroreal(hes_f, nd*3*nt)
         call zeroreal(hes_c, nd*3*nt)
         call r2d_directdh(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_f,
     2      grd_f, hes_f, thresh)
         call r2d_directdh_c(nd, sources, ns,
     1      dsloc, dvloc, targ, nt, pot_c,
     2      grd_c, hes_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = dabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         do j = 1, nt
            do k = 1, 3
               do ii = 1, nd
                  e = dabs(hes_f(ii,k,j)
     1                - hes_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('r2d_directdh ', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directcdp ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call r2d_directcdp(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_f, thresh)
         call r2d_directcdp_c(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         call report('r2d_directcdp', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directcdg ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grd_f, nd*2*nt)
         call zeroreal(grd_c, nd*2*nt)
         call r2d_directcdg(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_f, grd_f, thresh)
         call r2d_directcdg_c(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_c, grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = dabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('r2d_directcdg', nd, ns,
     1      nt, errmax, nfail)

c        ---- r2d_directcdh ----
         call zeroreal(pot_f, nd*nt)
         call zeroreal(pot_c, nd*nt)
         call zeroreal(grd_f, nd*2*nt)
         call zeroreal(grd_c, nd*2*nt)
         call zeroreal(hes_f, nd*3*nt)
         call zeroreal(hes_c, nd*3*nt)
         call r2d_directcdh(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_f, grd_f, hes_f, thresh)
         call r2d_directcdh_c(nd, sources, ns,
     1      chloc, dsloc, dvloc, targ, nt,
     2      pot_c, grd_c, hes_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do ii = 1, nd
               e = dabs(pot_f(ii,j) - pot_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
         enddo
         do j = 1, nt
            do k = 1, 2
               do ii = 1, nd
                  e = dabs(grd_f(ii,k,j)
     1                - grd_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         do j = 1, nt
            do k = 1, 3
               do ii = 1, nd
                  e = dabs(hes_f(ii,k,j)
     1                - hes_c(ii,k,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
         enddo
         call report('r2d_directcdh', nd, ns,
     1      nt, errmax, nfail)

         deallocate(pot_f, pot_c)
         deallocate(grd_f, grd_c)
         deallocate(hes_f, hes_c)
         deallocate(chloc, dsloc, dvloc)
      enddo

      if (nfail .gt. 0) then
         write(*,*) 'FAIL: ',nfail,
     1      ' rlapkernels2d tests failed'
         stop 1
      endif
      write(*,*) 'PASS: all rlapkernels2d cases match'

      end


c     ----------------------------------------------------------------
      subroutine zeroreal(a, n)
      implicit none
      integer n, i
      real *8 a(n)
      do i = 1, n
         a(i) = 0.0d0
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
