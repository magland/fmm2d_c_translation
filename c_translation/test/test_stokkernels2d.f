c     test_stokkernels2d.f - differential test for the C translation of
c     src/stokes/stokkernels2d.f.

      program test_stokkernels2d
      implicit none

      external st2ddirectstokg
      external st2ddirectstokg_c
      external st2ddirectstokstrsg
      external st2ddirectstokstrsg_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax, ndmax
      parameter (nsmax = 50, ntmax = 40, ndmax = 3)

      real *8 sources(2, nsmax), targ(2, ntmax)
      real *8 stoklet(ndmax, 2, nsmax)
      real *8 strslet(ndmax, 2, nsmax)
      real *8 strsvec(ndmax, 2, nsmax)

c     Per-test allocatable buffers
      real *8, allocatable :: pot_f(:,:,:)
      real *8, allocatable :: pot_c(:,:,:)
      real *8, allocatable :: pre_f(:,:)
      real *8, allocatable :: pre_c(:,:)
      real *8, allocatable :: grd_f(:,:,:,:)
      real *8, allocatable :: grd_c(:,:,:,:)
      real *8, allocatable :: stkloc(:,:,:)
      real *8, allocatable :: strloc(:,:,:)
      real *8, allocatable :: svloc(:,:,:)

      integer ns, nt
      integer nds(2)
      data nds / 1, 3 /

      integer ind, nd, i, j, ii
      integer k1, k2
      integer nfail, ifstok, istrs
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
c     fill stoklet, strslet, strsvec
      do i = 1, ns
         do ii = 1, ndmax
            stoklet(ii, 1, i) = hkrand(0)
            stoklet(ii, 2, i) = hkrand(0)
            strslet(ii, 1, i) = hkrand(0)
            strslet(ii, 2, i) = hkrand(0)
            strsvec(ii, 1, i) = hkrand(0)
            strsvec(ii, 2, i) = hkrand(0)
         enddo
      enddo

c     ============================================================
c     Loop over nd values
c     ============================================================
      do ind = 1, 2
         nd = nds(ind)

         allocate(stkloc(nd, 2, ns))
         allocate(strloc(nd, 2, ns))
         allocate(svloc(nd, 2, ns))
         do i = 1, ns
            do ii = 1, nd
               stkloc(ii, 1, i) = stoklet(ii, 1, i)
               stkloc(ii, 2, i) = stoklet(ii, 2, i)
               strloc(ii, 1, i) = strslet(ii, 1, i)
               strloc(ii, 2, i) = strslet(ii, 2, i)
               svloc(ii, 1, i) = strsvec(ii, 1, i)
               svloc(ii, 2, i) = strsvec(ii, 2, i)
            enddo
         enddo

         allocate(pot_f(nd, 2, nt))
         allocate(pot_c(nd, 2, nt))
         allocate(pre_f(nd, nt))
         allocate(pre_c(nd, nt))
         allocate(grd_f(nd, 2, 2, nt))
         allocate(grd_c(nd, 2, 2, nt))

c        ---- st2ddirectstokg ----
         call zeroreal(pot_f, nd*2*nt)
         call zeroreal(pot_c, nd*2*nt)
         call zeroreal(pre_f, nd*nt)
         call zeroreal(pre_c, nd*nt)
         call zeroreal(grd_f, nd*2*2*nt)
         call zeroreal(grd_c, nd*2*2*nt)
         call st2ddirectstokg(nd, sources,
     1      stkloc, ns, targ, nt, pot_f,
     2      pre_f, grd_f, thresh)
         call st2ddirectstokg_c(nd, sources,
     1      stkloc, ns, targ, nt, pot_c,
     2      pre_c, grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do k1 = 1, 2
               do ii = 1, nd
                  e = dabs(pot_f(ii,k1,j)
     1                - pot_c(ii,k1,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
            do ii = 1, nd
               e = dabs(pre_f(ii,j) - pre_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
            do k2 = 1, 2
               do k1 = 1, 2
                  do ii = 1, nd
                     e = dabs(grd_f(ii,k1,k2,j)
     1                   - grd_c(ii,k1,k2,j))
                     if (e .gt. errmax)
     1                  errmax = e
                  enddo
               enddo
            enddo
         enddo
         call report('st2dstokg    ', nd, ns,
     1      nt, errmax, nfail)

c        ---- st2ddirectstokstrsg (stoklet + stresslet) ----
         call zeroreal(pot_f, nd*2*nt)
         call zeroreal(pot_c, nd*2*nt)
         call zeroreal(pre_f, nd*nt)
         call zeroreal(pre_c, nd*nt)
         call zeroreal(grd_f, nd*2*2*nt)
         call zeroreal(grd_c, nd*2*2*nt)
         ifstok = 1
         istrs = 1
         call st2ddirectstokstrsg(nd, sources,
     1      ifstok, stkloc, istrs, strloc,
     2      svloc, ns, targ, nt, pot_f,
     3      pre_f, grd_f, thresh)
         call st2ddirectstokstrsg_c(nd,
     1      sources, ifstok, stkloc, istrs,
     2      strloc, svloc, ns, targ, nt,
     3      pot_c, pre_c, grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do k1 = 1, 2
               do ii = 1, nd
                  e = dabs(pot_f(ii,k1,j)
     1                - pot_c(ii,k1,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
            do ii = 1, nd
               e = dabs(pre_f(ii,j) - pre_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
            do k2 = 1, 2
               do k1 = 1, 2
                  do ii = 1, nd
                     e = dabs(grd_f(ii,k1,k2,j)
     1                   - grd_c(ii,k1,k2,j))
                     if (e .gt. errmax)
     1                  errmax = e
                  enddo
               enddo
            enddo
         enddo
         call report('st2dstokstrsg', nd, ns,
     1      nt, errmax, nfail)

c        ---- st2ddirectstokstrsg (stresslet only) ----
         call zeroreal(pot_f, nd*2*nt)
         call zeroreal(pot_c, nd*2*nt)
         call zeroreal(pre_f, nd*nt)
         call zeroreal(pre_c, nd*nt)
         call zeroreal(grd_f, nd*2*2*nt)
         call zeroreal(grd_c, nd*2*2*nt)
         ifstok = 0
         istrs = 1
         call st2ddirectstokstrsg(nd, sources,
     1      ifstok, stkloc, istrs, strloc,
     2      svloc, ns, targ, nt, pot_f,
     3      pre_f, grd_f, thresh)
         call st2ddirectstokstrsg_c(nd,
     1      sources, ifstok, stkloc, istrs,
     2      strloc, svloc, ns, targ, nt,
     3      pot_c, pre_c, grd_c, thresh)
         errmax = 0.0d0
         do j = 1, nt
            do k1 = 1, 2
               do ii = 1, nd
                  e = dabs(pot_f(ii,k1,j)
     1                - pot_c(ii,k1,j))
                  if (e .gt. errmax) errmax = e
               enddo
            enddo
            do ii = 1, nd
               e = dabs(pre_f(ii,j) - pre_c(ii,j))
               if (e .gt. errmax) errmax = e
            enddo
            do k2 = 1, 2
               do k1 = 1, 2
                  do ii = 1, nd
                     e = dabs(grd_f(ii,k1,k2,j)
     1                   - grd_c(ii,k1,k2,j))
                     if (e .gt. errmax)
     1                  errmax = e
                  enddo
               enddo
            enddo
         enddo
         call report('st2dstrs_only', nd, ns,
     1      nt, errmax, nfail)

         deallocate(pot_f, pot_c)
         deallocate(pre_f, pre_c)
         deallocate(grd_f, grd_c)
         deallocate(stkloc, strloc, svloc)
      enddo

      if (nfail .gt. 0) then
         write(*,*) 'FAIL: ',nfail,
     1      ' stokkernels2d tests failed'
         stop 1
      endif
      write(*,*) 'PASS: all stokkernels2d cases match'

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
