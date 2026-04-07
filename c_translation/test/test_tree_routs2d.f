c     test_tree_routs2d.f - differential test for the C translation of
c     src/common/tree_routs2d.f.
c
c     For each translated routine, calls both the Fortran reference
c     and the _c_ C port with identical inputs into separate output
c     buffers and asserts exact integer equality (these routines do
c     bookkeeping, not floating-point arithmetic, so the tolerance
c     is zero).
c
c     A real quad tree is built via pts_tree_mem / pts_tree_build to
c     feed realistic inputs to computecoll, computelists, updateflags,
c     computemnlists, and tree_copy. The two tree_refine_boxes_*
c     routines are exercised with small synthetic inputs.

      program test_tree_routs2d
      implicit none

      external pts_tree_mem, pts_tree_build
      external tree_refine_boxes, tree_refine_boxes_c
      external tree_copy, tree_copy_c
      external computecoll, computecoll_c
      external updateflags, updateflags_c
      external tree_refine_boxes_flag, tree_refine_boxes_flag_c
      external computemnlists, computemnlists_c
      external computelists, computelists_c
      external hkrand
      real *8 hkrand

      integer ns, nt
      parameter (ns = 200, nt = 0)
      real *8 sources(2, ns)
      real *8 targ(2, 1)
      integer ndiv, idivflag, nlmin, nlmax, ifunif, iper
      integer nlevels, nboxes, ltree
      integer, allocatable :: itree(:)
      integer iptr(8)
      real *8, allocatable :: tcenters(:, :), boxsize(:)

c     Local unpacked tree arrays
      integer, allocatable :: laddr(:)
      integer, allocatable :: iparent(:), nchild(:), ichild(:)

c     Colleague buffers
      integer, allocatable :: nnbors_f(:), nbors_f(:)
      integer, allocatable :: nnbors_c(:), nbors_c(:)

c     List buffers
      integer mnlist1_f, mnlist2_f, mnlist3_f, mnlist4_f
      integer mnlist1_c, mnlist2_c, mnlist3_c, mnlist4_c
      integer mn1, mn2, mn3, mn4
      integer, allocatable :: nlist1_f(:), list1_f(:)
      integer, allocatable :: nlist2_f(:), list2_f(:)
      integer, allocatable :: nlist3_f(:), list3_f(:)
      integer, allocatable :: nlist4_f(:), list4_f(:)
      integer, allocatable :: nlist1_c(:), list1_c(:)
      integer, allocatable :: nlist2_c(:), list2_c(:)
      integer, allocatable :: nlist3_c(:), list3_c(:)
      integer, allocatable :: nlist4_c(:), list4_c(:)

c     updateflags buffers
      integer, allocatable :: iflag_f(:), iflag_c(:)
      integer curlev

c     tree_refine_boxes buffers (synthetic small input)
      integer nb_s_p
      parameter (nb_s_p = 50)
      integer nb_s
      real *8 centers_f(2, nb_s_p), centers_c(2, nb_s_p)
      integer ilevel_f(nb_s_p), ilevel_c(nb_s_p)
      integer iparent_f(nb_s_p), iparent_c(nb_s_p)
      integer nchild_f(nb_s_p), nchild_c(nb_s_p)
      integer ichild_f(4, nb_s_p), ichild_c(4, nb_s_p)
      integer irefinebox(4)
      integer iflag_s_f(nb_s_p), iflag_s_c(nb_s_p)
      integer nbctr_f, nbctr_c
      integer nlctr, ifirstbox, nbloc
      real *8 bs

c     tree_copy buffers (reuse sizes via allocation)
      real *8, allocatable :: centers2_f(:, :), centers2_c(:, :)
      integer, allocatable :: ilevel2_f(:), ilevel2_c(:)
      integer, allocatable :: iparent2_f(:), iparent2_c(:)
      integer, allocatable :: nchild2_f(:), nchild2_c(:)
      integer, allocatable :: ichild2_f(:), ichild2_c(:)
      integer, allocatable :: ilevel_tree(:)

      real *8 dummy
      integer i, j, errcount, nfail

      nfail = 0
      nb_s = nb_s_p

      dummy = hkrand(1234)
      do i = 1, ns
         sources(1, i) = hkrand(0)
         sources(2, i) = hkrand(0)
      enddo

      ndiv = 8
      idivflag = 0
      nlmin = 0
      nlmax = 51
      ifunif = 0
      iper = 0
      ltree = 0
      nlevels = 0
      nboxes = 0
      call pts_tree_mem(sources, ns, targ, nt, idivflag, ndiv,
     1     nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree)

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(2, nboxes))
      call pts_tree_build(sources, ns, targ, nt, idivflag, ndiv,
     1     nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree,
     2     itree, iptr, tcenters, boxsize)

      write(*,*) 'Built tree: nlevels=', nlevels, ' nboxes=', nboxes

c     -----------------------------------------------------------------
c     Extract laddr, ilevel, iparent, nchild, ichild slices from itree
c     for convenient reuse below.
c     -----------------------------------------------------------------
      allocate(laddr(2*(nlevels+1)))
      allocate(ilevel_tree(nboxes))
      allocate(iparent(nboxes))
      allocate(nchild(nboxes))
      allocate(ichild(4*nboxes))

      do i = 1, 2*(nlevels+1)
         laddr(i) = itree(iptr(1) + i - 1)
      enddo
      do i = 1, nboxes
         ilevel_tree(i) = itree(iptr(2) + i - 1)
         iparent(i)     = itree(iptr(3) + i - 1)
         nchild(i)      = itree(iptr(4) + i - 1)
      enddo
      do i = 1, 4*nboxes
         ichild(i) = itree(iptr(5) + i - 1)
      enddo

c     =================================================================
c     computemnlists - trivial (just 4 constants).
c     =================================================================
      mnlist1_f = -1
      mnlist2_f = -1
      mnlist3_f = -1
      mnlist4_f = -1
      mnlist1_c = -1
      mnlist2_c = -1
      mnlist3_c = -1
      mnlist4_c = -1
      call computemnlists(nlevels, nboxes, itree, ltree,
     1     iptr, tcenters, boxsize, iper,
     2     mnlist1_f, mnlist2_f, mnlist3_f, mnlist4_f)
      call computemnlists_c(nlevels, nboxes, itree, ltree,
     1     iptr, tcenters, boxsize, iper,
     2     mnlist1_c, mnlist2_c, mnlist3_c, mnlist4_c)
      errcount = 0
      if (mnlist1_f .ne. mnlist1_c) errcount = errcount + 1
      if (mnlist2_f .ne. mnlist2_c) errcount = errcount + 1
      if (mnlist3_f .ne. mnlist3_c) errcount = errcount + 1
      if (mnlist4_f .ne. mnlist4_c) errcount = errcount + 1
      call report('computemnlists     ', errcount, nfail)

c     =================================================================
c     computecoll
c     =================================================================
      allocate(nnbors_f(nboxes), nbors_f(9*nboxes))
      allocate(nnbors_c(nboxes), nbors_c(9*nboxes))
      do i = 1, nboxes
         nnbors_f(i) = 0
         nnbors_c(i) = 0
      enddo
      do i = 1, 9*nboxes
         nbors_f(i) = 0
         nbors_c(i) = 0
      enddo
      call computecoll(nlevels, nboxes, laddr, boxsize, tcenters,
     1     iparent, nchild, ichild, iper, nnbors_f, nbors_f)
      call computecoll_c(nlevels, nboxes, laddr, boxsize, tcenters,
     1     iparent, nchild, ichild, iper, nnbors_c, nbors_c)
      errcount = 0
      do i = 1, nboxes
         if (nnbors_f(i) .ne. nnbors_c(i)) errcount = errcount + 1
      enddo
      do i = 1, 9*nboxes
         if (nbors_f(i) .ne. nbors_c(i)) errcount = errcount + 1
      enddo
      call report('computecoll        ', errcount, nfail)

c     =================================================================
c     computelists - use the Fortran colleague data (nnbors/nbors) and
c     feed it through itree by itself via the existing itree packing
c     (computelists reads nnbors/nbors from inside itree at iptr(6)
c     and iptr(7)). The tree built by pts_tree_build already has
c     these populated, so we just call both versions directly.
c     =================================================================
      mn1 = mnlist1_f
      mn2 = mnlist2_f
      mn3 = mnlist3_f
      mn4 = mnlist4_f

      allocate(nlist1_f(nboxes), list1_f(mn1*nboxes))
      allocate(nlist2_f(nboxes), list2_f(mn2*nboxes))
      allocate(nlist3_f(nboxes), list3_f(mn3*nboxes))
      allocate(nlist4_f(nboxes), list4_f(mn4*nboxes))
      allocate(nlist1_c(nboxes), list1_c(mn1*nboxes))
      allocate(nlist2_c(nboxes), list2_c(mn2*nboxes))
      allocate(nlist3_c(nboxes), list3_c(mn3*nboxes))
      allocate(nlist4_c(nboxes), list4_c(mn4*nboxes))
      do i = 1, nboxes
         nlist1_f(i) = 0
         nlist2_f(i) = 0
         nlist3_f(i) = 0
         nlist4_f(i) = 0
         nlist1_c(i) = 0
         nlist2_c(i) = 0
         nlist3_c(i) = 0
         nlist4_c(i) = 0
      enddo
      do i = 1, mn1*nboxes
         list1_f(i) = 0
         list1_c(i) = 0
      enddo
      do i = 1, mn2*nboxes
         list2_f(i) = 0
         list2_c(i) = 0
      enddo
      do i = 1, mn3*nboxes
         list3_f(i) = 0
         list3_c(i) = 0
      enddo
      do i = 1, mn4*nboxes
         list4_f(i) = 0
         list4_c(i) = 0
      enddo

      call computelists(nlevels, nboxes, itree, ltree, iptr,
     1     tcenters, boxsize, iper,
     2     mn1, nlist1_f, list1_f,
     3     mn2, nlist2_f, list2_f,
     4     mn3, nlist3_f, list3_f,
     5     mn4, nlist4_f, list4_f)
      call computelists_c(nlevels, nboxes, itree, ltree, iptr,
     1     tcenters, boxsize, iper,
     2     mn1, nlist1_c, list1_c,
     3     mn2, nlist2_c, list2_c,
     4     mn3, nlist3_c, list3_c,
     5     mn4, nlist4_c, list4_c)
      errcount = 0
      do i = 1, nboxes
         if (nlist1_f(i) .ne. nlist1_c(i)) errcount = errcount + 1
         if (nlist2_f(i) .ne. nlist2_c(i)) errcount = errcount + 1
         if (nlist3_f(i) .ne. nlist3_c(i)) errcount = errcount + 1
         if (nlist4_f(i) .ne. nlist4_c(i)) errcount = errcount + 1
      enddo
      do i = 1, mn1*nboxes
         if (list1_f(i) .ne. list1_c(i)) errcount = errcount + 1
      enddo
      do i = 1, mn2*nboxes
         if (list2_f(i) .ne. list2_c(i)) errcount = errcount + 1
      enddo
      do i = 1, mn3*nboxes
         if (list3_f(i) .ne. list3_c(i)) errcount = errcount + 1
      enddo
      do i = 1, mn4*nboxes
         if (list4_f(i) .ne. list4_c(i)) errcount = errcount + 1
      enddo
      call report('computelists       ', errcount, nfail)

c     =================================================================
c     updateflags - set a few boxes to iflag=3 and compare.
c     =================================================================
      allocate(iflag_f(nboxes), iflag_c(nboxes))
      do i = 1, nboxes
         iflag_f(i) = 0
         iflag_c(i) = 0
      enddo
c     Flag every 5th box at level 1 (if any) as 3, using the laddr
c     extracted above. updateflags uses boxsize(curlev) and
c     boxsize(curlev+1), so curlev must be strictly less than nlevels.
      curlev = 1
      if (nlevels .lt. 2) curlev = 0
      if (nlevels .lt. 1) then
         write(*,*) 'Tree too shallow to test updateflags'
         stop 2
      endif
      do i = laddr(2*curlev+1), laddr(2*curlev+2), 5
         iflag_f(i) = 3
         iflag_c(i) = 3
      enddo
      call updateflags(curlev, nboxes, nlevels, laddr, nchild, ichild,
     1     nnbors_f, nbors_f, tcenters, boxsize, iflag_f)
      call updateflags_c(curlev, nboxes, nlevels, laddr, nchild, ichild,
     1     nnbors_f, nbors_f, tcenters, boxsize, iflag_c)
      errcount = 0
      do i = 1, nboxes
         if (iflag_f(i) .ne. iflag_c(i)) errcount = errcount + 1
      enddo
      call report('updateflags        ', errcount, nfail)

c     =================================================================
c     tree_refine_boxes - synthetic small input.
c
c     We start with 4 existing boxes at level 0 (treated as the
c     "current" boxes). irefinebox = [1, 0, 1, 0] so boxes 1 and 3
c     are refined, each producing 4 new children that go at slots
c     nbctr+1..nbctr+8. nbctr is advanced by 2*4 = 8.
c     =================================================================
      do i = 1, nb_s
         centers_f(1, i) = 0.0d0
         centers_f(2, i) = 0.0d0
         centers_c(1, i) = 0.0d0
         centers_c(2, i) = 0.0d0
         ilevel_f(i) = 0
         ilevel_c(i) = 0
         iparent_f(i) = -1
         iparent_c(i) = -1
         nchild_f(i) = 0
         nchild_c(i) = 0
         do j = 1, 4
            ichild_f(j, i) = -1
            ichild_c(j, i) = -1
         enddo
      enddo

c     Set up 4 pre-existing boxes with distinct centers
      centers_f(1, 1) = 0.0d0
      centers_f(2, 1) = 0.0d0
      centers_f(1, 2) = 1.0d0
      centers_f(2, 2) = 0.0d0
      centers_f(1, 3) = 0.0d0
      centers_f(2, 3) = 1.0d0
      centers_f(1, 4) = 1.0d0
      centers_f(2, 4) = 1.0d0
      do i = 1, 4
         centers_c(1, i) = centers_f(1, i)
         centers_c(2, i) = centers_f(2, i)
         ilevel_f(i) = 2
         ilevel_c(i) = 2
         iparent_f(i) = 0
         iparent_c(i) = 0
         nchild_f(i) = 0
         nchild_c(i) = 0
      enddo

      irefinebox(1) = 1
      irefinebox(2) = 0
      irefinebox(3) = 1
      irefinebox(4) = 0

      ifirstbox = 1
      nbloc = 4
      nbctr_f = 4
      nbctr_c = 4
      nlctr = 3
      bs = 0.5d0

      call tree_refine_boxes(irefinebox, nb_s, ifirstbox, nbloc,
     1     centers_f, bs, nbctr_f, nlctr,
     2     ilevel_f, iparent_f, nchild_f, ichild_f)
      call tree_refine_boxes_c(irefinebox, nb_s, ifirstbox, nbloc,
     1     centers_c, bs, nbctr_c, nlctr,
     2     ilevel_c, iparent_c, nchild_c, ichild_c)

      errcount = 0
      if (nbctr_f .ne. nbctr_c) errcount = errcount + 1
      do i = 1, nb_s
         if (centers_f(1, i) .ne. centers_c(1, i))
     1        errcount = errcount + 1
         if (centers_f(2, i) .ne. centers_c(2, i))
     1        errcount = errcount + 1
         if (ilevel_f(i) .ne. ilevel_c(i)) errcount = errcount + 1
         if (iparent_f(i) .ne. iparent_c(i)) errcount = errcount + 1
         if (nchild_f(i) .ne. nchild_c(i)) errcount = errcount + 1
         do j = 1, 4
            if (ichild_f(j, i) .ne. ichild_c(j, i))
     1           errcount = errcount + 1
         enddo
      enddo
      call report('tree_refine_boxes  ', errcount, nfail)

c     =================================================================
c     tree_refine_boxes_flag - synthetic small input. Boxes 1 and 3
c     get iflag=1 so they refine and children inherit iflag=3.
c     Box 4 gets iflag=2 so it refines and children inherit iflag=0.
c     =================================================================
      do i = 1, nb_s
         centers_f(1, i) = 0.0d0
         centers_f(2, i) = 0.0d0
         centers_c(1, i) = 0.0d0
         centers_c(2, i) = 0.0d0
         ilevel_f(i) = 0
         ilevel_c(i) = 0
         iparent_f(i) = -1
         iparent_c(i) = -1
         nchild_f(i) = 0
         nchild_c(i) = 0
         iflag_s_f(i) = 0
         iflag_s_c(i) = 0
         do j = 1, 4
            ichild_f(j, i) = -1
            ichild_c(j, i) = -1
         enddo
      enddo

      centers_f(1, 1) = 0.0d0
      centers_f(2, 1) = 0.0d0
      centers_f(1, 2) = 1.0d0
      centers_f(2, 2) = 0.0d0
      centers_f(1, 3) = 0.0d0
      centers_f(2, 3) = 1.0d0
      centers_f(1, 4) = 1.0d0
      centers_f(2, 4) = 1.0d0
      do i = 1, 4
         centers_c(1, i) = centers_f(1, i)
         centers_c(2, i) = centers_f(2, i)
         ilevel_f(i) = 2
         ilevel_c(i) = 2
         iparent_f(i) = 0
         iparent_c(i) = 0
         nchild_f(i) = 0
         nchild_c(i) = 0
      enddo
      iflag_s_f(1) = 1
      iflag_s_f(2) = 0
      iflag_s_f(3) = 1
      iflag_s_f(4) = 2
      do i = 1, nb_s
         iflag_s_c(i) = iflag_s_f(i)
      enddo

      ifirstbox = 1
      nbloc = 4
      nbctr_f = 4
      nbctr_c = 4
      nlctr = 2
      bs = 0.5d0

      call tree_refine_boxes_flag(iflag_s_f, nb_s, ifirstbox, nbloc,
     1     centers_f, bs, nbctr_f, nlctr,
     2     ilevel_f, iparent_f, nchild_f, ichild_f)
      call tree_refine_boxes_flag_c(iflag_s_c, nb_s, ifirstbox, nbloc,
     1     centers_c, bs, nbctr_c, nlctr,
     2     ilevel_c, iparent_c, nchild_c, ichild_c)

      errcount = 0
      if (nbctr_f .ne. nbctr_c) errcount = errcount + 1
      do i = 1, nb_s
         if (centers_f(1, i) .ne. centers_c(1, i))
     1        errcount = errcount + 1
         if (centers_f(2, i) .ne. centers_c(2, i))
     1        errcount = errcount + 1
         if (ilevel_f(i) .ne. ilevel_c(i)) errcount = errcount + 1
         if (iparent_f(i) .ne. iparent_c(i)) errcount = errcount + 1
         if (nchild_f(i) .ne. nchild_c(i)) errcount = errcount + 1
         if (iflag_s_f(i) .ne. iflag_s_c(i)) errcount = errcount + 1
         do j = 1, 4
            if (ichild_f(j, i) .ne. ichild_c(j, i))
     1           errcount = errcount + 1
         enddo
      enddo
      call report('tree_refine_flag   ', errcount, nfail)

c     =================================================================
c     tree_copy - straight copy of all 5 fields.
c     =================================================================
      allocate(centers2_f(2, nboxes), centers2_c(2, nboxes))
      allocate(ilevel2_f(nboxes), ilevel2_c(nboxes))
      allocate(iparent2_f(nboxes), iparent2_c(nboxes))
      allocate(nchild2_f(nboxes), nchild2_c(nboxes))
      allocate(ichild2_f(4*nboxes), ichild2_c(4*nboxes))
      do i = 1, nboxes
         centers2_f(1, i) = -99.0d0
         centers2_f(2, i) = -99.0d0
         centers2_c(1, i) = -99.0d0
         centers2_c(2, i) = -99.0d0
         ilevel2_f(i) = -1
         ilevel2_c(i) = -1
         iparent2_f(i) = -1
         iparent2_c(i) = -1
         nchild2_f(i) = -1
         nchild2_c(i) = -1
      enddo
      do i = 1, 4*nboxes
         ichild2_f(i) = -1
         ichild2_c(i) = -1
      enddo

      call tree_copy(nboxes, tcenters, ilevel_tree, iparent, nchild,
     1     ichild, centers2_f, ilevel2_f, iparent2_f, nchild2_f,
     2     ichild2_f)
      call tree_copy_c(nboxes, tcenters, ilevel_tree, iparent, nchild,
     1     ichild, centers2_c, ilevel2_c, iparent2_c, nchild2_c,
     2     ichild2_c)
      errcount = 0
      do i = 1, nboxes
         if (centers2_f(1, i) .ne. centers2_c(1, i))
     1        errcount = errcount + 1
         if (centers2_f(2, i) .ne. centers2_c(2, i))
     1        errcount = errcount + 1
         if (ilevel2_f(i) .ne. ilevel2_c(i)) errcount = errcount + 1
         if (iparent2_f(i) .ne. iparent2_c(i))
     1        errcount = errcount + 1
         if (nchild2_f(i) .ne. nchild2_c(i)) errcount = errcount + 1
      enddo
      do i = 1, 4*nboxes
         if (ichild2_f(i) .ne. ichild2_c(i)) errcount = errcount + 1
      enddo
      call report('tree_copy          ', errcount, nfail)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all tree_routs2d cases match'

      end


c     ----------------------------------------------------------------
      subroutine report(name, errcount, nfail)
      implicit none
      character*(*) name
      integer errcount, nfail
      if (errcount .eq. 0) then
         write(*,1000) name, errcount
      else
         write(*,1001) name, errcount
         nfail = nfail + 1
      endif
 1000 format(' [ ok ] ',a19,'  errs=',i8)
 1001 format(' [FAIL] ',a19,'  errs=',i8)
      return
      end
