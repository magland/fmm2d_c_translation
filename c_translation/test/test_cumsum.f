c     test_cumsum.f - differential test for the C translation of cumsum.
c
c     Calls both cumsum (Fortran reference, from libfmm2d.a) and
c     cumsum_c_ (C translation, from build/cumsum.o) on a battery of
c     integer arrays and verifies that every element of the resulting
c     prefix-sum vectors matches exactly. Also exercises cumsum1 and
c     cumsum_para directly.

      program test_cumsum
      implicit none

      external cumsum, cumsum_c
      external cumsum1, cumsum1_c
      external cumsum_para, cumsum_para_c
      real *8 hkrand
      external hkrand

      integer nsizes
      parameter (nsizes = 4)
      integer sizes(nsizes)
      data sizes / 10, 100, 1000, 15000 /

      integer nmax
      parameter (nmax = 15000)
      integer a(nmax)
      integer b_f(nmax), b_c(nmax)
      integer b1_f(nmax), b1_c(nmax)
      integer bp_f(nmax), bp_c(nmax)

      integer ndwork
      parameter (ndwork = 200)
      integer d_f(ndwork), d_c(ndwork)
      integer nd

      integer is, n, i, nfail, nbad
      real *8 dummy

c     seed hkrand once
      dummy = hkrand(1234)

      nfail = 0
      nd = ndwork

      do is = 1, nsizes
         n = sizes(is)

c        fill input array with a mix of positive and negative integers
         do i = 1, n
            a(i) = int(hkrand(0)*10000.0d0) - 5000
         enddo

c        ---- test cumsum (wrapper) ----
         do i = 1, n
            b_f(i) = -999999
            b_c(i) = -999999
         enddo

         call cumsum(n, a, b_f)
         call cumsum_c(n, a, b_c)

         nbad = 0
         do i = 1, n
            if (b_f(i) .ne. b_c(i)) nbad = nbad + 1
         enddo

         if (nbad .ne. 0) then
            write(*,1000) 'cumsum     ', n, nbad
            nfail = nfail + 1
         else
            write(*,1001) 'cumsum     ', n
         endif

c        ---- test cumsum1 (serial) ----
         do i = 1, n
            b1_f(i) = -999999
            b1_c(i) = -999999
         enddo

         call cumsum1(n, a, b1_f)
         call cumsum1_c(n, a, b1_c)

         nbad = 0
         do i = 1, n
            if (b1_f(i) .ne. b1_c(i)) nbad = nbad + 1
         enddo

         if (nbad .ne. 0) then
            write(*,1000) 'cumsum1    ', n, nbad
            nfail = nfail + 1
         else
            write(*,1001) 'cumsum1    ', n
         endif

c        ---- test cumsum_para (parallel-but-serialized) ----
         do i = 1, n
            bp_f(i) = -999999
            bp_c(i) = -999999
         enddo
         do i = 1, ndwork
            d_f(i) = 0
            d_c(i) = 0
         enddo

         call cumsum_para(n, a, bp_f, nd, d_f)
         call cumsum_para_c(n, a, bp_c, nd, d_c)

         nbad = 0
         do i = 1, n
            if (bp_f(i) .ne. bp_c(i)) nbad = nbad + 1
         enddo

         if (nbad .ne. 0) then
            write(*,1000) 'cumsum_para', n, nbad
            nfail = nfail + 1
         else
            write(*,1001) 'cumsum_para', n
         endif

      enddo

 1000 format(' [FAIL] ',a11,' n=',i7,'  ',i7,' element(s) differ')
 1001 format(' [ ok ] ',a11,' n=',i7,'  all elements match')

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all cumsum cases match'

      end
