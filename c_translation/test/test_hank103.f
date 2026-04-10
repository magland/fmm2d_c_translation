c     test_hank103.f - differential test for the C translation of
c     src/common/hank103.f.
c
c     Tests hank103, hank103u, hank103r, hank103l, hank103a, hank103p,
c     hanks103, hanks104 by calling both Fortran and C versions on
c     identical inputs and comparing bit-for-bit.

      program test_hank103
      implicit none

      external hkrand
      real *8 hkrand

      external hank103,    hank103_c
      external hank103u,   hank103u_c
      external hank103r,   hank103r_c
      external hank103l,   hank103l_c
      external hank103a,   hank103a_c
      external hank103p,   hank103p_c
      external hanks103,   hanks103_c
      external hanks104,   hanks104_c

      integer nfail, i
      real *8 dummy, errmax, e
      complex *16 z, h0_f, h1_f, h0_c, h1_c
      integer ifexpon, ier_f, ier_c

c     For hanks103/104 tests
      integer nmax
      parameter (nmax = 30)
      complex *16 hanks_f(0:nmax), hanks_c(0:nmax)
      integer n, j

c     For hank103p test
      complex *16 poly(35), zarg, f_f, f_c
      integer m

      nfail = 0
      dummy = hkrand(1234)

c     ============================================================
c     Test hank103: multiple z values, both ifexpon=0 and 1
c     ============================================================
      do i = 1, 20
         z = dcmplx(10.0d0*hkrand(0)-5.0d0,
     1              10.0d0*hkrand(0)-5.0d0)
         do ifexpon = 0, 1
            h0_f = (0.0d0, 0.0d0)
            h1_f = (0.0d0, 0.0d0)
            h0_c = (0.0d0, 0.0d0)
            h1_c = (0.0d0, 0.0d0)

            call hank103(z, h0_f, h1_f, ifexpon)
            call hank103_c(z, h0_c, h1_c, ifexpon)

            errmax = 0.0d0
            e = cdabs(h0_f - h0_c)
            if (e .gt. errmax) errmax = e
            e = cdabs(h1_f - h1_c)
            if (e .gt. errmax) errmax = e

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hank103', i, ifexpon,
     1              errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hank103', 40

c     ============================================================
c     Test hank103u: z in upper half-plane
c     ============================================================
      do i = 1, 20
         z = dcmplx(20.0d0*hkrand(0)-10.0d0,
     1              10.0d0*hkrand(0))
         if (dimag(z) .lt. 0.01d0) z = z + (0.0d0, 0.1d0)
         do ifexpon = 0, 1
            h0_f = (0.0d0, 0.0d0)
            h1_f = (0.0d0, 0.0d0)
            h0_c = (0.0d0, 0.0d0)
            h1_c = (0.0d0, 0.0d0)

            call hank103u(z, ier_f, h0_f, h1_f,
     1           ifexpon)
            call hank103u_c(z, ier_c, h0_c, h1_c,
     1           ifexpon)

            errmax = 0.0d0
            e = cdabs(h0_f - h0_c)
            if (e .gt. errmax) errmax = e
            e = cdabs(h1_f - h1_c)
            if (e .gt. errmax) errmax = e
            if (ier_f .ne. ier_c) errmax = 999.0d0

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hank103u', i,
     1              ifexpon, errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hank103u', 40

c     ============================================================
c     Test hank103r: z in right lower quadrant
c     ============================================================
      do i = 1, 20
         z = dcmplx(10.0d0*hkrand(0),
     1              -10.0d0*hkrand(0))
         if (dble(z) .lt. 0.01d0) z = z + (0.1d0, 0.0d0)
         if (dimag(z) .gt. -0.01d0) z = z - (0.0d0, 0.1d0)
         do ifexpon = 0, 1
            h0_f = (0.0d0, 0.0d0)
            h1_f = (0.0d0, 0.0d0)
            h0_c = (0.0d0, 0.0d0)
            h1_c = (0.0d0, 0.0d0)

            call hank103r(z, ier_f, h0_f, h1_f,
     1           ifexpon)
            call hank103r_c(z, ier_c, h0_c, h1_c,
     1           ifexpon)

            errmax = 0.0d0
            e = cdabs(h0_f - h0_c)
            if (e .gt. errmax) errmax = e
            e = cdabs(h1_f - h1_c)
            if (e .gt. errmax) errmax = e
            if (ier_f .ne. ier_c) errmax = 999.0d0

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hank103r', i,
     1              ifexpon, errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hank103r', 40

c     ============================================================
c     Test hank103l: local regime (small |z|)
c     ============================================================
      do i = 1, 10
         z = dcmplx(0.5d0*hkrand(0), 0.5d0*hkrand(0))
         if (cdabs(z) .lt. 1.0d-10) z = (0.1d0, 0.1d0)
         do ifexpon = 0, 1
            h0_f = (0.0d0, 0.0d0)
            h1_f = (0.0d0, 0.0d0)
            h0_c = (0.0d0, 0.0d0)
            h1_c = (0.0d0, 0.0d0)

            call hank103l(z, h0_f, h1_f, ifexpon)
            call hank103l_c(z, h0_c, h1_c,
     1           ifexpon)

            errmax = 0.0d0
            e = cdabs(h0_f - h0_c)
            if (e .gt. errmax) errmax = e
            e = cdabs(h1_f - h1_c)
            if (e .gt. errmax) errmax = e

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hank103l', i,
     1              ifexpon, errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hank103l', 20

c     ============================================================
c     Test hank103a: asymptotic regime (large |z|, upper half)
c     ============================================================
      do i = 1, 10
         z = dcmplx(20.0d0+30.0d0*hkrand(0),
     1              5.0d0*hkrand(0))
         do ifexpon = 0, 1
            h0_f = (0.0d0, 0.0d0)
            h1_f = (0.0d0, 0.0d0)
            h0_c = (0.0d0, 0.0d0)
            h1_c = (0.0d0, 0.0d0)

            call hank103a(z, h0_f, h1_f, ifexpon)
            call hank103a_c(z, h0_c, h1_c,
     1           ifexpon)

            errmax = 0.0d0
            e = cdabs(h0_f - h0_c)
            if (e .gt. errmax) errmax = e
            e = cdabs(h1_f - h1_c)
            if (e .gt. errmax) errmax = e

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hank103a', i,
     1              ifexpon, errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hank103a', 20

c     ============================================================
c     Test hank103p: polynomial evaluation
c     ============================================================
      do i = 1, 35
         poly(i) = dcmplx(hkrand(0), hkrand(0))
      enddo
      m = 35
      zarg = dcmplx(0.1d0*hkrand(0), 0.1d0*hkrand(0))
      f_f = (0.0d0, 0.0d0)
      f_c = (0.0d0, 0.0d0)
      call hank103p(poly, m, zarg, f_f)
      call hank103p_c(poly, m, zarg, f_c)
      errmax = cdabs(f_f - f_c)
      if (errmax .gt. 0.0d0) then
         write(*,1002) 'hank103p', errmax
         nfail = nfail + 1
      endif
      write(*,1001) 'hank103p', 1

c     ============================================================
c     Test hanks103: H_0..H_n via optimized recursion
c     ============================================================
      do i = 1, 10
         z = dcmplx(5.0d0*hkrand(0)+0.5d0,
     1              5.0d0*hkrand(0)+0.5d0)
         n = 20
         do ifexpon = 0, 1
            do j = 0, n
               hanks_f(j) = (0.0d0, 0.0d0)
               hanks_c(j) = (0.0d0, 0.0d0)
            enddo

            call hanks103(z, hanks_f, n, ifexpon)
            call hanks103_c(z, hanks_c, n,
     1           ifexpon)

            errmax = 0.0d0
            do j = 0, n
               e = cdabs(hanks_f(j) - hanks_c(j))
               if (e .gt. errmax) errmax = e
            enddo

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hanks103', i,
     1              ifexpon, errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hanks103', 20

c     ============================================================
c     Test hanks104: H_0..H_n via simple recursion
c     ============================================================
      do i = 1, 10
         z = dcmplx(5.0d0*hkrand(0)+0.5d0,
     1              5.0d0*hkrand(0)+0.5d0)
         n = 20
         do ifexpon = 0, 1
            do j = 0, n
               hanks_f(j) = (0.0d0, 0.0d0)
               hanks_c(j) = (0.0d0, 0.0d0)
            enddo

            call hanks104(z, hanks_f, n, ifexpon)
            call hanks104_c(z, hanks_c, n,
     1           ifexpon)

            errmax = 0.0d0
            do j = 0, n
               e = cdabs(hanks_f(j) - hanks_c(j))
               if (e .gt. errmax) errmax = e
            enddo

            if (errmax .gt. 0.0d0) then
               write(*,1000) 'hanks104', i,
     1              ifexpon, errmax
               nfail = nfail + 1
            endif
         enddo
      enddo
      write(*,1001) 'hanks104', 20

c     ============================================================
c     Summary
c     ============================================================
 1000 format(' [FAIL] ',a,' case ',i3,' ifexpon=',i2,
     1   ' errmax=',1pe12.5)
 1001 format(' [ ok ] ',a,': ',i4,' cases passed')
 1002 format(' [FAIL] ',a,' errmax=',1pe12.5)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: hank103 all tests passed'

      end
