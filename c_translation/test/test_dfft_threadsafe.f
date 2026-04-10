c     test_dfft_threadsafe.f - differential test for dfft_threadsafe.c
c
c     Tests the complex FFT routines (zffti, zfftf, zfftb) by calling
c     both the Fortran original and the C translation on identical
c     inputs and comparing bit-for-bit.

      program test_dfft_threadsafe
      implicit none

      external hkrand
      real *8 hkrand

      external zffti, zffti_c
      external zfftf, zfftf_c
      external zfftb, zfftb_c

      integer nfail, ntot
      real *8 dummy

c     workspace and data arrays
      integer nmax
      parameter (nmax = 256)
      real *8 wsave_f(4*nmax+100), wsave_c(4*nmax+100)
      real *8 c_f(2*nmax), c_c(2*nmax)
      real *8 c_orig(2*nmax)

      integer n, i, ns
      integer test_sizes(12)
      real *8 errmax, e

      data test_sizes /2,3,4,5,6,8,12,15,16,30,64,256/

      nfail = 0
      ntot = 0
      dummy = hkrand(42)

c     ============================================================
c     Test zffti + zfftf: forward transform
c     ============================================================
      do ns = 1, 12
         n = test_sizes(ns)
         ntot = ntot + 1

c        Initialize both versions
         call zffti(n, wsave_f)
         call zffti_c(n, wsave_c)

c        Check twiddle factors (wa region) match bit-for-bit.
c        Skip ifac region — it stores integers via sequence association
c        (4-byte int in 8-byte double slots) so raw comparison is not
c        meaningful; the FFT output comparison tests it end-to-end.
         errmax = 0.0d0
         do i = 2*n+1, 4*n
            e = dabs(wsave_f(i) - wsave_c(i))
            if (e .gt. errmax) errmax = e
         enddo
         if (errmax .ne. 0.0d0) then
            write(*,*) '[FAIL] zffti twiddles n=', n,
     1           ' max diff=', errmax
            nfail = nfail + 1
         endif

c        Generate random input data
         do i = 1, 2*n
            c_orig(i) = hkrand(0)
         enddo

c        Copy to both work arrays
         do i = 1, 2*n
            c_f(i) = c_orig(i)
            c_c(i) = c_orig(i)
         enddo

c        Forward FFT with both versions
         call zfftf(n, c_f, wsave_f)
         call zfftf_c(n, c_c, wsave_c)

c        Check bit-for-bit equality
         errmax = 0.0d0
         do i = 1, 2*n
            e = dabs(c_f(i) - c_c(i))
            if (e .gt. errmax) errmax = e
         enddo
         if (errmax .ne. 0.0d0) then
            write(*,*) '[FAIL] zfftf n=', n,
     1           ' max diff=', errmax
            nfail = nfail + 1
         endif

c        ============================================================
c        Test zfftb: backward transform on the forward-transformed data
c        ============================================================
c        Copy forward-transformed result to both arrays
         do i = 1, 2*n
            c_f(i) = c_c(i)
         enddo

c        Backward FFT with both versions
         call zfftb(n, c_f, wsave_f)
         call zfftb_c(n, c_c, wsave_c)

c        Check bit-for-bit equality
         errmax = 0.0d0
         do i = 1, 2*n
            e = dabs(c_f(i) - c_c(i))
            if (e .gt. errmax) errmax = e
         enddo
         if (errmax .ne. 0.0d0) then
            write(*,*) '[FAIL] zfftb n=', n,
     1           ' max diff=', errmax
            nfail = nfail + 1
         endif

c        Verify round-trip: backward(forward(x)) = n*x
         errmax = 0.0d0
         do i = 1, 2*n
            e = dabs(c_c(i)/dble(n) - c_orig(i))
            if (e .gt. errmax) errmax = e
         enddo
         if (errmax .gt. 1.0d-12) then
            write(*,*) '[FAIL] round-trip n=', n,
     1           ' max err=', errmax
            nfail = nfail + 1
         endif
      enddo

c     ============================================================
c     Test n=1 edge case (should be no-op)
c     ============================================================
      ntot = ntot + 1
      c_f(1) = 3.14d0
      c_f(2) = 2.72d0
      c_c(1) = 3.14d0
      c_c(2) = 2.72d0
      call zffti(1, wsave_f)
      call zffti_c(1, wsave_c)
      call zfftf(1, c_f, wsave_f)
      call zfftf_c(1, c_c, wsave_c)
      if (c_f(1) .ne. c_c(1) .or. c_f(2) .ne. c_c(2)) then
         write(*,*) '[FAIL] zfftf n=1'
         nfail = nfail + 1
      endif

c     ============================================================
c     Summary
c     ============================================================
      if (nfail .ne. 0) then
         write(*,*) 'FAIL:', nfail, ' test(s) failed out of',ntot
         stop 1
      endif
      write(*,*) 'PASS: dfft_threadsafe all', ntot,
     1     ' test sizes match bit-for-bit'
      end
