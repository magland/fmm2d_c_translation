c     test_fmmcommon2d.f - differential test for the C translation of
c     selected routines from src/common/fmmcommon2d.f.
c
c     Calls both the Fortran reference (from libfmm2d.a) and the C
c     translation (from build/fmmcommon2d.o, exported with a _c_
c     suffix) on a battery of inputs and verifies bit-for-bit
c     equality of the outputs. The reordering routines do pure
c     copies (no arithmetic) so exact equality is appropriate;
c     init_carray uses only integer-valued additions so the result
c     is exactly representable in double precision and exact
c     equality is again appropriate.

      program test_fmmcommon2d
      implicit none

      external dreorderf, dreorderf_c
      external dreorderi, dreorderi_c
      external init_carray, init_carray_c
      real *8 hkrand
      external hkrand

      integer ndimmax, nmax
      parameter (ndimmax = 4, nmax = 1000)

      real *8 arr(ndimmax, nmax)
      real *8 sort_f(ndimmax, nmax), sort_c(ndimmax, nmax)
      integer iarr(nmax)

      integer ldcmax, csize
      parameter (ldcmax = 60)
      parameter (csize = (ldcmax+1)*(ldcmax+1))
      real *8 carray_f(csize), carray_c(csize)

      integer ndims(2), ns(2)
      data ndims / 2, 4 /
      data ns    / 100, 1000 /

      integer ldcs(3)
      data ldcs / 10, 30, 60 /

      integer is, idim_idx, ndim, n, ldc, csz
      integer i, j, itmp, itest
      integer nbad, nfail
      real *8 dummy

c     seed hkrand once
      dummy = hkrand(1234)

      nfail = 0

c     ---------------- dreorderf / dreorderi ----------------
      do idim_idx = 1, 2
         ndim = ndims(idim_idx)
         do is = 1, 2
            n = ns(is)

c           fill arr with random doubles
            do i = 1, n
               do j = 1, ndim
                  arr(j, i) = hkrand(0)
               enddo
            enddo

c           build a random permutation of 1..n by Fisher-Yates
            do i = 1, n
               iarr(i) = i
            enddo
            do i = n, 2, -1
               j = int(hkrand(0)*dble(i)) + 1
               if (j .gt. i) j = i
               itmp = iarr(i)
               iarr(i) = iarr(j)
               iarr(j) = itmp
            enddo

c           ---- dreorderf ----
            do i = 1, n
               do j = 1, ndim
                  sort_f(j, i) = -1.0d300
                  sort_c(j, i) = -1.0d300
               enddo
            enddo

            call dreorderf(ndim, n, arr, sort_f, iarr)
            call dreorderf_c(ndim, n, arr, sort_c, iarr)

            nbad = 0
            do i = 1, n
               do j = 1, ndim
                  if (sort_f(j, i) .ne. sort_c(j, i)) then
                     nbad = nbad + 1
                  endif
               enddo
            enddo

            if (nbad .ne. 0) then
               write(*,1000) 'dreorderf  ', ndim, n, nbad
               nfail = nfail + 1
            else
               write(*,1001) 'dreorderf  ', ndim, n
            endif

c           ---- dreorderi ----
            do i = 1, n
               do j = 1, ndim
                  sort_f(j, i) = -1.0d300
                  sort_c(j, i) = -1.0d300
               enddo
            enddo

            call dreorderi(ndim, n, arr, sort_f, iarr)
            call dreorderi_c(ndim, n, arr, sort_c, iarr)

            nbad = 0
            do i = 1, n
               do j = 1, ndim
                  if (sort_f(j, i) .ne. sort_c(j, i)) then
                     nbad = nbad + 1
                  endif
               enddo
            enddo

            if (nbad .ne. 0) then
               write(*,1000) 'dreorderi  ', ndim, n, nbad
               nfail = nfail + 1
            else
               write(*,1001) 'dreorderi  ', ndim, n
            endif

         enddo
      enddo

c     ---------------- init_carray ----------------
      do itest = 1, 3
         ldc = ldcs(itest)
         csz = (ldc+1)*(ldc+1)

         do i = 1, csz
            carray_f(i) = -1.0d300
            carray_c(i) = -1.0d300
         enddo

         call init_carray(carray_f, ldc)
         call init_carray_c(carray_c, ldc)

         nbad = 0
         do i = 1, csz
            if (carray_f(i) .ne. carray_c(i)) nbad = nbad + 1
         enddo

         if (nbad .ne. 0) then
            write(*,1002) 'init_carray', ldc, nbad
            nfail = nfail + 1
         else
            write(*,1003) 'init_carray', ldc
         endif
      enddo

 1000 format(' [FAIL] ',a11,' ndim=',i3,' n=',i6,
     1   '  ',i7,' element(s) differ')
 1001 format(' [ ok ] ',a11,' ndim=',i3,' n=',i6,
     1   '  all elements match')
 1002 format(' [FAIL] ',a11,' ldc=',i4,
     1   '  ',i7,' element(s) differ')
 1003 format(' [ ok ] ',a11,' ldc=',i4,
     1   '  all elements match')

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' test case(s) failed'
         stop 1
      endif
      write(*,*) 'PASS: all fmmcommon2d cases match'

      end
