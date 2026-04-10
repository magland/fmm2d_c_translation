c     test_cdjseval2d.f - differential test for the C translation of
c     src/common/cdjseval2d.f (jbessel2d).

      program test_cdjseval2d
      implicit none

      external jbessel2d, jbessel2d_c
      real *8 hkrand
      external hkrand

      integer ntmax
      parameter (ntmax = 100)

      complex *16 z
      real *8 rscale
      complex *16, allocatable :: fjs_f(:), fjs_c(:)
      complex *16, allocatable :: fjder_f(:), fjder_c(:)
      integer nterms, ifder, i, icase
      integer nfail
      real *8 dummy, errmax, e

      nfail = 0
      dummy = hkrand(1234)

c     Test various z values, rscales, and nterms
      do icase = 1, 20
         z = dcmplx(10.0d0*hkrand(0)-5.0d0,
     1              10.0d0*hkrand(0)-5.0d0)
         if (cdabs(z) .lt. 0.01d0)
     1      z = dcmplx(1.0d0, 1.0d0)
         nterms = 10 + int(80*hkrand(0))
         rscale = 1.0d0
         if (cdabs(z) .lt. 1.0d0) rscale = cdabs(z)

         do ifder = 0, 1
            allocate(fjs_f(0:nterms+1))
            allocate(fjs_c(0:nterms+1))
            allocate(fjder_f(0:nterms+1))
            allocate(fjder_c(0:nterms+1))

            do i = 0, nterms+1
               fjs_f(i) = (0.0d0, 0.0d0)
               fjs_c(i) = (0.0d0, 0.0d0)
               fjder_f(i) = (0.0d0, 0.0d0)
               fjder_c(i) = (0.0d0, 0.0d0)
            enddo

            call jbessel2d(nterms, z, rscale,
     1           fjs_f, ifder, fjder_f)
            call jbessel2d_c(nterms, z, rscale,
     1           fjs_c, ifder, fjder_c)

            errmax = 0.0d0
            do i = 0, nterms
               e = cdabs(fjs_f(i) - fjs_c(i))
               if (e .gt. errmax) errmax = e
            enddo

            if (ifder .eq. 1) then
               do i = 0, nterms
                  e = cdabs(fjder_f(i) - fjder_c(i))
                  if (e .gt. errmax) errmax = e
               enddo
            endif

            if (errmax .gt. 0.0d0) then
               write(*,1000) icase, ifder, errmax
               nfail = nfail + 1
            endif

            deallocate(fjs_f, fjs_c, fjder_f, fjder_c)
         enddo
      enddo

 1000 format(' [FAIL] jbessel2d case ',i3,
     1   ' ifder=',i2,' errmax=',1pe12.5)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: jbessel2d all 40 cases match'

      end
