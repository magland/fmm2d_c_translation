c     test_next235.f - differential test for next235.c

      program test_next235
      implicit none

      integer next235, next235_c
      external next235, next235_c

      integer nfail, i
      real *8 base
      integer res_f, res_c

      nfail = 0

      do i = 1, 50
         base = dble(i) * 2.0d0
         res_f = next235(base)
         res_c = next235_c(base)
         if (res_f .ne. res_c) then
            write(*,*) '[FAIL] base=', base,
     1           ' f=', res_f, ' c=', res_c
            nfail = nfail + 1
         endif
      enddo

      if (nfail .ne. 0) then
         write(*,*) 'FAIL:', nfail, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: next235 all 50 cases match'
      end
