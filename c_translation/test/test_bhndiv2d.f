c     test_bhndiv2d.f - differential test for the C translation of
c     bhndiv2d. Mirrors test_hndiv2d.f.

      program test_bhndiv2d
      implicit none

      external bhndiv2d
      external bhndiv2d_c

      integer ntests
      parameter (ntests = 12)

      real *8 eps_list(ntests)
      integer i, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg
      integer ndiv_f, idivflag_f, ndiv_c, idivflag_c
      integer nfail

      data eps_list /
     1   1.0d0,    0.5d0,    0.1d0,    0.05d0,
     2   0.01d0,   0.001d0,  1.0d-6,   1.0d-9,
     3   1.0d-12,  1.0d-15,  1.0d-16,  1.0d-20 /

      nfail = 0

      do i = 1, ntests
         ns        = 100 + i
         nt        = 200 + 7*i
         ifcharge  = mod(i, 2)
         ifdipole  = mod(i+1, 2)
         ifpgh     = 1 + mod(i, 3)
         ifpghtarg = 1 + mod(i+1, 3)

         ndiv_f    = -999
         idivflag_f = -999
         ndiv_c    = -999
         idivflag_c = -999

         call bhndiv2d(eps_list(i), ns, nt, ifcharge, ifdipole,
     1        ifpgh, ifpghtarg, ndiv_f, idivflag_f)
         call bhndiv2d_c(eps_list(i), ns, nt, ifcharge, ifdipole,
     1        ifpgh, ifpghtarg, ndiv_c, idivflag_c)

         if (ndiv_f .ne. ndiv_c .or. idivflag_f .ne. idivflag_c) then
            write(*,1000) i, eps_list(i),
     1           ndiv_f, ndiv_c, idivflag_f, idivflag_c
            nfail = nfail + 1
         else
            write(*,1001) i, eps_list(i), ndiv_f, idivflag_f
         endif
      enddo

 1000 format(' [FAIL] case ',i3,' eps=',1pe10.3,
     1   '  ndiv f/c=',i6,'/',i6,'  idivflag f/c=',i3,'/',i3)
 1001 format(' [ ok ] case ',i3,' eps=',1pe10.3,
     1   '  ndiv=',i6,'  idivflag=',i3)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' of ', ntests, ' cases differ'
         stop 1
      endif
      write(*,*) 'PASS: all ', ntests, ' cases match'

      end
