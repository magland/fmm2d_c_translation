c     test_mbhgreen2d.f - differential test for the C translation
c     of src/modified-biharmonic/mbhgreen2d.f

      program test_mbhgreen2d
      implicit none

      external modbhgreen_all
      external modbhgreen_all_c
      external modbhgreen
      external modbhgreen_c
      external modbhgreend1
      external modbhgreend1_c
      external modbhgreend2
      external modbhgreend2_c
      external difflogbk
      external difflogbk_c
      external diffslogbk
      external diffslogbk_c
      external diffslogbk_fast
      external diffslogbk_fast_c
      external diffszkik
      external diffszkik_c
      external diffszkik_fast
      external diffszkik_fast_c
      real *8 hkrand
      external hkrand

      integer nfail
      real *8 beta, r, dummy
      real *8 zx(2), zy(2), dir1(2), dir2(2)

c     outputs for modbhgreen_all
      real *8 pot_f, pot_c
      real *8 grad_f(2), grad_c(2)
      real *8 hess_f(3), hess_c(3)
      real *8 der3_f(4), der3_c(4)
      real *8 der4_f(5), der4_c(5)
      real *8 der5_f(6), der5_c(6)
      integer ifpot, ifgrad, ifhess
      integer ifder3, ifder4, ifder5

c     outputs for difflogbk
      real *8 g0_f, g0_c, g1_f, g1_c
      real *8 g2_f, g2_c, g3_f, g3_c
      integer if0, if1, if2, if3

c     outputs for diffslogbk / diffszkik
      integer nterms
      parameter (nterms = 5)
      real *8 diffs_f(0:nterms), diffs_c(0:nterms)
      real *8 ders_f(0:nterms), ders_c(0:nterms)
      real *8 kvec_f(0:nterms+4), kvec_c(0:nterms+4)
      real *8 ivec_f(0:nterms), ivec_c(0:nterms)
      real *8 rscale
      integer ifders

      real *8 errmax, e, tol
      integer i, j

      nfail = 0
      tol = 1.0d-13

c     seed
      dummy = hkrand(1234)

c     random source/target
      zx(1) = 1.0d0 + 0.5d0*hkrand(0)
      zx(2) = 2.0d0 + 0.5d0*hkrand(0)
      zy(1) = 0.0d0 + 0.5d0*hkrand(0)
      zy(2) = 0.0d0 + 0.5d0*hkrand(0)
      beta = 1.5d0 + hkrand(0)
      dir1(1) = hkrand(0)
      dir1(2) = hkrand(0)
      dir2(1) = hkrand(0)
      dir2(2) = hkrand(0)

c     ---- modbhgreen_all ----
      ifpot = 1
      ifgrad = 1
      ifhess = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 1
      call modbhgreen_all(beta,zx,zy,ifpot,
     1   pot_f,ifgrad,grad_f,ifhess,hess_f,
     2   ifder3,der3_f,ifder4,der4_f,
     3   ifder5,der5_f)
      call modbhgreen_all_c(beta,zx,zy,ifpot,
     1   pot_c,ifgrad,grad_c,ifhess,hess_c,
     2   ifder3,der3_c,ifder4,der4_c,
     3   ifder5,der5_c)
      errmax = 0.0d0
      e = dabs(pot_f - pot_c)
      if (e .gt. errmax) errmax = e
      do i = 1,2
         e = dabs(grad_f(i) - grad_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,3
         e = dabs(hess_f(i) - hess_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,4
         e = dabs(der3_f(i) - der3_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,5
         e = dabs(der4_f(i) - der4_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,6
         e = dabs(der5_f(i) - der5_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('green_all    ', errmax, nfail)

c     ---- modbhgreen ----
      call modbhgreen(beta,zx,zy,ifpot,
     1   pot_f,ifgrad,grad_f,ifhess,hess_f)
      call modbhgreen_c(beta,zx,zy,ifpot,
     1   pot_c,ifgrad,grad_c,ifhess,hess_c)
      errmax = 0.0d0
      e = dabs(pot_f - pot_c)
      if (e .gt. errmax) errmax = e
      do i = 1,2
         e = dabs(grad_f(i) - grad_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,3
         e = dabs(hess_f(i) - hess_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('green        ', errmax, nfail)

c     ---- modbhgreend1 ----
      call modbhgreend1(beta,zx,zy,ifpot,
     1   pot_f,ifgrad,grad_f,ifhess,hess_f,
     2   dir1)
      call modbhgreend1_c(beta,zx,zy,ifpot,
     1   pot_c,ifgrad,grad_c,ifhess,hess_c,
     2   dir1)
      errmax = 0.0d0
      e = dabs(pot_f - pot_c)
      if (e .gt. errmax) errmax = e
      do i = 1,2
         e = dabs(grad_f(i) - grad_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,3
         e = dabs(hess_f(i) - hess_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('greend1      ', errmax, nfail)

c     ---- modbhgreend2 ----
      call modbhgreend2(beta,zx,zy,ifpot,
     1   pot_f,ifgrad,grad_f,ifhess,hess_f,
     2   dir1,dir2)
      call modbhgreend2_c(beta,zx,zy,ifpot,
     1   pot_c,ifgrad,grad_c,ifhess,hess_c,
     2   dir1,dir2)
      errmax = 0.0d0
      e = dabs(pot_f - pot_c)
      if (e .gt. errmax) errmax = e
      do i = 1,2
         e = dabs(grad_f(i) - grad_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      do i = 1,3
         e = dabs(hess_f(i) - hess_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('greend2      ', errmax, nfail)

c     ---- difflogbk ----
      r = dsqrt((zx(1)-zy(1))**2+(zx(2)-zy(2))**2)
      if0 = 1
      if1 = 1
      if2 = 1
      if3 = 1
      call difflogbk(r,beta,if0,g0_f,if1,g1_f,
     1   if2,g2_f,if3,g3_f)
      call difflogbk_c(r,beta,if0,g0_c,if1,g1_c,
     1   if2,g2_c,if3,g3_c)
      errmax = 0.0d0
      e = dabs(g0_f-g0_c)
      if (e .gt. errmax) errmax = e
      e = dabs(g1_f-g1_c)
      if (e .gt. errmax) errmax = e
      e = dabs(g2_f-g2_c)
      if (e .gt. errmax) errmax = e
      e = dabs(g3_f-g3_c)
      if (e .gt. errmax) errmax = e
      call report('difflogbk    ', errmax, nfail)

c     ---- diffslogbk (small yh) ----
      rscale = 1.0d0
      r = 0.5d0
      call diffslogbk(r,beta,rscale,diffs_f,nterms)
      call diffslogbk_c(r,beta,rscale,diffs_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('diffslogbk_sm', errmax, nfail)

c     ---- diffslogbk (large yh) ----
      r = 5.0d0
      call diffslogbk(r,beta,rscale,diffs_f,nterms)
      call diffslogbk_c(r,beta,rscale,diffs_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('diffslogbk_lg', errmax, nfail)

c     ---- diffslogbk_fast (small yh, no ders) ----
      r = 0.5d0
      ifders = 0
      call diffslogbk_fast(r,beta,rscale,diffs_f,
     1   ifders,ders_f,kvec_f,nterms)
      call diffslogbk_fast_c(r,beta,rscale,diffs_c,
     1   ifders,ders_c,kvec_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('fastlogbk_sm ', errmax, nfail)

c     ---- diffslogbk_fast (large yh, with ders) ----
      r = 5.0d0
      ifders = 1
      call diffslogbk_fast(r,beta,rscale,diffs_f,
     1   ifders,ders_f,kvec_f,nterms)
      call diffslogbk_fast_c(r,beta,rscale,diffs_c,
     1   ifders,ders_c,kvec_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
         e = dabs(ders_f(i)-ders_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('fastlogbk_lg ', errmax, nfail)

c     ---- diffszkik (small yh) ----
      r = 0.5d0
      call diffszkik(r,beta,rscale,diffs_f,nterms)
      call diffszkik_c(r,beta,rscale,diffs_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('diffszkik_sm ', errmax, nfail)

c     ---- diffszkik (large yh) ----
      r = 5.0d0
      call diffszkik(r,beta,rscale,diffs_f,nterms)
      call diffszkik_c(r,beta,rscale,diffs_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('diffszkik_lg ', errmax, nfail)

c     ---- diffszkik_fast (small yh, with ders) ----
      r = 0.5d0
      ifders = 1
      call diffszkik_fast(r,beta,rscale,diffs_f,
     1   ifders,ders_f,ivec_f,nterms)
      call diffszkik_fast_c(r,beta,rscale,diffs_c,
     1   ifders,ders_c,ivec_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
         e = dabs(ders_f(i)-ders_c(i))
         if (e .gt. errmax) errmax = e
         e = dabs(ivec_f(i)-ivec_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('fastzkik_sm  ', errmax, nfail)

c     ---- diffszkik_fast (large yh, with ders) ----
      r = 5.0d0
      ifders = 1
      call diffszkik_fast(r,beta,rscale,diffs_f,
     1   ifders,ders_f,ivec_f,nterms)
      call diffszkik_fast_c(r,beta,rscale,diffs_c,
     1   ifders,ders_c,ivec_c,nterms)
      errmax = 0.0d0
      do i = 0,nterms
         e = dabs(diffs_f(i)-diffs_c(i))
         if (e .gt. errmax) errmax = e
         e = dabs(ders_f(i)-ders_c(i))
         if (e .gt. errmax) errmax = e
         e = dabs(ivec_f(i)-ivec_c(i))
         if (e .gt. errmax) errmax = e
      enddo
      call report('fastzkik_lg  ', errmax, nfail)

      if (nfail .gt. 0) then
         write(*,*) 'FAIL: ',nfail,
     1      ' mbhgreen2d tests failed'
         stop 1
      endif
      write(*,*) 'PASS: all mbhgreen2d cases match'

      end


c     -------------------------------------------------------
      subroutine report(name, errmax, nfail)
      implicit none
      character*(*) name
      integer nfail
      real *8 errmax, tol
      tol = 1.0d-13
      if (errmax .lt. tol) then
         write(*,1000) name, errmax
      else
         write(*,1001) name, errmax
         nfail = nfail + 1
      endif
 1000 format(' [ ok ] ',a13,
     1   '  errmax=',1pe11.3)
 1001 format(' [FAIL] ',a13,
     1   '  errmax=',1pe11.3)
      return
      end
