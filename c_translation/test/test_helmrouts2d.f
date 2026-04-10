c     test_helmrouts2d.f - differential test for helmrouts2d.c

      program test_helmrouts2d
      implicit none

      external h2dformmpc,  h2dformmpc_c
      external h2dformmpd,  h2dformmpd_c
      external h2dformmpcd, h2dformmpcd_c
      external h2dmpevalp,  h2dmpevalp_c
      external h2dmpevalg,  h2dmpevalg_c
      external h2dmpevalh,  h2dmpevalh_c
      external h2dformtac,  h2dformtac_c
      external h2dformtad,  h2dformtad_c
      external h2dformtacd, h2dformtacd_c
      external h2dtaevalp,  h2dtaevalp_c
      external h2dtaevalg,  h2dtaevalg_c
      external h2dtaevalh,  h2dtaevalh_c
      external h2dmpmp,     h2dmpmp_c
      external h2dlocloc,   h2dlocloc_c
      external h2dmploc,    h2dmploc_c
      real *8 hkrand
      external hkrand

      integer nsmax, ntmax
      parameter (nsmax=20, ntmax=15)
      integer nd, nterms
      parameter (nd=2, nterms=12)

      real *8 sources(2,nsmax), targ(2,ntmax)
      real *8 center(2), center2(2)
      complex *16 zk
      complex *16 charge(nd,nsmax), dipstr(nd,nsmax)
      real *8 dipvec(nd,2,nsmax), rscale, rscale2
      integer ns, nt, i, j, ii
      integer nfail
      real *8 dummy, errmax, e

      integer msize
      parameter (msize = nd*(2*nterms+1))
      complex *16 mp_f(nd,-nterms:nterms)
      complex *16 mp_c(nd,-nterms:nterms)
      complex *16 loc_f(nd,-nterms:nterms)
      complex *16 loc_c(nd,-nterms:nterms)
      complex *16 pot_f(nd,ntmax), pot_c(nd,ntmax)
      complex *16 grad_f(nd,2,ntmax), grad_c(nd,2,ntmax)
      complex *16 hess_f(nd,3,ntmax), hess_c(nd,3,ntmax)

      ns = nsmax
      nt = ntmax
      nfail = 0

      zk = dcmplx(3.0d0, 0.05d0)
      rscale = 1.0d0
      rscale2 = 1.0d0
      center(1) = 0.0d0
      center(2) = 0.0d0
      center2(1) = 3.0d0
      center2(2) = 0.0d0

      dummy = hkrand(1234)
      do i = 1, ns
         sources(1,i) = 0.1d0*(hkrand(0)-0.5d0)
         sources(2,i) = 0.1d0*(hkrand(0)-0.5d0)
      enddo
      do j = 1, nt
         targ(1,j) = 3.0d0+0.1d0*(hkrand(0)-0.5d0)
         targ(2,j) = 0.1d0*(hkrand(0)-0.5d0)
      enddo
      do i = 1, ns
         do ii = 1, nd
            charge(ii,i) = dcmplx(hkrand(0),hkrand(0))
            dipstr(ii,i) = dcmplx(hkrand(0),hkrand(0))
            dipvec(ii,1,i) = hkrand(0)
            dipvec(ii,2,i) = hkrand(0)
         enddo
      enddo

c     ---- h2dformmpc ----
      call zeromp(mp_f, msize)
      call zeromp(mp_c, msize)
      call h2dformmpc(nd, zk, rscale, sources, ns,
     1     charge, center, nterms, mp_f)
      call h2dformmpc_c(nd, zk, rscale, sources,
     1     ns, charge, center, nterms, mp_c)
      call chkmp(mp_f, mp_c, msize, 'h2dformmpc',
     1     nfail)

c     ---- h2dformmpd ----
      call zeromp(mp_f, msize)
      call zeromp(mp_c, msize)
      call h2dformmpd(nd, zk, rscale, sources, ns,
     1     dipstr, dipvec, center, nterms, mp_f)
      call h2dformmpd_c(nd, zk, rscale, sources,
     1     ns, dipstr, dipvec, center, nterms,
     1     mp_c)
      call chkmp(mp_f, mp_c, msize, 'h2dformmpd',
     1     nfail)

c     ---- h2dformmpcd ----
      call zeromp(mp_f, msize)
      call zeromp(mp_c, msize)
      call h2dformmpcd(nd, zk, rscale, sources,
     1     ns, charge, dipstr, dipvec, center,
     1     nterms, mp_f)
      call h2dformmpcd_c(nd, zk, rscale, sources,
     1     ns, charge, dipstr, dipvec, center,
     1     nterms, mp_c)
      call chkmp(mp_f, mp_c, msize,'h2dformmpcd',
     1     nfail)

c     Use the Fortran multipole for eval tests
c     ---- h2dmpevalp ----
      call zeromp(pot_f, nd*nt)
      call zeromp(pot_c, nd*nt)
      call h2dmpevalp(nd, zk, rscale, center,
     1     mp_f, nterms, targ, nt, pot_f)
      call h2dmpevalp_c(nd, zk, rscale, center,
     1     mp_f, nterms, targ, nt, pot_c)
      call chkmp(pot_f, pot_c, nd*nt,
     1     'h2dmpevalp', nfail)

c     ---- h2dmpevalg ----
      call zeromp(pot_f, nd*nt)
      call zeromp(pot_c, nd*nt)
      call zeromp(grad_f, nd*2*nt)
      call zeromp(grad_c, nd*2*nt)
      call h2dmpevalg(nd, zk, rscale, center,
     1     mp_f, nterms, targ, nt, pot_f, grad_f)
      call h2dmpevalg_c(nd, zk, rscale, center,
     1     mp_f, nterms, targ, nt, pot_c, grad_c)
      call chkmp(pot_f, pot_c, nd*nt,
     1     'h2dmpevalg:p', nfail)
      call chkmp(grad_f, grad_c, nd*2*nt,
     1     'h2dmpevalg:g', nfail)

c     ---- h2dmpevalh ----
      call zeromp(pot_f, nd*nt)
      call zeromp(pot_c, nd*nt)
      call zeromp(grad_f, nd*2*nt)
      call zeromp(grad_c, nd*2*nt)
      call zeromp(hess_f, nd*3*nt)
      call zeromp(hess_c, nd*3*nt)
      call h2dmpevalh(nd, zk, rscale, center,
     1     mp_f, nterms, targ, nt, pot_f,
     1     grad_f, hess_f)
      call h2dmpevalh_c(nd, zk, rscale, center,
     1     mp_f, nterms, targ, nt, pot_c,
     1     grad_c, hess_c)
      call chkmp(pot_f, pot_c, nd*nt,
     1     'h2dmpevalh:p', nfail)
      call chkmp(grad_f, grad_c, nd*2*nt,
     1     'h2dmpevalh:g', nfail)
      call chkmp(hess_f, hess_c, nd*3*nt,
     1     'h2dmpevalh:h', nfail)

c     ---- h2dformtac ----
      call zeromp(loc_f, msize)
      call zeromp(loc_c, msize)
      call h2dformtac(nd, zk, rscale, sources, ns,
     1     charge, center2, nterms, loc_f)
      call h2dformtac_c(nd, zk, rscale, sources,
     1     ns, charge, center2, nterms, loc_c)
      call chkmp(loc_f, loc_c, msize,
     1     'h2dformtac', nfail)

c     ---- h2dformtad ----
      call zeromp(loc_f, msize)
      call zeromp(loc_c, msize)
      call h2dformtad(nd, zk, rscale, sources, ns,
     1     dipstr, dipvec, center2, nterms, loc_f)
      call h2dformtad_c(nd, zk, rscale, sources,
     1     ns, dipstr, dipvec, center2, nterms,
     1     loc_c)
      call chkmp(loc_f, loc_c, msize,
     1     'h2dformtad', nfail)

c     ---- h2dformtacd ----
      call zeromp(loc_f, msize)
      call zeromp(loc_c, msize)
      call h2dformtacd(nd, zk, rscale, sources,
     1     ns, charge, dipstr, dipvec, center2,
     1     nterms, loc_f)
      call h2dformtacd_c(nd, zk, rscale, sources,
     1     ns, charge, dipstr, dipvec, center2,
     1     nterms, loc_c)
      call chkmp(loc_f, loc_c, msize,
     1     'h2dformtacd', nfail)

c     ---- h2dtaevalp ----
      call zeromp(pot_f, nd*nt)
      call zeromp(pot_c, nd*nt)
      call h2dtaevalp(nd, zk, rscale, center2,
     1     loc_f, nterms, targ, nt, pot_f)
      call h2dtaevalp_c(nd, zk, rscale, center2,
     1     loc_f, nterms, targ, nt, pot_c)
      call chkmp(pot_f, pot_c, nd*nt,
     1     'h2dtaevalp', nfail)

c     ---- h2dtaevalg ----
      call zeromp(pot_f, nd*nt)
      call zeromp(pot_c, nd*nt)
      call zeromp(grad_f, nd*2*nt)
      call zeromp(grad_c, nd*2*nt)
      call h2dtaevalg(nd, zk, rscale, center2,
     1     loc_f, nterms, targ, nt, pot_f, grad_f)
      call h2dtaevalg_c(nd, zk, rscale, center2,
     1     loc_f, nterms, targ, nt, pot_c, grad_c)
      call chkmp(pot_f, pot_c, nd*nt,
     1     'h2dtaevalg:p', nfail)
      call chkmp(grad_f, grad_c, nd*2*nt,
     1     'h2dtaevalg:g', nfail)

c     ---- h2dtaevalh ----
      call zeromp(pot_f, nd*nt)
      call zeromp(pot_c, nd*nt)
      call zeromp(grad_f, nd*2*nt)
      call zeromp(grad_c, nd*2*nt)
      call zeromp(hess_f, nd*3*nt)
      call zeromp(hess_c, nd*3*nt)
      call h2dtaevalh(nd, zk, rscale, center2,
     1     loc_f, nterms, targ, nt, pot_f,
     1     grad_f, hess_f)
      call h2dtaevalh_c(nd, zk, rscale, center2,
     1     loc_f, nterms, targ, nt, pot_c,
     1     grad_c, hess_c)
      call chkmp(pot_f, pot_c, nd*nt,
     1     'h2dtaevalh:p', nfail)
      call chkmp(grad_f, grad_c, nd*2*nt,
     1     'h2dtaevalh:g', nfail)
      call chkmp(hess_f, hess_c, nd*3*nt,
     1     'h2dtaevalh:h', nfail)

c     ---- h2dmpmp ----
      call zeromp(mp_c, msize)
      call zeromp(loc_f, msize)
      call zeromp(loc_c, msize)
      call h2dmpmp(nd, zk, rscale, center, mp_f,
     1     nterms, rscale2, center2, loc_f,
     1     nterms)
      call h2dmpmp_c(nd, zk, rscale, center, mp_f,
     1     nterms, rscale2, center2, loc_c,
     1     nterms)
      call chkmp(loc_f, loc_c, msize, 'h2dmpmp',
     1     nfail)

c     ---- h2dmploc ----
      call zeromp(loc_f, msize)
      call zeromp(loc_c, msize)
      call h2dmploc(nd, zk, rscale, center, mp_f,
     1     nterms, rscale2, center2, loc_f,
     1     nterms)
      call h2dmploc_c(nd, zk, rscale, center,
     1     mp_f, nterms, rscale2, center2, loc_c,
     1     nterms)
      call chkmp(loc_f, loc_c, msize, 'h2dmploc',
     1     nfail)

c     ---- h2dlocloc ----
c     First form a proper local expansion
      call zeromp(loc_f, msize)
      call h2dformtac(nd, zk, rscale, sources, ns,
     1     charge, center2, nterms, loc_f)
      call zeromp(mp_f, msize)
      call zeromp(mp_c, msize)
      call h2dlocloc(nd, zk, rscale, center2,
     1     loc_f, nterms, rscale2, center, mp_f,
     1     nterms)
      call h2dlocloc_c(nd, zk, rscale, center2,
     1     loc_f, nterms, rscale2, center, mp_c,
     1     nterms)
      call chkmp(mp_f, mp_c, msize, 'h2dlocloc',
     1     nfail)

      if (nfail .ne. 0) then
         write(*,*) 'FAIL: ', nfail, ' checks differ'
         stop 1
      endif
      write(*,*) 'PASS: helmrouts2d all tests'

      end

      subroutine zeromp(buf, n)
      implicit none
      integer n, i
      complex *16 buf(n)
      do i = 1, n
         buf(i) = (0.0d0, 0.0d0)
      enddo
      return
      end

      subroutine chkmp(f, c, n, name, nfail)
      implicit none
      integer n, nfail, i
      complex *16 f(n), c(n)
      character*(*) name
      real *8 errmax, e
      errmax = 0.0d0
      do i = 1, n
         e = cdabs(f(i) - c(i))
         if (e .gt. errmax) errmax = e
      enddo
      if (errmax .gt. 0.0d0) then
         write(*,1000) name, errmax
         nfail = nfail + 1
      else
         write(*,1001) name
      endif
 1000 format(' [FAIL] ',a,' errmax=',1pe12.5)
 1001 format(' [ ok ] ',a)
      return
      end
