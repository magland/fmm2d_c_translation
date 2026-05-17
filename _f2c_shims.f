c     _f2c_shims.f -- helper routines for f2c-translated fmm2d build.
c
c     Provides Fortran-77 implementations of F90 intrinsics that f2c
c     does not (yet) understand when applied to particular argument
c     forms, plus a `second()` timer that returns CPU seconds.
c
c     Sources are pre-processed by makefile.f2c (sed) to rewrite the
c     offending call sites into calls to these shims.
c
c----------------------------------------------------------------------
c
c     imaxvalh(a, n) -- equivalent of maxval(a(1:n)) for integer arrays.
c
        integer function imaxvalh(a, n)
        integer a(*), n, i, m
        if (n .le. 0) then
          imaxvalh = 0
          return
        endif
        m = a(1)
        do i = 2, n
          if (a(i) .gt. m) m = a(i)
        enddo
        imaxvalh = m
        return
        end
c
c----------------------------------------------------------------------
c
c     second() -- LAPACK-style CPU timer returning REAL*8 seconds.
c     Implemented in terms of the Fortran 90 cpu_time intrinsic, which
c     f2c does recognise (it expands to a libf2c clock_() call).
c
        real *8 function second()
        real *8 t
        call cpu_time(t)
        second = t
        return
        end
