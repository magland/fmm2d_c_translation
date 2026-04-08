// native: fmm2d
// wasm: fmm2d
//
// lfmm2d_call(nd, eps, ns, sources, ifcharge, charges, ifdipole,
//             dipstr, dipvec, ifpgh, nt, targ, ifpghtarg) -> complex tensor
//
// Low-level shim that calls lfmm2d_w in fmm2d.{wasm,so,dylib}.
// The matlab-facing lfmm2d.m overrides matlab/lfmm2d.m and calls this
// with column-major buffers. charges and dipstr are complex tensors;
// sources, dipvec, and targ are real tensors.
//
// All output buffers are always allocated full size; the underlying
// Fortran routine only writes to the ones the requested ifpgh /
// ifpghtarg level needs.
//
// Layout differences vs rfmm2d:
//   - charge / dipstr are complex (each element is a (re, im) pair)
//   - dipvec is REAL with shape (2*nd, ns) — same as rfmm2d
//   - pot has shape (nd, ns) complex; grad (nd, 2, ns) complex;
//     hess (nd, 3, ns) complex.
//
// The shim returns a single packed complex tensor (real and imag parts
// each length n_total) with this layout:
//   [ pot(nd*ns) | grad(2*nd*ns) | hess(3*nd*ns) |
//     pottarg(nd*ntuse) | gradtarg(2*nd*ntuse) | hesstarg(3*nd*ntuse) ]
// where ntuse = max(nt, 1). The matlab wrapper slices and reshapes.
//
// gfortran complex*16 ABI: complex arrays are stored as interleaved
// (re, im, re, im, ...) doubles in linear memory. The JS shim
// interleaves on input and deinterleaves on output.
register({
  resolve: function (argTypes, nargout) {
    if (argTypes.length !== 13) {
      return null;
    }
    return {
      outputTypes: [{ kind: "tensor", isComplex: true }],
      apply: function (args, nargout) {
        var nd = args[0];
        var eps = args[1];
        var ns = args[2];
        var sources = args[3];
        var ifcharge = args[4];
        var charges = args[5];
        var ifdipole = args[6];
        var dipstr = args[7];
        var dipvec = args[8];
        var ifpgh = args[9];
        var nt = args[10];
        var targ = args[11];
        var ifpghtarg = args[12];

        var ntuse = nt > 0 ? nt : 1;

        // Element counts (1 element = 1 complex number).
        var n_pot = nd * ns;
        var n_grad = 2 * nd * ns;
        var n_hess = 3 * nd * ns;
        var n_pottarg = nd * ntuse;
        var n_gradtarg = 2 * nd * ntuse;
        var n_hesstarg = 3 * nd * ntuse;
        var n_total = n_pot + n_grad + n_hess + n_pottarg + n_gradtarg + n_hesstarg;

        // Input element counts.
        var n_sources = 2 * ns;          // real
        var n_charges = nd * ns;         // complex
        var n_dipstr = nd * ns;          // complex
        var n_dipvec = 2 * nd * ns;      // real (same as rfmm2d)
        var n_targ = 2 * ntuse;          // real

        // Helper: pull a real flat array from a real tensor.
        function flatReal(tensor, n) {
          var out = new Float64Array(n);
          if (tensor && tensor.data) {
            var src = tensor.data;
            var k = Math.min(src.length, n);
            for (var i = 0; i < k; i++) out[i] = src[i];
          }
          return out;
        }

        // Helper: interleave a complex tensor into (re, im, re, im, ...)
        // doubles. n is the number of complex elements; output length is 2*n.
        function flatComplex(tensor, n) {
          var out = new Float64Array(2 * n);
          if (tensor && tensor.data) {
            var re = tensor.data;
            var im = tensor.imag;
            var k = Math.min(re.length, n);
            for (var i = 0; i < k; i++) {
              out[2 * i] = re[i];
              out[2 * i + 1] = im ? im[i] : 0;
            }
          }
          return out;
        }

        var sourcesArr = flatReal(sources, n_sources);
        var chargesArr = flatComplex(charges, n_charges);
        var dipstrArr = flatComplex(dipstr, n_dipstr);
        var dipvecArr = flatReal(dipvec, n_dipvec);
        var targArr = flatReal(targ, n_targ);

        // Helper: deinterleave a (re, im, re, im, ...) Float64Array
        // segment of length 2*n into separate re/im arrays at offset
        // off in the output Float* arrays.
        function deinterleave(view, base, n, packed_re, packed_im, off) {
          for (var i = 0; i < n; i++) {
            packed_re[off + i] = view[base + 2 * i];
            packed_im[off + i] = view[base + 2 * i + 1];
          }
        }

        if (native) {
          var fn = native.func(
            "int lfmm2d_w(int nd, double eps, int ns, double *sources," +
              " int ifcharge, double *charge, int ifdipole," +
              " double *dipstr, double *dipvec," +
              " int ifpgh, double *pot, double *grad, double *hess," +
              " int nt, double *targ, int ifpghtarg," +
              " double *pottarg, double *gradtarg, double *hesstarg)"
          );

          var pot = new Float64Array(2 * n_pot);
          var grad = new Float64Array(2 * n_grad);
          var hess = new Float64Array(2 * n_hess);
          var pottarg = new Float64Array(2 * n_pottarg);
          var gradtarg = new Float64Array(2 * n_gradtarg);
          var hesstarg = new Float64Array(2 * n_hesstarg);

          var ier = fn(
            nd, eps, ns, sourcesArr,
            ifcharge, chargesArr, ifdipole, dipstrArr, dipvecArr,
            ifpgh, pot, grad, hess,
            nt, targArr, ifpghtarg, pottarg, gradtarg, hesstarg
          );
          if (ier !== 0) {
            throw new RuntimeError("lfmm2d_w failed with error code " + ier);
          }

          var packed_re = new FloatXArray(n_total);
          var packed_im = new FloatXArray(n_total);
          var off = 0;
          deinterleave(pot, 0, n_pot, packed_re, packed_im, off); off += n_pot;
          deinterleave(grad, 0, n_grad, packed_re, packed_im, off); off += n_grad;
          deinterleave(hess, 0, n_hess, packed_re, packed_im, off); off += n_hess;
          deinterleave(pottarg, 0, n_pottarg, packed_re, packed_im, off); off += n_pottarg;
          deinterleave(gradtarg, 0, n_gradtarg, packed_re, packed_im, off); off += n_gradtarg;
          deinterleave(hesstarg, 0, n_hesstarg, packed_re, packed_im, off);
          return RTV.tensor(packed_re, [n_total, 1], packed_im);
        }

        // ── WASM path ─────────────────────────────────────────────────
        var BYTES = 8;
        var exports = wasm.exports;
        var mem = exports.memory;

        var sources_ptr = exports.my_malloc(n_sources * BYTES);
        var charges_ptr = exports.my_malloc(2 * n_charges * BYTES);
        var dipstr_ptr = exports.my_malloc(2 * n_dipstr * BYTES);
        var dipvec_ptr = exports.my_malloc(n_dipvec * BYTES);
        var targ_ptr = exports.my_malloc(n_targ * BYTES);

        var pot_ptr = exports.my_malloc(2 * n_pot * BYTES);
        var grad_ptr = exports.my_malloc(2 * n_grad * BYTES);
        var hess_ptr = exports.my_malloc(2 * n_hess * BYTES);
        var pottarg_ptr = exports.my_malloc(2 * n_pottarg * BYTES);
        var gradtarg_ptr = exports.my_malloc(2 * n_gradtarg * BYTES);
        var hesstarg_ptr = exports.my_malloc(2 * n_hesstarg * BYTES);

        var view = new Float64Array(mem.buffer);
        view.set(sourcesArr, sources_ptr / BYTES);
        view.set(chargesArr, charges_ptr / BYTES);
        view.set(dipstrArr, dipstr_ptr / BYTES);
        view.set(dipvecArr, dipvec_ptr / BYTES);
        view.set(targArr, targ_ptr / BYTES);

        var ier = exports.lfmm2d_w(
          nd, eps, ns, sources_ptr,
          ifcharge, charges_ptr, ifdipole, dipstr_ptr, dipvec_ptr,
          ifpgh, pot_ptr, grad_ptr, hess_ptr,
          nt, targ_ptr, ifpghtarg, pottarg_ptr, gradtarg_ptr, hesstarg_ptr
        );

        view = new Float64Array(mem.buffer);

        var packed_re = new FloatXArray(n_total);
        var packed_im = new FloatXArray(n_total);
        var off = 0;
        deinterleave(view, pot_ptr / BYTES, n_pot, packed_re, packed_im, off); off += n_pot;
        deinterleave(view, grad_ptr / BYTES, n_grad, packed_re, packed_im, off); off += n_grad;
        deinterleave(view, hess_ptr / BYTES, n_hess, packed_re, packed_im, off); off += n_hess;
        deinterleave(view, pottarg_ptr / BYTES, n_pottarg, packed_re, packed_im, off); off += n_pottarg;
        deinterleave(view, gradtarg_ptr / BYTES, n_gradtarg, packed_re, packed_im, off); off += n_gradtarg;
        deinterleave(view, hesstarg_ptr / BYTES, n_hesstarg, packed_re, packed_im, off);

        exports.my_free(sources_ptr);
        exports.my_free(charges_ptr);
        exports.my_free(dipstr_ptr);
        exports.my_free(dipvec_ptr);
        exports.my_free(targ_ptr);
        exports.my_free(pot_ptr);
        exports.my_free(grad_ptr);
        exports.my_free(hess_ptr);
        exports.my_free(pottarg_ptr);
        exports.my_free(gradtarg_ptr);
        exports.my_free(hesstarg_ptr);

        if (ier !== 0) {
          throw new RuntimeError("lfmm2d_w failed with error code " + ier);
        }

        return RTV.tensor(packed_re, [n_total, 1], packed_im);
      },
    };
  },
});
