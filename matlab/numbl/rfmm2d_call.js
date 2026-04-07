// native: fmm2d
// wasm: fmm2d
//
// rfmm2d_call(nd, eps, ns, sources, ifcharge, charges, ifdipole,
//             dipstr, dipvec, ifpgh, nt, targ, ifpghtarg) -> tensor
//
// Low-level shim that calls rfmm2d_w in fmm2d.{wasm,so,dylib}. Currently
// the only fmm2d entry point exposed via numbl. The matlab-facing
// rfmm2d.m overrides matlab/rfmm2d.m and calls this with flat,
// column-major buffers.
//
// All output buffers are always allocated full size (matching the
// existing matlab/rfmm2d.m pattern); the underlying Fortran routine
// only writes to the ones the requested ifpgh / ifpghtarg level needs.
//
// The shim returns a single packed Float64 tensor with this layout:
//   [ pot(nd*ns) |
//     grad(2*nd*ns) |
//     hess(3*nd*ns) |
//     pottarg(nd*ntuse) |
//     gradtarg(2*nd*ntuse) |
//     hesstarg(3*nd*ntuse) ]
// where ntuse = max(nt, 1). The matlab wrapper slices and reshapes.
register({
  resolve: function (argTypes, nargout) {
    if (argTypes.length !== 13) {
      return null;
    }
    return {
      outputTypes: [{ kind: "tensor", isComplex: false }],
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

        // Output sizes (col-major flat lengths). Always full-size,
        // matching the Fortran callee's expected dims.
        var n_pot = nd * ns;
        var n_grad = 2 * nd * ns;
        var n_hess = 3 * nd * ns;
        var n_pottarg = nd * ntuse;
        var n_gradtarg = 2 * nd * ntuse;
        var n_hesstarg = 3 * nd * ntuse;
        var n_total = n_pot + n_grad + n_hess + n_pottarg + n_gradtarg + n_hesstarg;

        // Input sizes.
        var n_sources = 2 * ns;
        var n_charges = nd * ns;
        var n_dipstr = nd * ns;
        var n_dipvec = 2 * nd * ns;
        var n_targ = 2 * ntuse;

        // Helpers to pull a flat Float64Array out of a tensor argument.
        // Inputs that are "not used" will still be passed as full-size
        // zero arrays by the .m wrapper, so we don't need to handle null
        // here, but we do tolerate length mismatches by zero-padding.
        function flat(tensor, n) {
          var out = new Float64Array(n);
          if (tensor && tensor.data) {
            var src = tensor.data;
            var k = Math.min(src.length, n);
            for (var i = 0; i < k; i++) out[i] = src[i];
          }
          return out;
        }

        var sourcesArr = flat(sources, n_sources);
        var chargesArr = flat(charges, n_charges);
        var dipstrArr = flat(dipstr, n_dipstr);
        var dipvecArr = flat(dipvec, n_dipvec);
        var targArr = flat(targ, n_targ);

        if (native) {
          var fn = native.func(
            "int rfmm2d_w(int nd, double eps, int ns, double *sources," +
              " int ifcharge, double *charge, int ifdipole," +
              " double *dipstr, double *dipvec," +
              " int ifpgh, double *pot, double *grad, double *hess," +
              " int nt, double *targ, int ifpghtarg," +
              " double *pottarg, double *gradtarg, double *hesstarg)"
          );

          var pot = new Float64Array(n_pot);
          var grad = new Float64Array(n_grad);
          var hess = new Float64Array(n_hess);
          var pottarg = new Float64Array(n_pottarg);
          var gradtarg = new Float64Array(n_gradtarg);
          var hesstarg = new Float64Array(n_hesstarg);

          var ier = fn(
            nd, eps, ns, sourcesArr,
            ifcharge, chargesArr, ifdipole, dipstrArr, dipvecArr,
            ifpgh, pot, grad, hess,
            nt, targArr, ifpghtarg, pottarg, gradtarg, hesstarg
          );
          if (ier !== 0) {
            throw new RuntimeError("rfmm2d_w failed with error code " + ier);
          }

          var packed = new FloatXArray(n_total);
          var off = 0;
          packed.set(pot, off); off += n_pot;
          packed.set(grad, off); off += n_grad;
          packed.set(hess, off); off += n_hess;
          packed.set(pottarg, off); off += n_pottarg;
          packed.set(gradtarg, off); off += n_gradtarg;
          packed.set(hesstarg, off);
          return RTV.tensor(packed, [n_total, 1]);
        }

        // ── WASM path ─────────────────────────────────────────────────
        var BYTES = 8;
        var exports = wasm.exports;
        var mem = exports.memory;

        // Allocate every input and output buffer in linear memory.
        var sources_ptr = exports.my_malloc(n_sources * BYTES);
        var charges_ptr = exports.my_malloc(n_charges * BYTES);
        var dipstr_ptr = exports.my_malloc(n_dipstr * BYTES);
        var dipvec_ptr = exports.my_malloc(n_dipvec * BYTES);
        var targ_ptr = exports.my_malloc(n_targ * BYTES);

        var pot_ptr = exports.my_malloc(n_pot * BYTES);
        var grad_ptr = exports.my_malloc(n_grad * BYTES);
        var hess_ptr = exports.my_malloc(n_hess * BYTES);
        var pottarg_ptr = exports.my_malloc(n_pottarg * BYTES);
        var gradtarg_ptr = exports.my_malloc(n_gradtarg * BYTES);
        var hesstarg_ptr = exports.my_malloc(n_hesstarg * BYTES);

        // Copy inputs into WASM memory. Re-acquire the view after every
        // my_malloc call in case ALLOW_MEMORY_GROWTH grew the buffer.
        var view = new Float64Array(mem.buffer);
        view.set(sourcesArr, sources_ptr / BYTES);
        view.set(chargesArr, charges_ptr / BYTES);
        view.set(dipstrArr, dipstr_ptr / BYTES);
        view.set(dipvecArr, dipvec_ptr / BYTES);
        view.set(targArr, targ_ptr / BYTES);
        // Output buffers don't need to be zeroed; the FMM writes them
        // (or leaves them alone if the corresponding flag is off and
        // the matlab wrapper ignores them).

        var ier = exports.rfmm2d_w(
          nd, eps, ns, sources_ptr,
          ifcharge, charges_ptr, ifdipole, dipstr_ptr, dipvec_ptr,
          ifpgh, pot_ptr, grad_ptr, hess_ptr,
          nt, targ_ptr, ifpghtarg, pottarg_ptr, gradtarg_ptr, hesstarg_ptr
        );

        // Re-acquire the view after the call too.
        view = new Float64Array(mem.buffer);

        var packed = new FloatXArray(n_total);
        var off = 0;
        packed.set(view.subarray(pot_ptr / BYTES, pot_ptr / BYTES + n_pot), off);
        off += n_pot;
        packed.set(view.subarray(grad_ptr / BYTES, grad_ptr / BYTES + n_grad), off);
        off += n_grad;
        packed.set(view.subarray(hess_ptr / BYTES, hess_ptr / BYTES + n_hess), off);
        off += n_hess;
        packed.set(view.subarray(pottarg_ptr / BYTES, pottarg_ptr / BYTES + n_pottarg), off);
        off += n_pottarg;
        packed.set(view.subarray(gradtarg_ptr / BYTES, gradtarg_ptr / BYTES + n_gradtarg), off);
        off += n_gradtarg;
        packed.set(view.subarray(hesstarg_ptr / BYTES, hesstarg_ptr / BYTES + n_hesstarg), off);

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
          throw new RuntimeError("rfmm2d_w failed with error code " + ier);
        }

        return RTV.tensor(packed, [n_total, 1]);
      },
    };
  },
});
