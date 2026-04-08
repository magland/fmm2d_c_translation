// native: fmm2d
// wasm: fmm2d
//
// stfmm2d_call(nd, eps, ns, sources, ifstoklet, stoklet, ifstrslet,
//              strslet, strsvec, ifppreg, nt, targ, ifppregtarg) -> tensor
//
// Low-level shim that calls stfmm2d_w in fmm2d.{wasm,so,dylib}.
// Stokes FMM. All inputs/outputs are REAL doubles (no complex).
//
// Output buffers are always allocated full size; the underlying
// Fortran routine only writes to the ones the requested ifppreg /
// ifppregtarg level needs.
//
// The shim returns a single packed Float64 tensor with this layout:
//   [ pot(2*nd*ns) | pre(nd*ns) | grad(4*nd*ns) |
//     pottarg(2*nd*ntuse) | pretarg(nd*ntuse) | gradtarg(4*nd*ntuse) ]
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
        var ifstoklet = args[4];
        var stoklet = args[5];
        var ifstrslet = args[6];
        var strslet = args[7];
        var strsvec = args[8];
        var ifppreg = args[9];
        var nt = args[10];
        var targ = args[11];
        var ifppregtarg = args[12];

        var ntuse = nt > 0 ? nt : 1;

        // Output sizes (col-major flat lengths). Always full-size.
        // pot/pottarg are (nd, 2, *) -> 2*nd per element
        // pre/pretarg are (nd, *)    -> nd per element
        // grad/gradtarg are (nd, 2, 2, *) -> 4*nd per element
        var n_pot = 2 * nd * ns;
        var n_pre = nd * ns;
        var n_grad = 4 * nd * ns;
        var n_pottarg = 2 * nd * ntuse;
        var n_pretarg = nd * ntuse;
        var n_gradtarg = 4 * nd * ntuse;
        var n_total = n_pot + n_pre + n_grad
                    + n_pottarg + n_pretarg + n_gradtarg;

        // Input sizes (real flat lengths).
        var n_sources = 2 * ns;
        var n_stoklet = 2 * nd * ns;
        var n_strslet = 2 * nd * ns;
        var n_strsvec = 2 * nd * ns;
        var n_targ = 2 * ntuse;

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
        var stokletArr = flat(stoklet, n_stoklet);
        var strsletArr = flat(strslet, n_strslet);
        var strsvecArr = flat(strsvec, n_strsvec);
        var targArr = flat(targ, n_targ);

        if (native) {
          var fn = native.func(
            "int stfmm2d_w(int nd, double eps, int ns, double *source," +
              " int ifstoklet, double *stoklet," +
              " int ifstrslet, double *strslet, double *strsvec," +
              " int ifppreg, double *pot, double *pre, double *grad," +
              " int nt, double *targ, int ifppregtarg," +
              " double *pottarg, double *pretarg, double *gradtarg)"
          );

          var pot = new Float64Array(n_pot);
          var pre = new Float64Array(n_pre);
          var grad = new Float64Array(n_grad);
          var pottarg = new Float64Array(n_pottarg);
          var pretarg = new Float64Array(n_pretarg);
          var gradtarg = new Float64Array(n_gradtarg);

          var ier = fn(
            nd, eps, ns, sourcesArr,
            ifstoklet, stokletArr,
            ifstrslet, strsletArr, strsvecArr,
            ifppreg, pot, pre, grad,
            nt, targArr, ifppregtarg, pottarg, pretarg, gradtarg
          );
          if (ier !== 0) {
            throw new RuntimeError("stfmm2d_w failed with error code " + ier);
          }

          var packed = new FloatXArray(n_total);
          var off = 0;
          packed.set(pot, off); off += n_pot;
          packed.set(pre, off); off += n_pre;
          packed.set(grad, off); off += n_grad;
          packed.set(pottarg, off); off += n_pottarg;
          packed.set(pretarg, off); off += n_pretarg;
          packed.set(gradtarg, off);
          return RTV.tensor(packed, [n_total, 1]);
        }

        // ── WASM path ─────────────────────────────────────────────────
        var BYTES = 8;
        var exports = wasm.exports;
        var mem = exports.memory;

        var sources_ptr = exports.my_malloc(n_sources * BYTES);
        var stoklet_ptr = exports.my_malloc(n_stoklet * BYTES);
        var strslet_ptr = exports.my_malloc(n_strslet * BYTES);
        var strsvec_ptr = exports.my_malloc(n_strsvec * BYTES);
        var targ_ptr = exports.my_malloc(n_targ * BYTES);

        var pot_ptr = exports.my_malloc(n_pot * BYTES);
        var pre_ptr = exports.my_malloc(n_pre * BYTES);
        var grad_ptr = exports.my_malloc(n_grad * BYTES);
        var pottarg_ptr = exports.my_malloc(n_pottarg * BYTES);
        var pretarg_ptr = exports.my_malloc(n_pretarg * BYTES);
        var gradtarg_ptr = exports.my_malloc(n_gradtarg * BYTES);

        var view = new Float64Array(mem.buffer);
        view.set(sourcesArr, sources_ptr / BYTES);
        view.set(stokletArr, stoklet_ptr / BYTES);
        view.set(strsletArr, strslet_ptr / BYTES);
        view.set(strsvecArr, strsvec_ptr / BYTES);
        view.set(targArr, targ_ptr / BYTES);

        var ier = exports.stfmm2d_w(
          nd, eps, ns, sources_ptr,
          ifstoklet, stoklet_ptr,
          ifstrslet, strslet_ptr, strsvec_ptr,
          ifppreg, pot_ptr, pre_ptr, grad_ptr,
          nt, targ_ptr, ifppregtarg, pottarg_ptr, pretarg_ptr, gradtarg_ptr
        );

        view = new Float64Array(mem.buffer);

        var packed = new FloatXArray(n_total);
        var off = 0;
        packed.set(view.subarray(pot_ptr / BYTES,
                                 pot_ptr / BYTES + n_pot), off);
        off += n_pot;
        packed.set(view.subarray(pre_ptr / BYTES,
                                 pre_ptr / BYTES + n_pre), off);
        off += n_pre;
        packed.set(view.subarray(grad_ptr / BYTES,
                                 grad_ptr / BYTES + n_grad), off);
        off += n_grad;
        packed.set(view.subarray(pottarg_ptr / BYTES,
                                 pottarg_ptr / BYTES + n_pottarg), off);
        off += n_pottarg;
        packed.set(view.subarray(pretarg_ptr / BYTES,
                                 pretarg_ptr / BYTES + n_pretarg), off);
        off += n_pretarg;
        packed.set(view.subarray(gradtarg_ptr / BYTES,
                                 gradtarg_ptr / BYTES + n_gradtarg), off);

        exports.my_free(sources_ptr);
        exports.my_free(stoklet_ptr);
        exports.my_free(strslet_ptr);
        exports.my_free(strsvec_ptr);
        exports.my_free(targ_ptr);
        exports.my_free(pot_ptr);
        exports.my_free(pre_ptr);
        exports.my_free(grad_ptr);
        exports.my_free(pottarg_ptr);
        exports.my_free(pretarg_ptr);
        exports.my_free(gradtarg_ptr);

        if (ier !== 0) {
          throw new RuntimeError("stfmm2d_w failed with error code " + ier);
        }

        return RTV.tensor(packed, [n_total, 1]);
      },
    };
  },
});
