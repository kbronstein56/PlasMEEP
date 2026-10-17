# Free-91 IndexError (recorded, not prioritized)

## Symptom

```
File multimode_opt_core.py, gradient_q_adjoint
    o_masks[m],
IndexError: list index out of range
```

during `refine_free91_multimode.py` first `eval_and_grad`.

## Cause (trivial)

`refine_free91_multimode.py` builds **91 singleton orbits** and passes them as `q`/`orbits`
to `eval_and_grad`.

Inside `eval_and_grad`, orbit masks are rebuilt via:

```python
_, _, o_masks = orbits_and_masks(rod_masks=st0["rod_masks"])
```

`orbits_and_masks()` **always** constructs the C6 tied orbits (~16), ignoring the
caller’s singleton list. So `len(o_masks)≈16` while `gradient_q_adjoint` loops
`m = 0 .. 90` → IndexError on `o_masks[m]`.

## Fix (deferred)

Build `o_masks` from the **provided** `orbits` argument (or add
`orbits_and_masks(..., orbits=...)`), then re-run free-91 only after power
balance + M-opt + M→F are done.

## Status

**Deprioritized** per campaign instructions. No compute spent beyond this note.
