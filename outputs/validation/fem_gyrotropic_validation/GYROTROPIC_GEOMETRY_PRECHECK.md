# Gyrotropic geometry precheck

No production code was changed.

## Answer

The production circulator material region uses plasma radius **0.230 a = 4.60 mm**.

`sixport_common.r_plasma = 0.250 a` is not the radius of that material cylinder. It is stored as trainable-rod metadata and is not read back when the bulb is built.

## Where 0.250 a comes from

`scripts/validation/sixport_common.py` sets

```
r_plasma = 0.005 / a
```

with `a = 0.020 m`. That is 5.00 mm, and the comment says "~5 mm plasma radius".

`build` of the six-port device passes this value to `PMMI.Rod_Array_Hexagon_train(..., r=r_plasma, bulbs=False)`. `Add_Rod_train` only appends `{"r": r, "center": center}` to `train_elems`. The only later use of `train_elems` is `assert len(train_elems) == 91`. Nothing reads `train_elems[i]["r"]` into a Meep cylinder.

## Where 0.230 a comes from

The mounted bulb calls `PlasMEEP.Add_Bulb(..., profile=0)` in `plasmeep/lib.py`. Profile 0 appends three cylinders:

- quartz at `r_bulb_outer`
- vacuum at `r_bulb_inner`
- plasma at `4.6 * r_bulb_inner / 6.5`

`r_bulb_inner = 0.0065 / a = 0.325` (13 mm inner diameter). Then

```
4.6 * 0.325 / 6.5 = 0.230
```

which is 4.60 mm. The same factor is used by the quartz-only dielectric fill in `sixport_common._add_bulb_dielectric_fill`. The B=0 campaign's `R_CORE = 0.230` is this production core, not a separate invention.

## What the production geometry uses

Configuration `r_plasma` → array placement metadata only.

Configuration `(r_bulb_inner, r_bulb_outer)` → `Add_Bulb` → the plasma `mp.Cylinder` of radius `4.6 * r_bulb_inner / 6.5 = 0.230 a`.

The intended physical core for the circulator bulb is **0.230 a**. The gyrotropic ladder keeps that radius when it reaches a bulb geometry. `0.250 a` is not substituted for it.
