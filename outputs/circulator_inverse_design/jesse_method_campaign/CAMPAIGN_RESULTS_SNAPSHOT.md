# Campaign results snapshot

Generated from artifacts in `jesse_method_campaign`.

## Forward-mode vs adjoint

- Adjoint wall: **126.1 s**, LU=1, RHS=6+6
- Forward full est: **149.7 s**, sens RHS=96, rel err vs adj=1.999e-14, pass=True
- Scaling C6 (n=16): adj~48s fwd~70s → prefer **adjoint**
- Scaling C3 (n=31): adj~83s fwd~126s → prefer **adjoint**
- Scaling free91 (n=91): adj~221s fwd~352s → prefer **adjoint**

## Mode A vs B baselines

| design | J_A | J_B | des0 | des1 | iso TE0 | iso tot |
|--------|----:|----:|-----:|-----:|--------:|--------:|
| uniform | -0.4876 | -0.4610 | 0.0108 | 0.0193 | 7.84 | 1.82 |
| opt_uniform_C_warm_M | 0.1936 | 0.2415 | 0.5810 | 0.0270 | 9.10 | 8.47 |
| opt_uniform_C | 0.7082 | 0.7221 | 0.8041 | 0.0041 | 19.00 | 18.67 |

## Spectral 200MHz
- n=11 worst iso=-24.1755926092997 dB @ 5.8 GHz

## Spectral 100MHz
- n=21 worst iso=-33.071581038049814 dB @ 5.9 GHz

## Spectral 50MHz
- n=41 worst iso=-52.82874138800315 dB @ 5.75 GHz

## Grid comparison
```json
{
  "200_vs_100": {
    "d_worst_iso_dB": 8.895988428750115,
    "d_worst_IL_dB": 2.687343480710446,
    "d_worst_RL_dB": 0.05007565507218237,
    "d_worst_iso_location_GHz": 0.10000000000000053,
    "pass_iso": false,
    "pass_IL": false,
    "pass_RL": true,
    "pass_location": false,
    "pass_all": false
  },
  "100_vs_50": {
    "d_worst_iso_dB": 19.757160349953338,
    "d_worst_IL_dB": 14.584359529389637,
    "d_worst_RL_dB": 0.30470895384399516,
    "d_worst_iso_location_GHz": 0.15000000000000036,
    "pass_iso": false,
    "pass_IL": false,
    "pass_RL": false,
    "pass_location": false,
    "pass_all": false
  }
}
```

## Rotation source check
- C6 max rel err=0.050163875859541414, pass=False

## Soft-min E2E FD
- J_BB=-0.643303294635704, weight FD pass=True, overall pass=True

## Mode A/B short opt

- A: J_best=0.7274568295397404, nit=5
- B: J_best=0.7548485585136101, nit=5

## C3 vs C6 short opt
- map |Δs|=1.1102230246251565e-16
- C6: n_var=16 J_best=0.7274568295397404
- C3: n_var=31 J_best=0.4441399147431724

