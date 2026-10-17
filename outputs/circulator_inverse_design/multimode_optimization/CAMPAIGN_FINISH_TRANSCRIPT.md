# Campaign finish transcript

```
===== FINISH SF CAMPAIGN 2026-10-15T18:48:56-07:00 =====
BEST_JSON=outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json
===== C6 opt report =====
wrote C6 report 3 runs
===== PHASE 6: full12 validation =====
Using MPI version 4.1, 1 processes
validating outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json grade=C
{
  "OPTIMIZED_FULL12_ONSAGER": "PASS",
  "r_F": 1.9832093813819298e-13,
  "OPTIMIZED_POWER_BALANCE": "FAIL",
  "max_residual": 0.1361038109251833,
  "elapsed_s": 4190.408073186874
}

Elapsed run time = 4191.2885 s
===== PHASE 7: fine-mesh C→M then M→F =====
Using MPI version 4.1, 1 processes
eval C ...
eval M ...
{
  "fine_ok": true,
  "elapsed_s": 8018.879281997681
}

Elapsed run time = 8019.5573 s
Using MPI version 4.1, 1 processes
eval M ...
eval F ...
{
  "fine_ok": true,
  "elapsed_s": 8514.3892288208
}

Elapsed run time = 8516.0894 s
===== PHASE 9: frequency characterization =====
Using MPI version 4.1, 1 processes
f=3.700 GHz n_prop/port≈1 near_cutoff=False ...
f=3.750 GHz n_prop/port≈2 near_cutoff=True ...
f=3.800 GHz n_prop/port≈2 near_cutoff=False ...
f=3.850 GHz n_prop/port≈2 near_cutoff=False ...
f=3.900 GHz n_prop/port≈2 near_cutoff=False ...
f=3.950 GHz n_prop/port≈2 near_cutoff=False ...
f=4.000 GHz n_prop/port≈2 near_cutoff=False ...
{
  "TE1_cutoff_GHz": 3.747405725,
  "recommendation": "NEEDS_CUTOFF_DERIVATION",
  "elapsed_s": 2538.82035946846
}

Elapsed run time = 2539.0424 s
===== PHASE 8: free-91 short =====
Using MPI version 4.1, 1 processes
Traceback (most recent call last):
  File "/home/daq_user/PlasMEEP/scripts/circulator_inverse_design/multimode_optimization/refine_free91_multimode.py", line 156, in <module>
    main()
    ~~~~^^
  File "/home/daq_user/PlasMEEP/scripts/circulator_inverse_design/multimode_optimization/refine_free91_multimode.py", line 88, in main
    J0, g0, met0 = _eval(s0)
                   ~~~~~^^^^
  File "/home/daq_user/PlasMEEP/scripts/circulator_inverse_design/multimode_optimization/refine_free91_multimode.py", line 55, in _eval
    bun, g, s_out = eval_and_grad(
                    ~~~~~~~~~~~~~^
        args.grade, s, orbits, F_HZ, args.B, with_absorption=False
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    )
    ^
  File "/home/daq_user/PlasMEEP/scripts/circulator_inverse_design/multimode_optimization/multimode_opt_core.py", line 319, in eval_and_grad
    gq = gradient_q_adjoint(bundle, orbits, o_masks, s)
  File "/home/daq_user/PlasMEEP/scripts/circulator_inverse_design/multimode_optimization/multimode_opt_core.py", line 266, in gradient_q_adjoint
    o_masks[m],
    ~~~~~~~^^^
IndexError: list index out of range

Elapsed run time = 413.0947 s
free91 skipped/failed
===== Final status =====
wrote outputs/circulator_inverse_design/multimode_optimization/FINAL_MULTIMODE_OPT_STATUS.md
{
  "MULTIMODE_C6_OPT": "PASS",
  "OPTIMIZED_FULL12_ONSAGER": "PASS",
  "OPTIMIZED_POWER_BALANCE": "FAIL",
  "MULTIMODE_FINE_MESH": "PASS",
  "FREE91_MULTIMODE": "NOT_ATTEMPTED",
  "BROADBAND_READY": "NEEDS_CUTOFF_DERIVATION"
}
===== SF CAMPAIGN COMPLETE 2026-10-16T01:23:44-07:00 =====
```
