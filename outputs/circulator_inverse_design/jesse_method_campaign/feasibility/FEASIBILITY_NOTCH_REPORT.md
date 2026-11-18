# Notch character (dense audit, A2 ≤ −5 regions)

Stage-6: NOT STARTED

## A_4p9_5p1_C6_B0.0500_stage5
- band 4.9–5.1 GHz | B=0.05 | C6/16 | init=Stage-5 q_best
- audit worst A2=-7.5401 @ 5.0344 GHz; IL=28.11 dB; iso=-2.20 dB
- hidden vs opt-grid: gap=-4.1852 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **1** A2≤−5 region(s):
  - 5.0344–5.0344 GHz (~3.1 MHz, 1 pts); min A2=-7.5401 @ 5.0344 GHz

## A_4p9_5p1_C6_B0.0500_uniform
- band 4.9–5.1 GHz | B=0.05 | C6/16 | init=uniform q=1
- audit worst A2=-0.6416 @ 5.0906 GHz; IL=9.10 dB; iso=8.10 dB
- hidden vs opt-grid: gap=-0.0539 (audit−opt; negative ⇒ denser audit worse); hidden_failure=False
- **no** contiguous A2≤−5 regions on the dense audit grid

## B_4p75_5p25_C3_B0.0500_Bs5
- band 4.75–5.25 GHz | B=0.05 | C3/31 | init=file:B_4p75_5p25_C6_B0.0500_stage5_opt.json
- audit worst A2=-5.0208 @ 4.9094 GHz; IL=15.37 dB; iso=-3.07 dB
- hidden vs opt-grid: gap=-1.4677 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **1** A2≤−5 region(s):
  - 4.9094–4.9094 GHz (~3.1 MHz, 1 pts); min A2=-5.0208 @ 4.9094 GHz

## B_4p75_5p25_C3_B0.0500_stage5
- band 4.75–5.25 GHz | B=0.05 | C3/31 | init=Stage-5 q_best
- audit worst A2=-5.0008 @ 5.2266 GHz; IL=14.05 dB; iso=-5.45 dB
- hidden vs opt-grid: gap=-0.0327 (audit−opt; negative ⇒ denser audit worse); hidden_failure=False
- **1** A2≤−5 region(s):
  - 5.2266–5.2281 GHz (~3.1 MHz, 2 pts); min A2=-5.0008 @ 5.2266 GHz

## B_4p75_5p25_C6_B0.0500_Auni
- band 4.75–5.25 GHz | B=0.05 | C6/16 | init=file:A_4p9_5p1_C6_B0.0500_uniform_opt.json
- audit worst A2=-11.4619 @ 5.0328 GHz; IL=30.06 dB; iso=-16.11 dB
- hidden vs opt-grid: gap=-4.2027 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **5** A2≤−5 region(s):
  - 4.7516–4.7656 GHz (~15.6 MHz, 13 pts); min A2=-9.3293 @ 4.7578 GHz
  - 4.7844–4.7906 GHz (~7.8 MHz, 6 pts); min A2=-7.3550 @ 4.7891 GHz
  - 5.0312–5.0359 GHz (~6.2 MHz, 5 pts); min A2=-11.4619 @ 5.0328 GHz
  - 5.0578–5.0625 GHz (~6.2 MHz, 5 pts); min A2=-7.8892 @ 5.0609 GHz
  - 5.0750–5.1031 GHz (~29.7 MHz, 26 pts); min A2=-7.6659 @ 5.0844 GHz

## B_4p75_5p25_C6_B0.0500_stage5
- band 4.75–5.25 GHz | B=0.05 | C6/16 | init=Stage-5 q_best
- audit worst A2=-4.9639 @ 4.9094 GHz; IL=15.20 dB; iso=-3.02 dB
- hidden vs opt-grid: gap=-1.3925 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **no** contiguous A2≤−5 regions on the dense audit grid

## C_4p5_5p5_C6_B0.0500_stage5
- band 4.5–5.5 GHz | B=0.05 | C6/16 | init=Stage-5 q_best
- audit worst A2=-9.6490 @ 5.0844 GHz; IL=29.76 dB; iso=-10.18 dB
- hidden vs opt-grid: gap=-2.6497 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **8** A2≤−5 region(s):
  - 4.5625–4.5703 GHz (~9.4 MHz, 8 pts); min A2=-5.8991 @ 4.5672 GHz
  - 4.5734–4.5781 GHz (~6.2 MHz, 5 pts); min A2=-6.0208 @ 4.5734 GHz
  - 4.6813–4.6813 GHz (~3.1 MHz, 1 pts); min A2=-5.2868 @ 4.6813 GHz
  - 4.7359–4.7391 GHz (~4.7 MHz, 3 pts); min A2=-6.1782 @ 4.7375 GHz
  - 4.8031–4.8156 GHz (~14.1 MHz, 12 pts); min A2=-5.6692 @ 4.8063 GHz
  - 4.8438–4.8438 GHz (~3.1 MHz, 1 pts); min A2=-6.2791 @ 4.8438 GHz
  - 5.0750–5.0969 GHz (~23.4 MHz, 21 pts); min A2=-9.6490 @ 5.0844 GHz
  - 5.1156–5.1391 GHz (~25.0 MHz, 21 pts); min A2=-7.0408 @ 5.1297 GHz

## C_4p5_5p5_C6_B0.1000_Bs5
- band 4.5–5.5 GHz | B=0.1 | C6/16 | init=file:B_4p75_5p25_C6_B0.0500_stage5_opt.json
- audit worst A2=-5.3331 @ 4.6094 GHz; IL=14.16 dB; iso=-6.64 dB
- hidden vs opt-grid: gap=-0.0262 (audit−opt; negative ⇒ denser audit worse); hidden_failure=False
- **1** A2≤−5 region(s):
  - 4.6031–4.6250 GHz (~23.4 MHz, 21 pts); min A2=-5.3331 @ 4.6094 GHz

## C_4p5_5p5_C6_B0.1000_stage5
- band 4.5–5.5 GHz | B=0.1 | C6/16 | init=Stage-5 q_best
- audit worst A2=-8.7517 @ 4.5266 GHz; IL=25.49 dB; iso=-9.37 dB
- hidden vs opt-grid: gap=-4.0097 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **2** A2≤−5 region(s):
  - 4.5266–4.5297 GHz (~4.7 MHz, 4 pts); min A2=-8.7517 @ 4.5266 GHz
  - 4.5344–4.5344 GHz (~3.1 MHz, 1 pts); min A2=-5.3225 @ 4.5344 GHz

## D_4p25_5p75_C6_B0.0500_Bs5
- band 4.25–5.75 GHz | B=0.05 | C6/16 | init=file:B_4p75_5p25_C6_B0.0500_stage5_opt.json
- audit worst A2=-13.6664 @ 5.1641 GHz; IL=33.32 dB; iso=-23.80 dB
- hidden vs opt-grid: gap=-5.1493 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **5** A2≤−5 region(s):
  - 4.2516–4.2734 GHz (~23.4 MHz, 20 pts); min A2=-9.6831 @ 4.2594 GHz
  - 4.9625–4.9672 GHz (~6.2 MHz, 5 pts); min A2=-5.5310 @ 4.9672 GHz
  - 5.1594–5.1656 GHz (~7.8 MHz, 6 pts); min A2=-13.6664 @ 5.1641 GHz
  - 5.2281–5.2344 GHz (~7.8 MHz, 7 pts); min A2=-8.8154 @ 5.2328 GHz
  - 5.2438–5.2609 GHz (~18.7 MHz, 17 pts); min A2=-9.6698 @ 5.2469 GHz

## E_4p0_6p0_C6_B0.0500_stage5
- band 4.0–6.0 GHz | B=0.05 | C6/16 | init=Stage-5 campaign (production)
- audit worst A2=-8.7411 @ 4.7344 GHz; IL=25.04 dB; iso=-11.03 dB
- hidden vs opt-grid: gap=-1.1557 (audit−opt; negative ⇒ denser audit worse); hidden_failure=True
- **14** A2≤−5 region(s):
  - 4.1625–4.1719 GHz (~12.5 MHz, 4 pts); min A2=-6.2018 @ 4.1688 GHz
  - 4.1875–4.1906 GHz (~6.2 MHz, 2 pts); min A2=-7.2193 @ 4.1875 GHz
  - 4.2281–4.3125 GHz (~87.5 MHz, 28 pts); min A2=-8.2793 @ 4.2281 GHz
  - 4.3187–4.3656 GHz (~50.0 MHz, 16 pts); min A2=-7.6024 @ 4.3406 GHz
  - 4.5656–4.5812 GHz (~18.7 MHz, 6 pts); min A2=-6.6093 @ 4.5719 GHz
  - 4.7313–4.7500 GHz (~21.9 MHz, 7 pts); min A2=-8.7411 @ 4.7344 GHz
  - 4.8063–4.8125 GHz (~9.4 MHz, 3 pts); min A2=-5.4172 @ 4.8094 GHz
  - 4.8375–4.8375 GHz (~3.1 MHz, 1 pts); min A2=-5.8007 @ 4.8375 GHz
  - 5.0750–5.0969 GHz (~25.0 MHz, 8 pts); min A2=-8.2380 @ 5.0844 GHz
  - 5.1156–5.1437 GHz (~31.2 MHz, 10 pts); min A2=-7.6170 @ 5.1344 GHz
  - 5.5563–5.5781 GHz (~25.0 MHz, 8 pts); min A2=-8.0460 @ 5.5719 GHz
  - 5.5844–5.5906 GHz (~9.4 MHz, 3 pts); min A2=-6.5362 @ 5.5875 GHz
  - 5.7750–5.7812 GHz (~9.4 MHz, 3 pts); min A2=-6.5477 @ 5.7781 GHz
  - 5.8125–5.8125 GHz (~3.1 MHz, 1 pts); min A2=-5.1916 @ 5.8125 GHz

