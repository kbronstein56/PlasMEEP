#!/usr/bin/env python3
"""Deterministic Stage-5 frequency set generator.

Authorized rule from BB_STAGE5_SAMPLING_STRATEGY.md:
  uniform 6.25 MHz grid on [4.0, 6.0] GHz inclusive → n_freq = 321.

Do not modify this rule during Stage-5 optimization.
"""
from __future__ import annotations


def stage5_freqs_GHz(f0: float = 4.0, f1: float = 6.0, df: float = 0.00625) -> list[float]:
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def main():
    freqs = stage5_freqs_GHz()
    gaps = [freqs[i + 1] - freqs[i] for i in range(len(freqs) - 1)]
    print(f"n_freq = {len(freqs)}")
    print(f"f_min = {freqs[0]}  f_max = {freqs[-1]}")
    print(f"df_MHz = {min(gaps)*1000:.6f} .. {max(gaps)*1000:.6f}")
    print(f"max_gap_MHz = {max(gaps)*1000:.6f}")
    assert len(freqs) == 321
    assert abs(max(gaps) - 0.00625) < 1e-12
    assert freqs[0] == 4.0 and freqs[-1] == 6.0
    assert len(set(freqs)) == 321
    print("STAGE-5 FREQUENCY SET: AUTHORIZED (do not alter during opt)")


if __name__ == "__main__":
    main()
