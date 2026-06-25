# Module ④: when is a key-schedule modification cube-attackable?

The thesis attacks **one** modified ARADI key schedule. This module turns that single
example into a **characterisation**: which key-schedule modifications are practically and
uniquely cube-attackable, and why?

> **Scope:** reduced-round MODIFIED schemes only. Full ARADI is **not** attacked; recovery is
> of the seed, not the 256-bit master key (the thesis modification's seed→key map is the
> non-invertible projection `K_base = ⊕Kᵢ`).

## The principle

The attack brute-forces the seed half `l` over `2^h` candidates and filters them with an
AABB-style test of strength `f` bits, using `K` independent cubes. Two quantities govern it:

- **computational cost** ≈ `2^h` (brute force per half),
- **surviving candidates** ≈ `1 + 2^(h − f·K)` (true key + false positives).

So a modification is **weak** (practically + uniquely attackable) when it compresses the
relevant round-key entropy to a small half-width `h` that is brute-forceable **and** stays
below the filter strength `f·K` with a margin. It is **strong** when `h` stays large — the
original ARADI keeps the round-key pair at 192-bit entropy, so there is no feasible `h` and
the attack remains at ≈ `2^141`. The round function is unchanged throughout, so the AABB
filter is always available; only the key schedule varies.

## Measured (`python characterize.py`, 5 trials/point)

| h | f | K | avg survivors | idealised 1+2^(h−fK) | unique recovery | h ≤ f·K |
|---|---|---|---|---|---|---|
| 8 | 32 | 1 | 1.00 | 1.00 | 5/5 | yes |
| 12 | 32 | 1 | 1.00 | 1.00 | 5/5 | yes |
| **16** | **32** | **1** | **1.00** | **1.00** | **5/5** | **yes (thesis case)** |
| 8 | 16 | 1 | 1.00 | 1.00 | 5/5 | yes |
| 12 | 16 | 1 | 1.40 | 1.06 | 3/5 | yes (small margin) |
| 16 | 16 | 1 | 5.40 | 2.00 | 2/5 | boundary (margin 0) |
| 16 | 16 | 2 | 1.00 | 1.00 | 5/5 | yes |
| 8 | 8 | 1 | 3.40 | 2.00 | 0/5 | boundary |
| 12 | 8 | 1 | 22.80 | 17.00 | 0/5 | **NO** |
| 16 | 8 | 1 | 579.40 | 257.00 | 0/5 | **NO** |

**Reading.** The true half is *always* found. Unique recovery needs a margin `f·K − h > 0`:
survivors grow roughly as `2^(h−f·K)`, but the measured counts run a small factor (~2×) above
the idealised model — the byte-equality events are not fully independent — so a few bits of
margin are needed in practice. The thesis case (`h=16, f=32, K=1`) has a 16-bit margin, which
is why it recovers with **0 false positives, 100/100**. At the boundary uniqueness degrades;
beyond it false positives flood; multi-cube (`K=2`) widens the margin back. Full ARADI has no
feasible `h`, so the attack stays at ≈ `2^141`.

## Reproduce

```bash
python characterize.py    # imports the attack primitives from ../python/multicube_attack.py
```
