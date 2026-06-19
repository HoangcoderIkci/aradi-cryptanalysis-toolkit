# ARADI cryptanalysis toolkit

Reproducible tooling for cube cryptanalysis of the **ARADI** low-latency block cipher
(NSA, 2024): a verified C reference implementation, an MILP verification of the **AABB
cube-distinguisher**, and a practical cube attack on a **modified 6-round key schedule**.

> ⚠️ **Scope:** this does **not** break full ARADI and does **not** recover the 256-bit
> master key. Full 16-round ARADI stays at ≈ 2¹⁴¹ (practically infeasible). The practical
> attack targets a **modified, reduced 6-round** key schedule and recovers only a 32-bit
> pair `(u, ℓ)`. See [`SCOPE.md`](SCOPE.md) — please read it before citing any result.

---

## Tóm tắt (Tiếng Việt)

Bộ công cụ tái lập được cho **kubическая (cube) cryptanalysis** của ARADI:
- **`impl/`** — cài đặt ARADI bằng C, **khớp test-vector chính thức của NSA**.
- **`python/`** — kiểm chứng bằng **MILP (PuLP+CBC)** tính chất **AABB cube-distinguisher**
  (xuất hiện sau **5 vòng**), và **tấn công cube thực nghiệm** lên **biến thể 6 vòng đã sửa
  key schedule** (phục hồi cặp 32-bit `(u, ℓ)` trong vài giây, **100/100** master key ngẫu nhiên).
- Phạm vi trung thực: **không phá** ARADI đầy đủ (~2¹⁴¹), **không** lấy master key 256-bit
  (projection `K_base = ⊕Kᵢ` không khả nghịch). Xem [`SCOPE.md`](SCOPE.md).

---

## Headline result (modified 6-round scheme)

| Metric | Value |
|---|---|
| Target | **modified** 6-round ARADI key schedule (not full 16-round) |
| Recovered | 32-bit pair `(u, ℓ)` — **not** the 256-bit master key |
| Success rate | **100/100** random master keys (thesis original: 20/20) |
| False positives (Mode A, full AABB) | **0** |
| Time (Phase 2, Python+NumPy, single thread) | ~ a few seconds / key |
| Full 16-round ARADI | analysis gives ≈ **2¹⁴¹** → practically infeasible (not broken) |

## Quickstart

```bash
make test-c      # build the C impl and check it against the NSA test vector (exit 0 = PASS)
make reproduce   # run the headline experiment: 100 random keys on modified 6-round ARADI
make milp        # MILP verification of the AABB cube-distinguisher
```
Requirements: a C compiler (`gcc`/`clang`), Python 3.10+, `numpy`, `pulp` (bundles CBC).
Install Python deps with `pip install -r requirements.txt`.

## Repository layout

```
impl/         C reference implementation (aradi.c/.h) + NSA test vector check
python/        Python reference + MILP model + cube attack + experiment runners
  aradi_ref.py          reference impl used by the attack
  aradi_milp.py         MILP model of the cube-distinguisher (PuLP+CBC)
  multicube_attack.py   the cube attack + multi-cube false-positive optimisation
  run_multicube_100.py  headline experiment (100 keys, prints HEADLINE_VERDICT)
  run_6round.py / run_5round.py / run_classification.py / run_experimental.py
generalize/    [roadmap] generalising AABB/cube to the Toffoli/χ cipher family — see README
ml/            [roadmap] neural distinguisher vs the AABB distinguisher — see README
SCOPE.md       what this work does and does NOT claim (read first)
```

## Roadmap (work in progress)

- **`generalize/`** — is the AABB property a *structural* feature of Toffoli-gate / χ-based
  S-boxes (Simon, Xoodoo, Ascon/Keccak-χ)? Turn a single-cipher result into a class insight.
- **`ml/`** — train a Gohr-style neural distinguisher on round-reduced ARADI and compare it
  head-to-head with the analytic AABB distinguisher on the same testbed.

## Background

Distilled from a B.Sc. thesis (2026) in applied cryptography, *"Analysis of a modification
of the ARADI algorithm"* (cube cryptanalysis). Builds on the public ARADI analyses by
Bellini et al. (ePrint 2024/1559), Kim et al. (2024/1772) and the NSA specification
(2024/1240).

## License & citation

MIT (see [`LICENSE`](LICENSE)). If you use this, please cite via [`CITATION.cff`](CITATION.cff).
