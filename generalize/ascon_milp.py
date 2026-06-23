# -*- coding: utf-8 -*-
"""Milestone 3 of module 2: division-property degree UPPER bounds for Ascon.

A bit-based (conventional / two-subset) division-property MILP (PuLP+CBC) for the
Ascon permutation, using the same COPY / XOR / AND modelling rules as the ARADI
model in python/aradi_milp.py. It computes a SOUND (not necessarily tight) upper
bound on the algebraic degree of the r-round permutation, to pair with the lower
bounds from degree_growth.py.

The Ascon S-box is modelled directly from its algebraic normal form; the linear
layer is the per-word COPY+XOR of two rotations. Round constants are affine and
do not affect the degree, so they are skipped (standard convention).

Provenance of the parameters:
  - S-box LUT and linear-layer rotations: NIST SP 800-232 (final standard);
    the S-box LUT is its Table 6, the linear layer its rotation equations. The
    LUT here is also pinned in code: verify_anf() regenerates it from the
    canonical chi-based bit-sliced instruction sequence.
  - Algebraic-degree bound (cross-validation anchor): the Ascon v1.2 design
    specification (Dobraunig-Eichlseder-Mendel-Schlaeffer, CAESAR/NIST-LWC
    submission; J. Cryptology 34(3), 2021), Table 16: deg(p^R) <= 2,4,8,16,32,
    64,128,256,298,312,317,319 for R=1..12 (so the bound equals 2^R only up to
    R=8). NIST SP 800-232 itself contains no algebraic-degree analysis.

We assert the MILP reproduces the upper side 2,4,8 for R=1,2,3, AND a
discrimination test (a linear S-box must give degree 1; an injected degree-3
monomial must give degree 3) so the 2^R match is meaningful, not a vacuous
over-count.

Honest scope: structural degree measurement, NOT an attack on Ascon.
"""
from __future__ import annotations

import time

import pulp

# ---------------------------------------------------------------------------
# Ascon S-box: LUT (NIST SP 800-232 Table 6) and ANF (eq. (6)), x0 = MSB.
# ---------------------------------------------------------------------------
ASCON_SBOX = [0x4, 0xb, 0x1f, 0x14, 0x1a, 0x15, 0x9, 0x2, 0x1b, 0x5, 0x8, 0x12,
              0x1d, 0x3, 0x6, 0x1c, 0x1e, 0x13, 0x7, 0xe, 0x0, 0xd, 0x11, 0x18,
              0x10, 0xc, 0x1, 0x19, 0x16, 0xa, 0xf, 0x17]

# ANF of each output bit y0..y4 as a list of monomials (tuples of input indices
# into x0..x4); y2 also has a constant 1 (affine, ignored by division property
# but kept so the ANF can be self-checked against the LUT).
ASCON_ANF = {
    0: [(0, 1), (1, 2), (1, 4), (0,), (1,), (2,), (3,)],
    1: [(1, 2), (1, 3), (2, 3), (0,), (1,), (2,), (3,), (4,)],
    2: [(3, 4), (1,), (2,), (4,)],                      # + constant 1
    3: [(0, 3), (0, 4), (0,), (1,), (2,), (3,), (4,)],
    4: [(0, 1), (1, 4), (1,), (3,), (4,)],
}
ASCON_CONST = {2: 1}  # y2 has the affine +1

ROT = [(19, 28), (61, 39), (1, 6), (10, 17), (7, 41)]
WORD_BITS = 64
NWORDS = 5


def _sbox_from_anf(v: int) -> int:
    """Evaluate the ANF on input integer v (x0=MSB) -> output integer (y0=MSB)."""
    x = [(v >> (4 - i)) & 1 for i in range(5)]
    out = 0
    for j in range(5):
        bit = ASCON_CONST.get(j, 0)
        for mon in ASCON_ANF[j]:
            term = 1
            for idx in mon:
                term &= x[idx]
            bit ^= term
        out |= bit << (4 - j)
    return out


def _sbox_from_instructions(v: int) -> int:
    """Canonical Ascon S-box from the bit-sliced chi instruction sequence
    (independent of the hardcoded LUT) -> pins the table in code."""
    n = lambda b: b ^ 1
    x0, x1, x2, x3, x4 = (v >> 4) & 1, (v >> 3) & 1, (v >> 2) & 1, (v >> 1) & 1, v & 1
    x0 ^= x4; x4 ^= x3; x2 ^= x1
    t0 = n(x0) & x1; t1 = n(x1) & x2; t2 = n(x2) & x3; t3 = n(x3) & x4; t4 = n(x4) & x0
    x0 ^= t1; x1 ^= t2; x2 ^= t3; x3 ^= t4; x4 ^= t0
    x1 ^= x0; x0 ^= x4; x3 ^= x2; x2 = n(x2)
    return (x0 << 4) | (x1 << 3) | (x2 << 2) | (x3 << 1) | x4


def verify_anf() -> None:
    """Pin the S-box two independent ways: (1) the LUT must equal the canonical
    chi instruction sequence, (2) the ANF must reproduce that LUT."""
    assert [_sbox_from_instructions(v) for v in range(32)] == ASCON_SBOX, \
        "ASCON LUT does not match the canonical chi instruction sequence"
    assert [_sbox_from_anf(v) for v in range(32)] == ASCON_SBOX, \
        "ASCON ANF does not reproduce the S-box LUT"


# ---------------------------------------------------------------------------
# MILP primitives (same rules as python/aradi_milp.py)
# ---------------------------------------------------------------------------
_CNT = {"n": 0}


def _bvar() -> pulp.LpVariable:
    _CNT["n"] += 1
    return pulp.LpVariable(f"v{_CNT['n']}", cat=pulp.LpBinary)


def _copy(model, a, outs):
    for b in outs:
        model += a >= b
    model += pulp.lpSum(outs) >= a


def _xor(model, ins, out):
    model += out == pulp.lpSum(ins)


def _and(model, ins, out):
    for b in ins:
        model += out == b


# ---------------------------------------------------------------------------
# S-box layer (from ANF): COPY each input to its uses, AND for quadratic
# monomials, XOR monomial-copies into each output bit. `anf` is a parameter so
# the discrimination test can substitute a linear or higher-degree S-box.
# ---------------------------------------------------------------------------
def model_sbox_column(model, x_in, anf=ASCON_ANF):
    """x_in: 5 binary vars (x0..x4) -> 5 output vars (y0..y4), division property."""
    uses = [0] * 5
    for j in range(5):
        for mon in anf[j]:
            for idx in mon:
                uses[idx] += 1
    pools = []
    for idx in range(5):
        n = uses[idx]
        if n == 0:
            pools.append([])
        elif n == 1:
            pools.append([x_in[idx]])
        else:
            cs = [_bvar() for _ in range(n)]
            _copy(model, x_in[idx], cs)
            pools.append(cs)
    ptr = [0] * 5

    def take(idx):
        c = pools[idx][ptr[idx]]
        ptr[idx] += 1
        return c

    y_out = []
    for j in range(5):
        term_vars = []
        for mon in anf[j]:
            if len(mon) == 1:
                term_vars.append(take(mon[0]))
            else:
                copies = [take(idx) for idx in mon]
                t = _bvar()
                _and(model, copies, t)
                term_vars.append(t)
        yj = _bvar()
        _xor(model, term_vars, yj)
        y_out.append(yj)
    return y_out


def model_sbox_layer(model, state):
    out = [[None] * WORD_BITS for _ in range(NWORDS)]
    for p in range(WORD_BITS):
        col_out = model_sbox_column(model, [state[w][p] for w in range(NWORDS)])
        for w in range(NWORDS):
            out[w][p] = col_out[w]
    return out


# ---------------------------------------------------------------------------
# Linear layer: per word, y[p] = x[p] ^ x[p+r1] ^ x[p+r2] (right rotation).
# Each input bit feeds 3 output positions -> COPY into 3; each output = XOR of 3.
# ---------------------------------------------------------------------------
def model_linear_layer(model, state):
    out = []
    for w in range(NWORDS):
        r1, r2 = ROT[w]
        contributions = [[] for _ in range(WORD_BITS)]
        for p in range(WORD_BITS):
            cs = [_bvar() for _ in range(3)]
            _copy(model, state[w][p], cs)
            contributions[p].append(cs[0])
            contributions[(p + r1) % WORD_BITS].append(cs[1])
            contributions[(p + r2) % WORD_BITS].append(cs[2])
        out.append([_xor_new(model, contributions[q]) for q in range(WORD_BITS)])
    return out


def _xor_new(model, ins):
    o = _bvar()
    _xor(model, ins, o)
    return o


# ---------------------------------------------------------------------------
# Degree upper bound for r rounds (all input bits free = full-permutation degree).
# One target output bit suffices: the round function is identical on every bit
# column and the linear layer is a circular rotation, so by symmetry the degree
# bound is the same for each output bit of a word; we target bit 0 of word 0.
# ---------------------------------------------------------------------------
def degree_upper_bound(rounds: int, target_word: int = 0, target_bit: int = 0,
                       time_limit=None, msg=False):
    _CNT["n"] = 0
    model = pulp.LpProblem("ASCON_BDP", pulp.LpMaximize)
    state = [[_bvar() for _ in range(WORD_BITS)] for _ in range(NWORDS)]
    inputs = [state[w][p] for w in range(NWORDS) for p in range(WORD_BITS)]
    for _ in range(rounds):
        state = model_sbox_layer(model, state)   # round constant skipped (affine)
        state = model_linear_layer(model, state)
    for w in range(NWORDS):
        for p in range(WORD_BITS):
            model += state[w][p] == (1 if (w == target_word and p == target_bit) else 0)
    model += pulp.lpSum(inputs)
    t = time.time()
    kw = {"msg": msg}
    if time_limit:
        kw["timeLimit"] = time_limit
    status = model.solve(pulp.PULP_CBC_CMD(**kw))
    el = time.time() - t
    if pulp.LpStatus[status] == "Infeasible":
        return -1, el
    val = pulp.value(model.objective)
    return (int(round(val)) if val is not None else -2), el


# ---------------------------------------------------------------------------
# Discrimination / negative control: the S-box model must report the TRUE degree,
# not just over-count to the trivial bound. Test on a single S-box column.
# ---------------------------------------------------------------------------
def _sbox_column_degree(anf):
    best = 0
    for tb in range(5):
        _CNT["n"] = 0
        m = pulp.LpProblem("col", pulp.LpMaximize)
        xin = [_bvar() for _ in range(5)]
        yout = model_sbox_column(m, xin, anf)
        for j in range(5):
            m += yout[j] == (1 if j == tb else 0)
        m += pulp.lpSum(xin)
        m.solve(pulp.PULP_CBC_CMD(msg=False))
        v = pulp.value(m.objective)
        if v is not None:
            best = max(best, int(round(v)))
    return best


def discrimination_test():
    real = _sbox_column_degree(ASCON_ANF)                    # expect 2
    linear = _sbox_column_degree({j: [(j,)] for j in range(5)})  # expect 1
    deg3_anf = dict(ASCON_ANF)
    deg3_anf[0] = ASCON_ANF[0] + [(0, 1, 2)]                 # inject a cubic monomial
    deg3 = _sbox_column_degree(deg3_anf)                      # expect 3
    assert (real, linear, deg3) == (2, 1, 3), \
        f"discrimination failed: real={real} linear={linear} deg3={deg3}"
    return real, linear, deg3


def main():
    verify_anf()
    print("Self-check: LUT == chi instruction sequence, and ANF reproduces LUT  [OK]")
    real, linear, deg3 = discrimination_test()
    print(f"Discrimination: S-box column degree -> real={real}, linear={linear}, "
          f"degree-3-injected={deg3}  [OK] (model is not a vacuous over-count)")
    print("\nAscon permutation -- algebraic-degree UPPER bound via division property")
    print("(round constants skipped: affine, degree-neutral)")
    print(f"  {'rounds':>6} | {'deg upper bound':>16} | {'2^r (spec upper)':>16} | {'time s':>7}")
    print("  " + "-" * 60)
    for r in (1, 2, 3):
        ub, el = degree_upper_bound(r, time_limit=120)
        exp = 2 ** r
        tag = "OK" if ub == exp else "MISMATCH"
        print(f"  {r:>6} | {ub:>16} | {exp:>16} | {el:>7.1f}  {tag}")
        assert ub == exp, f"round {r}: MILP upper bound {ub} != 2^{r} = {exp}"
    print("\nMILP reproduces the upper side 2,4,8 for R=1..3 (Ascon v1.2 spec Table 16:")
    print("deg(p^R) <= 2^R for R<=8). The discrimination test above shows the model")
    print("detects degrees BELOW the trivial bound, so the match is meaningful.")
    print("NOT an attack on Ascon.")


if __name__ == "__main__":
    main()
