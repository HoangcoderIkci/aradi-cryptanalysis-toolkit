"""
MILP-based three-subset bit-based division property (3SBDP) model for the
ARADI block cipher.

The implementation follows the modelling rules of Hao et al. (EUROCRYPT 2020)
and Algorithms 1-2 of Kim et al. "Byte-wise equal property of ARADI"
(IACR ePrint 2024/1772).

State layout (matches main.tex / paper notation):
    A round state is the 128-bit value (W, X, Y, Z) where W is most significant
    and Z is least significant.  Inside each 32-bit word, bit index 0 in this
    file (the "logical" index used everywhere) refers to the SAME bit as
    Kim et al.'s W^r_i, X^r_i, Y^r_i, Z^r_i, i.e. the bits inside a word are
    numbered 0..31 from the *most* significant end (W_0 = MSB of W).

Bitwise modelling rules used (Hao et al. EUROCRYPT 2020 / Bellini et al. 2024):
    * XOR    : (b_1, ..., b_n) -> a   ==>   a = b_1 + ... + b_n
    * AND    : (b_1, ..., b_n) -> a   ==>   a = b_i for every i
                                            (correct for three-subset rule)
    * COPY   : a -> (b_1, ..., b_n)   ==>   a >= b_i for every i,
                                            sum(b_i) >= a
    * XOR with constant 1 : a -> b    ==>   b >= a

Solver: PuLP with CBC (free, ships with PuLP).

Public API:
    build_aradi_milp(num_rounds, cube_indices_w/_x/_y/_z, target_bit, target_word)
    compute_degree_upper_bound(num_rounds, cube_indices_w/_x/_y/_z,
                               target_bit, target_word, target_funcs=None,
                               time_limit=None)

The cube-index sets are sets of integers in [0, 31].  `target_word` is one of
'W','X','Y','Z'.  `target_funcs` lets the caller stack several output bits via
a single XOR (used for X^r ⊕ Z^r in Kim et al.).
"""

from __future__ import annotations

import time
from typing import Iterable, List, Optional, Sequence, Tuple

import pulp

# -----------------------------------------------------------------------------
# ARADI parameters
# -----------------------------------------------------------------------------

LINEAR_PARAMS = [(11, 8, 14), (10, 9, 11), (9, 4, 14), (8, 9, 7)]

WORD_BITS = 32
HALF_BITS = 16
STATE_BITS = 128

# In every word we number "logical bit 0" as the MSB (matches the indexing
# used in main.tex and in Kim et al.).  Conversions to/from the "u" (high half)
# and "l" (low half) are therefore:
#     u_j = bit j     (0 <= j < 16)
#     l_j = bit 16+j  (0 <= j < 16)


# -----------------------------------------------------------------------------
# Variable factory
# -----------------------------------------------------------------------------

_VAR_COUNTER = {"n": 0}


def _fresh_name(prefix: str) -> str:
    _VAR_COUNTER["n"] += 1
    return f"{prefix}_{_VAR_COUNTER['n']}"


def _binary_var(model: pulp.LpProblem, name: Optional[str] = None) -> pulp.LpVariable:
    if name is None:
        name = _fresh_name("v")
    v = pulp.LpVariable(name, cat=pulp.LpBinary)
    return v


def _binary_vars(model: pulp.LpProblem, n: int,
                 prefix: str = "v") -> List[pulp.LpVariable]:
    return [_binary_var(model, f"{prefix}_{_VAR_COUNTER['n']:d}_{i}")
            for i in range(n)]


# -----------------------------------------------------------------------------
# Basic 3SBDP modelling primitives
# -----------------------------------------------------------------------------

def model_copy(model: pulp.LpProblem,
               a: pulp.LpVariable,
               outs: Sequence[pulp.LpVariable]) -> None:
    """COPY rule: a -> (b_1, ..., b_n).  a >= b_i for every i; sum(b_i) >= a."""
    for b in outs:
        model += a >= b
    model += pulp.lpSum(outs) >= a


def model_xor(model: pulp.LpProblem,
              ins: Sequence[pulp.LpVariable],
              out: pulp.LpVariable) -> None:
    """XOR rule: a = b_1 + ... + b_n."""
    model += out == pulp.lpSum(ins)


def model_xor_const1(model: pulp.LpProblem,
                     a: pulp.LpVariable,
                     b: pulp.LpVariable) -> None:
    """XOR with the constant 1: b >= a."""
    model += b >= a


def model_and(model: pulp.LpProblem,
              ins: Sequence[pulp.LpVariable],
              out: pulp.LpVariable) -> None:
    """Bitwise AND rule: a = b_i for every i."""
    for b in ins:
        model += out == b


# -----------------------------------------------------------------------------
# Toffoli gate (used inside the S-box)
# -----------------------------------------------------------------------------

def model_toffoli(model: pulp.LpProblem,
                  a_in: pulp.LpVariable,
                  b_in: pulp.LpVariable,
                  c_in: pulp.LpVariable
                  ) -> Tuple[pulp.LpVariable, pulp.LpVariable, pulp.LpVariable]:
    """
    Single Toffoli gate (a, b, c) -> (a, b, c XOR a AND b).
    Returns the three output variables.

    Implementation:
        * COPY a_in into two copies (a_out, a_for_and)
        * COPY b_in into two copies (b_out, b_for_and)
        * AND  (a_for_and, b_for_and) -> t   (three-subset AND)
        * XOR  (c_in, t) -> c_out
    """
    a_out = _binary_var(model, _fresh_name("tof_a"))
    a_for_and = _binary_var(model, _fresh_name("tof_aA"))
    b_out = _binary_var(model, _fresh_name("tof_b"))
    b_for_and = _binary_var(model, _fresh_name("tof_bA"))
    t = _binary_var(model, _fresh_name("tof_t"))
    c_out = _binary_var(model, _fresh_name("tof_c"))

    model_copy(model, a_in, [a_out, a_for_and])
    model_copy(model, b_in, [b_out, b_for_and])
    model_and(model, [a_for_and, b_for_and], t)
    model_xor(model, [c_in, t], c_out)
    return a_out, b_out, c_out


# -----------------------------------------------------------------------------
# S-box layer:  4 Toffoli gates over (w, x, y, z), 32 lanes
# -----------------------------------------------------------------------------

def model_sbox(model: pulp.LpProblem,
               w_in: List[pulp.LpVariable], x_in: List[pulp.LpVariable],
               y_in: List[pulp.LpVariable], z_in: List[pulp.LpVariable]
               ) -> Tuple[List[pulp.LpVariable], List[pulp.LpVariable],
                          List[pulp.LpVariable], List[pulp.LpVariable]]:
    """
    Model the ARADI S-box layer (4 Toffoli gates, applied bit-by-bit).

    Order from the NSA spec / main.tex:
        x ^= w & y      ==>  Tof(w, y, x)
        z ^= x & y      ==>  Tof(x, y, z)
        y ^= w & z      ==>  Tof(w, z, y)
        w ^= x & z      ==>  Tof(x, z, w)

    Inputs are length-32 lists of binary MILP variables.  Outputs are again
    length-32 lists of fresh binary variables.
    """
    assert len(w_in) == 32 == len(x_in) == len(y_in) == len(z_in)
    w_out, x_out, y_out, z_out = [], [], [], []
    for bit in range(WORD_BITS):
        w, y, x = model_toffoli(model, w_in[bit], y_in[bit], x_in[bit])
        x, y, z = model_toffoli(model, x, y, z_in[bit])
        w, z, y = model_toffoli(model, w, z, y)
        x, z, w = model_toffoli(model, x, z, w)
        w_out.append(w); x_out.append(x); y_out.append(y); z_out.append(z)
    return w_out, x_out, y_out, z_out


# -----------------------------------------------------------------------------
# Linear layer L_i for one 32-bit word
# -----------------------------------------------------------------------------

def _u_index(j: int) -> int:
    """Logical bit index of bit j (LSB=0) of the upper 16-bit half u."""
    return HALF_BITS + j  # bits 16..31


def _l_index(j: int) -> int:
    """Logical bit index of bit j (LSB=0) of the lower 16-bit half l."""
    return j  # bits 0..15


def model_linear_word(model: pulp.LpProblem,
                      in_word: List[pulp.LpVariable],
                      round_idx: int
                      ) -> List[pulp.LpVariable]:
    """
    Model L_i acting on one 32-bit word.

    L_i is linear, computed as
        (u, l) -> (u XOR rotl16(u, a) XOR rotl16(l, c),
                   l XOR rotl16(l, a) XOR rotl16(u, b)).

    For division-property propagation we have to COPY each input bit into
    however many output bits it feeds.  In this linear layer every u-bit feeds
    exactly three output bits (one u-output via the identity term, one u-output
    via the rotl(u, a) term, one l-output via the rotl(u, b) term), and the
    same for l-bits.  Each output bit is then the XOR of exactly three
    incoming copies.
    """
    a, b, c = LINEAR_PARAMS[round_idx & 3]

    # Gather, for every input bit, the list of "output bit slots" it
    # contributes to.  We then COPY the input bit into len(slots) copies and
    # XOR three such copies into every output bit.

    # contributions[out_bit_index] = list of variables to XOR
    contributions: List[List[pulp.LpVariable]] = [[] for _ in range(WORD_BITS)]

    def alloc_copies(src: pulp.LpVariable, n: int) -> List[pulp.LpVariable]:
        if n == 1:
            return [src]
        outs = [_binary_var(model, _fresh_name("lin_c")) for _ in range(n)]
        model_copy(model, src, outs)
        return outs

    # ----- contributions from u-bits -----
    for j in range(HALF_BITS):
        u_src = in_word[_u_index(j)]
        # u-output identity: out u-bit at position j
        out_id = _u_index(j)
        # rotl(u, a): out u-bit at position (j + a) mod 16
        out_a = _u_index((j + a) % HALF_BITS)
        # rotl(u, b): out l-bit at position (j + b) mod 16
        out_b = _l_index((j + b) % HALF_BITS)

        copies = alloc_copies(u_src, 3)
        contributions[out_id].append(copies[0])
        contributions[out_a].append(copies[1])
        contributions[out_b].append(copies[2])

    # ----- contributions from l-bits -----
    for j in range(HALF_BITS):
        l_src = in_word[_l_index(j)]
        # l-output identity: out l-bit at position j
        out_id = _l_index(j)
        # rotl(l, a): out l-bit at position (j + a) mod 16
        out_a = _l_index((j + a) % HALF_BITS)
        # rotl(l, c): out u-bit at position (j + c) mod 16
        out_c = _u_index((j + c) % HALF_BITS)

        copies = alloc_copies(l_src, 3)
        contributions[out_id].append(copies[0])
        contributions[out_a].append(copies[1])
        contributions[out_c].append(copies[2])

    # XOR every contribution list into a fresh output bit.
    out_word: List[pulp.LpVariable] = []
    for ob in range(WORD_BITS):
        out_var = _binary_var(model, _fresh_name("lin_o"))
        model_xor(model, contributions[ob], out_var)
        out_word.append(out_var)

    return out_word


def model_linear(model: pulp.LpProblem,
                 w_in: List[pulp.LpVariable], x_in: List[pulp.LpVariable],
                 y_in: List[pulp.LpVariable], z_in: List[pulp.LpVariable],
                 round_idx: int
                 ) -> Tuple[List[pulp.LpVariable], List[pulp.LpVariable],
                            List[pulp.LpVariable], List[pulp.LpVariable]]:
    """Apply L_i to each of W, X, Y, Z (32-bit lanes)."""
    w_out = model_linear_word(model, w_in, round_idx)
    x_out = model_linear_word(model, x_in, round_idx)
    y_out = model_linear_word(model, y_in, round_idx)
    z_out = model_linear_word(model, z_in, round_idx)
    return w_out, x_out, y_out, z_out


# -----------------------------------------------------------------------------
# Round-key XOR
# -----------------------------------------------------------------------------

def model_roundkey_xor(model: pulp.LpProblem,
                       state_in: List[pulp.LpVariable],
                       key_vars: List[pulp.LpVariable]
                       ) -> List[pulp.LpVariable]:
    """
    Model XOR of a state half (or full state) with a fresh key half.

    The cube-attack division-property convention is to treat each round-key
    bit as a free input variable (its own fresh COPY).  Concretely, for each
    bit we have:
        state_in[i] (data side, may be zero) and key_vars[i] (key side)
        =>  out[i] = state_in[i] + key_vars[i]
    This is the standard XOR rule.
    """
    assert len(state_in) == len(key_vars)
    out = []
    for s_in, k_in in zip(state_in, key_vars):
        o = _binary_var(model, _fresh_name("rk"))
        model_xor(model, [s_in, k_in], o)
        out.append(o)
    return out


# -----------------------------------------------------------------------------
# Full ARADI MILP model
# -----------------------------------------------------------------------------

def _set_initial_division_property(model: pulp.LpProblem,
                                   w0: List[pulp.LpVariable],
                                   x0: List[pulp.LpVariable],
                                   y0: List[pulp.LpVariable],
                                   z0: List[pulp.LpVariable],
                                   cube_w: Iterable[int],
                                   cube_x: Iterable[int],
                                   cube_y: Iterable[int],
                                   cube_z: Iterable[int]) -> None:
    """
    Fix the initial 3SBDP vector k.  Cube bits are unrestricted (the solver
    sees them as binary), non-cube bits are forced to zero (constants).

    The MILP feasibility / objective is then taken over all valid trails.
    """
    cube = {
        'W': set(cube_w), 'X': set(cube_x),
        'Y': set(cube_y), 'Z': set(cube_z),
    }
    for word_letter, var_list in zip("WXYZ", (w0, x0, y0, z0)):
        for bit in range(WORD_BITS):
            if bit not in cube[word_letter]:
                model += var_list[bit] == 0


def build_aradi_milp(num_rounds: int,
                     cube_indices_w: Iterable[int] = (),
                     cube_indices_x: Iterable[int] = (),
                     cube_indices_y: Iterable[int] = (),
                     cube_indices_z: Iterable[int] = (),
                     target_bit_position: int = 0,
                     target_word: str = 'X',
                     target_funcs: Optional[Sequence[Tuple[str, int]]] = None,
                     ) -> Tuple[pulp.LpProblem, dict]:
    """
    Build a full r-round ARADI 3SBDP MILP model.

    Parameters
    ----------
    num_rounds : int
        Number of rounds r.
    cube_indices_{w,x,y,z} : iterable of int (in 0..31)
        The cube-index sets I_W, I_X, I_Y, I_Z.
    target_bit_position : int
        Bit index (0..31) of the *single* target bit, when target_funcs is None.
    target_word : 'W'|'X'|'Y'|'Z'
        Which output word the target bit lives in.
    target_funcs : optional list of (word_letter, bit_index)
        If supplied, the target is the XOR of these output bits.  Used to
        model X^r XOR Z^r.

    Returns
    -------
    (model, ctx)
        ctx contains 'cube_vars' (the cube-side input bits whose sum is
        maximised) and 'state_out' (the four output word variables).
    """
    global _VAR_COUNTER
    _VAR_COUNTER["n"] = 0

    model = pulp.LpProblem("ARADI_3SBDP", pulp.LpMaximize)

    # Step 1: create the initial 128 data-side bits w0/x0/y0/z0.
    w_state = _binary_vars(model, WORD_BITS, "W0")
    x_state = _binary_vars(model, WORD_BITS, "X0")
    y_state = _binary_vars(model, WORD_BITS, "Y0")
    z_state = _binary_vars(model, WORD_BITS, "Z0")

    initial_state = {'W': w_state, 'X': x_state,
                     'Y': y_state, 'Z': z_state}

    _set_initial_division_property(
        model, w_state, x_state, y_state, z_state,
        cube_indices_w, cube_indices_x, cube_indices_y, cube_indices_z)

    # Step 2: r rounds of (KeyXOR, S-box, Linear)
    for j in range(num_rounds):
        # Fresh 128-bit round-key variables.
        rk_w = _binary_vars(model, WORD_BITS, f"rkW{j}")
        rk_x = _binary_vars(model, WORD_BITS, f"rkX{j}")
        rk_y = _binary_vars(model, WORD_BITS, f"rkY{j}")
        rk_z = _binary_vars(model, WORD_BITS, f"rkZ{j}")

        w_state = model_roundkey_xor(model, w_state, rk_w)
        x_state = model_roundkey_xor(model, x_state, rk_x)
        y_state = model_roundkey_xor(model, y_state, rk_y)
        z_state = model_roundkey_xor(model, z_state, rk_z)

        w_state, x_state, y_state, z_state = model_sbox(
            model, w_state, x_state, y_state, z_state)
        w_state, x_state, y_state, z_state = model_linear(
            model, w_state, x_state, y_state, z_state, j)

    # Step 3: final post-whitening key XOR.
    rk_w = _binary_vars(model, WORD_BITS, f"rkW{num_rounds}")
    rk_x = _binary_vars(model, WORD_BITS, f"rkX{num_rounds}")
    rk_y = _binary_vars(model, WORD_BITS, f"rkY{num_rounds}")
    rk_z = _binary_vars(model, WORD_BITS, f"rkZ{num_rounds}")
    w_state = model_roundkey_xor(model, w_state, rk_w)
    x_state = model_roundkey_xor(model, x_state, rk_x)
    y_state = model_roundkey_xor(model, y_state, rk_y)
    z_state = model_roundkey_xor(model, z_state, rk_z)

    state_out = {'W': w_state, 'X': x_state, 'Y': y_state, 'Z': z_state}

    # Step 4: Set the target bit(s) to 1, all the others to 0.
    if target_funcs is None:
        target_funcs = [(target_word, target_bit_position)]
    target_set = set(target_funcs)
    for letter in "WXYZ":
        for bit in range(WORD_BITS):
            if (letter, bit) in target_set:
                model += state_out[letter][bit] == 1
            else:
                model += state_out[letter][bit] == 0

    # Step 5: Objective -- maximise number of active cube-side input bits.
    # The cube-side input bits are exactly the *initial* state bits whose
    # index belongs to one of the four cube-index sets.
    cube_input_vars = []
    for letter, idx_set in zip("WXYZ", (cube_indices_w, cube_indices_x,
                                         cube_indices_y, cube_indices_z)):
        for bit in idx_set:
            cube_input_vars.append(initial_state[letter][bit])

    model += pulp.lpSum(cube_input_vars)

    ctx = {
        'cube_vars': cube_input_vars,
        'num_vars': sum(1 for _ in model.variables()),
    }
    return model, ctx


# -----------------------------------------------------------------------------
# Solving
# -----------------------------------------------------------------------------

def compute_degree_upper_bound(num_rounds: int,
                               cube_indices_w: Iterable[int] = (),
                               cube_indices_x: Iterable[int] = (),
                               cube_indices_y: Iterable[int] = (),
                               cube_indices_z: Iterable[int] = (),
                               target_bit_position: int = 0,
                               target_word: str = 'X',
                               target_funcs: Optional[Sequence[Tuple[str, int]]] = None,
                               time_limit: Optional[int] = None,
                               msg: bool = False
                               ) -> Tuple[int, float, int]:
    """
    Build the model, solve it with CBC, and return
        (degree_upper_bound, solve_time_seconds, num_vars).

    If the MILP is infeasible the routine returns -1, signifying that no valid
    division trail reaches the chosen target bit (the bit is *constant zero*
    on the cube).
    """
    t_start = time.time()
    model, ctx = build_aradi_milp(
        num_rounds=num_rounds,
        cube_indices_w=cube_indices_w,
        cube_indices_x=cube_indices_x,
        cube_indices_y=cube_indices_y,
        cube_indices_z=cube_indices_z,
        target_bit_position=target_bit_position,
        target_word=target_word,
        target_funcs=target_funcs,
    )
    solver_kwargs = {"msg": msg}
    if time_limit is not None:
        solver_kwargs["timeLimit"] = time_limit
    solver = pulp.PULP_CBC_CMD(**solver_kwargs)
    status = model.solve(solver)
    elapsed = time.time() - t_start

    status_str = pulp.LpStatus[status]
    if status_str == "Infeasible":
        return -1, elapsed, ctx['num_vars']
    if status_str not in ("Optimal",):
        # E.g. timed out: return current best (CBC's incumbent is queryable).
        if model.objective is not None and pulp.value(model.objective) is not None:
            return int(round(pulp.value(model.objective))), elapsed, ctx['num_vars']
        return -2, elapsed, ctx['num_vars']

    obj = pulp.value(model.objective)
    if obj is None:
        return -1, elapsed, ctx['num_vars']
    return int(round(obj)), elapsed, ctx['num_vars']
