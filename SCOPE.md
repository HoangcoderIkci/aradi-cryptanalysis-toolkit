# Scope & limitations (read this first)

This repository must never be described as "breaking ARADI". The honest scope is:

**What this work DOES:**
- An independent C implementation of ARADI that matches the official NSA test vector
  (ePrint 2024/1240).
- An independent verification (via MILP, PuLP+CBC) of the **AABB cube-distinguisher**
  that holds after **5 rounds** of ARADI.
- A two-stage key-recovery analysis showing that, for the **full 16-round ARADI**, the
  attack stays at complexity **≈ 2¹⁴¹** — i.e. **practically infeasible** (full ARADI is
  **NOT broken**).
- A practical cube attack on a **MODIFIED, 6-round key schedule**, which recovers a
  **32-bit pair `(u, ℓ)`** in a few seconds, 100/100 success on random master keys.

**What this work does NOT do — never claim otherwise:**
- It does **NOT** break full 16-round ARADI.
- It does **NOT** present a new attack on the real, unmodified ARADI.
- It does **NOT** recover the 256-bit master key. The modified schedule projects the key
  through `K_base = K₀ ⊕ … ⊕ K₇` (32 bits), a **non-invertible** map — recovering the
  master key from `(u, ℓ)` is **information-theoretically impossible**, not merely hard.

**Vietnamese / Tiếng Việt:** Tấn công áp dụng cho **biến thể 6 vòng ĐÃ SỬA key schedule**;
ARADI đầy đủ vẫn an toàn (~2¹⁴¹) — **không phá ARADI**, **không khôi phục master key 256-bit**
(chỉ phục hồi cặp 32-bit `(u, ℓ)`). Sự rõ ràng về phạm vi này là chủ ý.
