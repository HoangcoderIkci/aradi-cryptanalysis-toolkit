# Neural distinguisher vs the AABB distinguisher  (roadmap — module ⑤)

**Question.** A Gohr-style neural distinguisher (CRYPTO 2019) learns statistical structure
that analytic methods sometimes miss. On the **same** ARADI testbed, how does a trained
neural distinguisher compare to the **analytic AABB cube-distinguisher** we already verified?

**Plan.**
1. Use the C impl ([`../impl/`](../impl)) / reference impl
   ([`../python/aradi_ref.py`](../python/aradi_ref.py)) as a fast data generator for
   round-reduced ARADI (real-vs-random pairs).
2. Train a small ResNet/MLP distinguisher per round count; measure accuracy vs #rounds and
   vs data complexity.
3. Put it **head-to-head** with AABB on the same rounds: where does each win? Does ML see
   anything the analytic distinguisher does not (and vice-versa)?
4. Write up the comparison — *formal + ML on one cipher*.

**Honest expectation.** ARADI's 128-bit block is much larger than Speck32; the neural
distinguisher may **not** beat AABB's round reach. That negative/comparative result is
itself worth reporting. Never describe ML here as "breaking" ARADI.

**Status:** not started. Owner: —.  Likely deps: `torch` (added later, kept out of core
`requirements.txt`).
