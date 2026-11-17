The Ozaki residue division uses shifts, so clang-cl no longer emits `__udivti3` and `__umodti3` for a 128-bit operator that lld-link does not pull into the DLL.
