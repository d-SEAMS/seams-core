A ThreadSanitizer build gives plumed, gmx and lmp the unsanitized plugin, because those hosts cannot map `libclang_rt.tsan`. Interpreter shims map that runtime before they dlopen `libseams`.
