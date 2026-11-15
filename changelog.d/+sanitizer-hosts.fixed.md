A sanitized build does not preload AddressSanitizer into Julia, R, or `nm`. Those hosts abort in their own startup. The runtime is mapped when they load `libseams`.
