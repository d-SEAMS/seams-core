`-Dcapi_sanitize=address,undefined` builds the C API suite under ASan and UBSan at once, and every UBSan build also checks float-to-int overflow, which GCC leaves out of `-fsanitize=undefined`.
