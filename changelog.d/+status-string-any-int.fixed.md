`seams_status_string` reads its argument as an `int`, so a code outside `seams_status` returns "unknown status" without undefined behaviour in the C++ library; clang's UBSan enum check aborted on it.
