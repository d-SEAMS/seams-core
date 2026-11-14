The Python and Lua bindings read neighbour lists through `seams_neighbours_copy` instead of the deprecated `seams_neighbours_csr` borrow, which stays only for ABI 1 callers.
