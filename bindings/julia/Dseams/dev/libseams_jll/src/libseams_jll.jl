# Stand-in for libseams_jll in a seams-core checkout, until the JLL is in
# the General registry. It has no binaries: Dseams then opens
# SEAMS_LANG_LIB, the libseams_lang preference, or libseams_lang on the
# loader path.
module libseams_jll
is_available() = false
end
