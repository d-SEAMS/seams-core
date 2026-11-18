Split docs reader text checks from live package-registry validation. Meson
`docs_reader_text` now stays offline while still parsing install claims, and CI
runs live registry claim checks in a dedicated docs workflow job.
