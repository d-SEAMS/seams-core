The in-tree Python package now has the documented Frame handle.
`pydseams.read`, `from_ase`, `to_ase`, `to_solvis`, `chill_plus`, `cages`,
`cages_by_signature`, `steinhardt`, `rdf_2d`, `monolayer_rings`,
`find_prisms`, and `pydseams.band.ring_census` call `libseams_lang`.
The tutorials keep those calls. `python_frame_parity` compares them with
pydseamslib 2.11.1 on the same dumps.
