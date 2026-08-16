Runtime knobs (cutoff, *k*, occupancy, residency) come from a
twelve-factor stack: code defaults, an optional dotenv file
(``SEAMS_CONFIG`` or ``./seams.env``), then the environment, then
CLI flags. ``seams --print-config`` dumps the resolved table.
The 2020 ``conf.yaml`` driver is not this layer.
