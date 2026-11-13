-- Copyright (c) d-SEAMS developers
-- Read one frame and print the CHILL+ clathrate count.
local dseams = require("dseams")

dseams.set_num_threads(1)
local cloud = dseams.read(arg[1], {all = true})
local chill = dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})
print(string.format("lua clathrate %d", chill.counts.clathrate))
cloud:close()
