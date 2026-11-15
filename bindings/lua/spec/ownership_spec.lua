-- Copyright (c) d-SEAMS developers
-- busted: errors are raised, and a closed frame stays closed.
local dseams = require("dseams")

describe("dseams ownership", function()
  it("rejects a cell with no volume", function()
    local sys = dseams.core.system()
    local ok, err = pcall(function()
      sys:set_ortho({0, 0, 0}, 0, 0, 0, 0, 1, 1)
    end)
    assert.is_false(ok)
    assert.is_truthy(tostring(err):find("box edge lengths must be positive", 1, true))
    sys:close()
  end)

  it("labels the eight-atom cell cubic", function()
    local a = 2.75 * 4.0 / math.sqrt(3.0)
    local fcc = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}
    local basis = {{0, 0, 0}, {0.25, 0.25, 0.25}}
    local xyz = {}
    for _, site in ipairs(fcc) do
      for _, shift in ipairs(basis) do
        xyz[#xyz + 1] = (site[1] + shift[1]) * a
        xyz[#xyz + 1] = (site[2] + shift[2]) * a
        xyz[#xyz + 1] = (site[3] + shift[3]) * a
      end
    end
    local cloud = dseams.system(xyz, {0, 0, 0, a, a, a, 0, 0, 0})
    local labels = dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})
    assert.are.equal(8, labels.counts.cubic)
    cloud:close()
  end)
end)
