print("classifyTemplates / SOAP / Voronoi on the loaded frame")

for frame = targetFrame, finalFrame, frameGap do
  resCloud = readFrameOnlyOne(trajectory, frame, resCloud, oxygenAtomType, false, {0, 0, 0}, {0, 0, 0})
  nList = neighborList(cutoffRadius, resCloud, oxygenAtomType)
  local hits = classifyTemplates(resCloud, nList, 12)
  local counts = {}
  for i = 1, #hits do
    local name = hits[i].name
    counts[name] = (counts[name] or 0) + 1
  end
  for name, n in pairs(counts) do
    print(string.format("template %s %d", name, n))
  end
  local spec = soapSpectrum(resCloud, 0, nList, 3, 6, cutoffRadius)
  print(string.format("soap0_len %d", #spec))
  local q = steinhardtQlVoronoi(resCloud, cutoffRadius, 6)
  print(string.format("voronoi_q6_n %d", #q.ql))
end
