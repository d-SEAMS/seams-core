"""Second worked example of the manuscript: one deleted bond in a prism tube.

Builds the bond graph of the periodic hexagonal ice nanotube of
worked_example.py (six sites per layer, twelve layers), lists every primitive
ring of at most six members by brute force (all simple cycles, then the
Franzblau shortest-path test), deletes the bond between site 0 of layer 0 and
site 0 of layer 1, and lists the rings again. Writes one tab-separated table
of the counts the text quotes. Standard library only.

Usage: python3 worked_update.py OUT.tsv
"""

import collections
import sys

SIDES = 6
LAYERS = 12
# A cycle of a simple graph has at least MIN_SIZE members.
MIN_SIZE = 3
MAX_SIZE = 6


def tube(sides, layers):
    edges = set()
    for layer in range(layers):
        for vertex in range(sides):
            here = layer * sides + vertex
            around = layer * sides + (vertex + 1) % sides
            up = ((layer + 1) % layers) * sides + vertex
            edges.add(frozenset((here, around)))
            edges.add(frozenset((here, up)))
    return sides * layers, edges


def adjacency(n, edges):
    rows = collections.defaultdict(set)
    for edge in edges:
        a, b = tuple(edge)
        rows[a].add(b)
        rows[b].add(a)
    return [sorted(rows[v]) for v in range(n)]


def distances(rows, source, depth):
    seen = {source: 0}
    frontier = [source]
    for level in range(1, depth + 1):
        nxt = []
        for v in frontier:
            for w in rows[v]:
                if w not in seen:
                    seen[w] = level
                    nxt.append(w)
        frontier = nxt
    return seen


def cycles(rows, max_size):
    found = set()

    def grow(path, members):
        head = path[-1]
        for w in rows[head]:
            if w == path[0] and len(path) >= MIN_SIZE:
                found.add(canonical(path))
            elif w > path[0] and w not in members and len(path) < max_size:
                path.append(w)
                members.add(w)
                grow(path, members)
                members.discard(w)
                path.pop()

    for start in range(len(rows)):
        grow([start], {start})
    return found


def canonical(path):
    k = path.index(min(path))
    forward = path[k:] + path[:k]
    backward = [forward[0]] + forward[1:][::-1]
    return tuple(min(forward, backward))


def primitive(rows, ring):
    size = len(ring)
    for i in range(size):
        row = distances(rows, ring[i], size // 2)
        for j in range(i + 1, size):
            along = min(j - i, size - (j - i))
            if row.get(ring[j], size) < along:
                return False
    return True


def census(rows):
    return sorted(r for r in cycles(rows, MAX_SIZE) if primitive(rows, r))


def by_size(rings):
    count = collections.Counter(len(r) for r in rings)
    return ",".join(str(count.get(s, 0)) for s in range(MIN_SIZE, MAX_SIZE + 1))


def main():
    n, edges = tube(SIDES, LAYERS)
    before = census(adjacency(n, edges))
    cut = frozenset((0, SIDES))
    after = census(adjacency(n, edges - {cut}))
    changed = set(cut)
    touched = [r for r in before if changed & set(r)]
    lost = sorted(set(before) - set(after))
    gained = sorted(set(after) - set(before))
    carried = [r for r in before if not changed & set(r)]
    assert all(r in after for r in carried)
    assert all(changed & set(r) for r in lost + gained)
    rows = [
        ("sites", n),
        ("rings_before_3_to_6", by_size(before)),
        ("rings_after_3_to_6", by_size(after)),
        ("rings_before", len(before)),
        ("rings_after", len(after)),
        ("changed_atoms", len(changed)),
        ("rings_through_changed", len(touched)),
        ("rings_carried", len(carried)),
        ("rings_lost", len(lost)),
        ("rings_gained", len(gained)),
        ("lost_sizes", ",".join(str(len(r)) for r in lost)),
        ("gained_sizes", ",".join(str(len(r)) for r in gained)),
    ]
    text = "quantity\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in rows)
    with open(sys.argv[1], "w") as handle:
        handle.write(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
