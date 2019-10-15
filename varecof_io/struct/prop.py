""" propargyl test
"""

import automol

ICH = 'InChI=1S/C3H3/c1-3-2/h1H,2H2'
GEO = automol.inchi.geometry(ICH)
print(automol.geom.string(GEO))

GRAPH = automol.geom.graph(GEO)
SITES = automol.graph.resonance_dominant_radical_atom_keys(GRAPH)
print(GRAPH)
print(SITES)

GRAPH2 = ({0: ('C', 1, None), 1: ('C', 2, None), 2: ('C', 0, None)}, {frozenset({0, 2}): (1, None), frozenset({1, 2}): (1, None)})
SITES2 = automol.graph.resonance_dominant_radical_atom_keys(GRAPH2)
print(GRAPH2)
print(SITES2)


