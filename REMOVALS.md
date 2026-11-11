These deletions remove library guards.
They do not remove tests.
Each guard is a check inside src/cage_affiliation.cpp.
The caller chain is the reason the failing side has no input.
basalPair tested a ring against itself before shareAtoms.
A ring shares every atom with itself, so shareAtoms of that pair returns
false.
batchSweep draws the partner from basalCandidates.
stackingPlanes draws the partner from the same list.
That list does not contain the asking ring.
The self-pair test had no caller that could succeed.
The self-candidate continue compared a candidate with the asking ring.
An on-ring neighbour is skipped before ringsOf.
The asking ring is never offered back.
The continue had no input.
basalCandidates returned an empty list for an index outside the cache.
batchSweep calls it with an index from 0 through the ring count minus 1.
stackingPlanes calls it with an index from that same range.
prepare sizes the cache to the ring count.
An index outside the cache had no caller.
The prismatic loop tested that a ring index lay inside the ring list.
prismaticOf calls findPrismatic.
findPrismatic appends indices taken from the ring index.
Those indices already lie in the ring list.
An index outside that list had no caller.
batch took a context pointer and built a context when the pointer was
null.
update calls batch twice.
Both calls passed a null pointer.
The function now builds the context itself.
A non-null pointer had no caller.
The constructor allocated the implementation with new.
operator new throws when allocation fails.
The constructor now calls make_unique.
The null result is not a value the constructor stores.
The hot path tested batchHot.
It also tested warm.
It also tested primed.
The empty-frame path clears batchHot.
The same path clears warm.
The same path clears primed.
batchHot is set only after fastAffiliation.
fastAffiliation builds the warm state.
primed is still set when batchHot is set.
The two extra flags added no reachable case.
The seeded rewrite is a separate commit.
It names the strict flag and the permissive flag before the assignment.
A kept atom is marked when either flag is set.
That rewrite does not delete a range test.
Coverage of this file before these deletions is in the pre-removal log.
Coverage after these deletions is in the post-removal log.
We recorded one commit for each library guard.
