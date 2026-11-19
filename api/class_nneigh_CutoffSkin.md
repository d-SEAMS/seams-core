<!-- https://docs.dseams.info/api/class_nneigh_CutoffSkin.md -->

# class nneigh::CutoffSkin {#class-nneigh-cutoffskin}

## Overview {#overview}

In-situ cutoff list. [More…](#details-classnneigh-1-1-cutoff-skin)

## Detailed Documentation {#detailed-documentation}

In-situ cutoff list.

Candidates are built at `cutoff + skin` and reused while every atom and the box have moved by at most `skin/2`. The filtered list matches a fresh build at `cutoff`.
