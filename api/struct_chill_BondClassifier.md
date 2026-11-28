<!-- https://docs.dseams.info/api/struct_chill_BondClassifier.md -->

# struct chill::BondClassifier {#struct-chill-bondclassifier}

## Overview {#overview}

One rule set for classifying bond correlations $c_{ij}$. [More…](#details-structchill-1-1-bond-classifier)

## Detailed Documentation {#detailed-documentation}

One rule set for classifying bond correlations $c_{ij}$.

A bond is staggered when $c_{ij} \le$ staggeredMax, eclipsed when it falls in \[eclipsedMin, eclipsedMax\], and out of range otherwise. The CHILL and CHILL+ water rules are two instances; other hydrogen-bonded or tetrahedral materials register their own thresholds and coordination.

### Fields {#fields}

Inclusive staggered edge, c_ij <= staggeredMax.

Inclusive lower edge of the eclipsed band.

Inclusive upper edge of the eclipsed band.
