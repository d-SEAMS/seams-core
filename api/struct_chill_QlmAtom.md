<!-- https://docs.dseams.info/api/struct_chill_QlmAtom.md -->

# struct chill::QlmAtom {#struct-chill-qlmatom}

## Overview {#overview}

This is the local orientational bond order parameter $q_{lm}$, of length $2l+1$. [More…](#details-structchill-1-1-qlm-atom)

## Detailed Documentation {#detailed-documentation}

This is the local orientational bond order parameter $q_{lm}$, of length $2l+1$.

This complex vector is averaged over the four nearest neighbours, according to the following equation:

$$
q_{lm}(i) = \frac{1}{N_b(i)} \Sigma_{j=1}^{N_b(i)} Y_{lm}(r_{ij})
$$

Here, $N_b(i)=4$ is the number of nearest neighbours for the molecule $i$. This struct contains specifically:

- A complex vector of length $2l+1$, calculated according to the equation above
