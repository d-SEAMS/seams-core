In OpenMP builds the threaded neighbour rows gave an atom missing from `idIndexMap` its neighbours' IDs as a row header, and gave its partners -1. That atom now keeps an empty row.
