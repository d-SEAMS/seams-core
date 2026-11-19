<!-- https://docs.dseams.info/api/global.md -->

# Global Namespace {#global-namespace}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Macros {#macros}

ABI generation of this header.

Bumped only on a breaking change.

Marks a symbol exported by libseams.

The library is compiled with hidden visibility and SEAMS_BUILDING defined, so only declarations carrying this macro are exported. Windows consumers get dllimport. Define SEAMS_STATIC when linking a static archive.

Marks an entry point kept for ABI 1 callers.

`msg` names the replacement.
