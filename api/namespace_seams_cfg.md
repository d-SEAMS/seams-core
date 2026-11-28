<!-- https://docs.dseams.info/api/namespace_seams_cfg.md -->

# namespace seams::cfg {#namespace-seams-cfg}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Functions {#global-functions}

Scan argv for config PATH / config=PATH.

Empty if absent.

Apply a dotenv file: set KEY=VAL only when KEY is not already in the environment.

Throws if required and the file is missing.

Load file (if any) then read the environment into a [Runtime](struct_seams_cfg_Runtime.md#doxid-structseams-1-1cfg-1-1-runtime).

Write [Runtime](struct_seams_cfg_Runtime.md#doxid-structseams-1-1cfg-1-1-runtime) fields into the environment so libraries that only call getenv (linkcell, OpenMP offload) see the CLI result.
