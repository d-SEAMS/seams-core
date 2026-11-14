The backtrace signal handler now outlives static destruction. A sanitizer abort at exit hit the freed handler stack and hung `lammps_hook` until its timeout.
