<!-- https://docs.dseams.info/api/enum_nucleation_CommittorFate.md -->

# enum nucleation::CommittorFate {#enum-nucleation-committorfate}

Outcome of one replica. committed means the final size is still strictly above that replica’s own grow-basin barrier. truncated means the series ends on the cap 2N while that cap is also its maximum: the run was cut there and is not a committor sample. The barrier is never 2N itself.
