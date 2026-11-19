<!-- https://docs.dseams.info/api/struct_chill_QlCutoffStages.md -->

# struct chill::QlCutoffStages {#struct-chill-qlcutoffstages}

Wall seconds of the last steinhardtQlCutoff call on this thread: the cutoff list, then q_l and any average from it. `set` is false when that call took a path that does not time the two apart (orders 8 and 12, or a cell the solid-harmonic kernel declines).
