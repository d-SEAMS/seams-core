The nucleation notebook rule declares ``figshare-incremental.json``
as an output. A later ``test -s`` stamp cannot see a side-effect
file because Snakemake deletes declared outputs before the job.
