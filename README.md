# batch_launcher
Python script to fragment tasks for parallel execution

batch_launcher understands many basic bioinformatics file types and can spawn parallel jobs on a single machine using local threads.

Local execution is the default. The `-X sge` and `-X slurm` options are no longer supported.

## Testing
### Set up conda env
`conda env create -p ./test.env -f test_data/conda.env.yaml`

### Activate environment
`conda activate ./test.env`

### BATS tests
`bats tests/bats`

### Python unit tests
`python -m pytest`