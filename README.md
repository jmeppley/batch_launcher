# batch_launcher
Python script to fragment tasks for parallel execution

batch_launcher understands many basic bioinformatics file types and can spawn parallel jobs on a single machine or accross a cluster using SGE/gridengine or SLURM.

## Testing
### Set up conda env
`conda env create -p ./test.env -f test_data/conda.env.yaml`

### Activate environment
`conda activate ./test.env`

### BATS tests
`bats test_data/bats`