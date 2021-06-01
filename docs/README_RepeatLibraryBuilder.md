# Repeat Library Builder pipeline

## Quickstart (NBIS Staff)

```bash
module load Singularity
nextflow run -profile nbis,singularity /path/to/RepeatLibraryBuilder.nf
```

Or:
```bash
nextflow run -profile nbis,conda /path/to/RepeatLibraryBuilder.nf
```


## Usage

### Parameters

- General:

Parameters to the workflow can be provided either using `--parameter` notation or via a config file as follows:

`params.config`:
```
// Workflow parameters

// Nextflow parameters
resume = true
workDir = '/path/to/temporary/workspace'
conda.cacheDir = "$HOME/.nextflow/conda"
singularity.cacheDir = "$HOME/.nextflow/singularity"
```

Run nextflow with config file:
```bash
# Open screen terminal
screen -S my_nextflow_analysis
# Load Nextflow environment with conda
conda activate nextflow-env
# Load Singularity for Nextflow to use -profile singularity
module load Singularity
# Run Nextflow analysis
nextflow run -c params.config -profile nbis,singularity /path/to/RepeatLibraryBuilder.nf
```

### Workflow Stages

1.
