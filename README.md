# fastqc-rs (reduced version)
This fork is prepared for benchmarking the base sequence UDF from [polars-bio](https://github.com/MikolajSzawerda/polars-bio/).

#### Source

Download the source code and within the root directory of source run

    cargo install

## Usage

```
cargo run -- -q path/to/my_sequence.fastq > report.html
```

Arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -q --fastq 	       |	-           |The path to the FASTQ file to use
| -k          | 5           |The length k of k-mers for k-mer counting
| -s --summary          | -           |Creates an output file for usage with [MultiQC](https://multiqc.info) under the given path
