# mop-rs
Inspired by [mop](https://github.com/RILAB/mop), mop-rs is a simple tool to identify sites with sufficient depth for genotyping.

Currently, mop-rs only supports a single [d4](https://github.com/38/d4-format) file as input. Per sample d4 files can be created using [mosdepth](https://github.com/brentp/mosdepth), and then merged using `d4tools merge`.

## Installation
Binaries are available on the releases page. 

Alternatively, you can clone the repo and build from source using `cargo`:
```
cargo build --release -p moprs
```

## Usage
```
Calculates callable sites from depth statistics.

Usage: moprs [OPTIONS] <--d4 <D4_FILE>>

Options:
      --d4 <D4_FILE>
          The path to the D4 file
  -m, --min-depth <MIN_DEPTH>
          Minimum depth to consider site callable per individual [default: 0]
  -M, --max-depth <MAX_DEPTH>
          Maximum depth to consider site callable per individual [default: inf]
  -d, --depth-proportion <DEPTH_PROPORTION>
          Proporition of samples passing thresholds at site to consider callable [default: 1]
  -u, --min-mean-depth <MEAN_DEPTH_MIN>
          Minimum mean depth across all samples at site to consider callable [default: 0]
  -c, --output-counts
          Output number of individuals callable at site. EXPERIMENTAL v0.1.0
  -t, --threads <THREADS>
          Number of threads to use [default: 1]
  -h, --help
          Print help
  -V, --version
          Print version
```