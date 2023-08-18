# silverbook_in_rust
[![CI](https://github.com/lvflwr/silverbook_in_rust/actions/workflows/ci.yml/badge.svg)](https://github.com/lvflwr/silverbook_in_rust/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/lvflwr/silverbook_in_rust/graph/badge.svg?token=JRFZFVPQST)](https://codecov.io/gh/lvflwr/silverbook_in_rust)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A sample implementation in Rust of "Numerical Methods for Computational Fluid Dynamics" by Kozo Fujii (藤井孝蔵著『流体力学の数値計算法』), commonly known as the Silver Book (銀本).


## Contents
This repository contains the following:
- A set of sample code of the book (under `./section_*/package_name/*`),
- Input files to run the code (under `./inputs/section_*/package_name/*`),
- Scripts to visualize the results of the code (under `./plots/section_*/package_name/*`).

The code is organized into one package per topic.

`section_*` (e.g., `section_1`) in the path corresponds to the chapter number of the book.


## Usage
You can easily run the sample code in the following steps.

### Install Rust
Install Rust.

See, for example, [the official manuals](https://www.rust-lang.org/tools/install).

### Clone this repogitory
Run the following commands.
```shell
cd /path/to/your/work/directory
git clone git@github.com:lvflwr/silverbook_in_rust.git
cd silverbook_in_rust
```

### Run the sample code
Run the following command.
```shell
cargo run --example example_name

# e.g.
cargo run --example solve_transport_eq_by_good_upwind_method
```

The `example_name` corresponds to a file name under `./section_*/package_name/examples/*`.

The output files are generated under `./outputs/section_*/package_name/example_name/*`.

You can change the input parameters by editing the input files under `./inputs/section_*/package_name/example_name/*`.


## Visualization
You can use some scripts to visualize the results.

### Install gnuplot
In order to visualize the results, you need to install gnuplot.

See, for example, [the official manuals](http://www.gnuplot.info/download.html).

### Run the gnuplot scripts
Run the following command to generate the figures.
```shell
gnuplot plots/section_*/package_name/example_name/script_name.gp

# e.g.
gnuplot plots/section_1/bad_upwind/solve_transport_eq_by_good_upwind_method/plot_solution.gp
```

The figures are generated under `./outputs/section_*/package_name/example_name/*`.


## Documentation
If you want an overview of the sample code, including formulations, schemes, etc., you can refer to the documentation.

The documentation also includes the input and output formats for executing the sample code.

### Generate and see the documentation
Run the following commands.
```shell
cargo doc --examples --no-deps --open
```

You can see the documentation in your browser.
