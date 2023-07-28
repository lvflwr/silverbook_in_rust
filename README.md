# silverbook_in_rust
A sample implementation in Rust of "Numerical Methods for Computational Fluid Dynamics" by Kozo Fujii (藤井孝蔵著『流体力学の数値計算法』), commonly known as the Silver Book (銀本).


## Contents
This repository contains the following:
- A set of sample code of the book (under `./section_*/package_name/*`),
- Input files to run the code (under `./inputs/section_*/package_name/*`),
- Scripts to visualize the results of the code (under `./plots/section_*/package_name/*`).

Each chunk of code is implemented as a single package.

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
cargo run -p package_name --bin binary_crate_name

# e.g.
cargo run -p bad_upwind --bin exec_good_upwind
```

The output files are generated under `./outputs/section_*/package_name/*`.

You can change the input parameters by editing the input files under `./inputs/section_*/package_name/*`.


## Visualization
We provide brief scripts to visualize the results.

### Install gnuplot
In order to visualize the results, you need to install gnuplot.

See, for example, [the official manuals](http://www.gnuplot.info/download.html).

### Run the gnuplot scripts
Run the following command to generate the figures.
```shell
gnuplot plots/section_*/package_name/script_name.gp

# e.g.
gnuplot plots/section_1/bad_upwind/exec_good_upwind.gp
```

The figures are generated under `./outputs/section_*/package_name/*`.


## Documentation
If you want an overview of the sample code, including formulations, schemes, etc., you can refer to the documentation.

The documentation also includes the input and output formats for executing the sample code.

### Generate and see the documentation
Run the following command.
```shell
cargo doc --no-deps --open
```

You can see the documentation in your browser.


## License
The MIT license.

See [LICENSE](LICENSE) for details.
