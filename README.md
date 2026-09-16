# EnKF-C (Ensemble Kalman Filter in C)

**EnKF-C** provides a compact generic framework for off-line data assimilation (DA) into large-scale layered geophysical models with the ensemble Kalman filter (EnKF).

Following are its other main features:
* Coded in C for GNU/Linux platform;
* Model-agnostic;
* Can conduct DA in EnKF, ensemble optimal interpolation (EnOI), or hybrid EnKF/EnOI modes;
* Permits multiple model grids;
* Can handle rectangular, curvilinear, or unstructured horizontal grids, z, sigma or hybrid vertical grids.

**EnKF-C** is coded for simplicity, scalability and robustness. To handle as large systems as possible it uses shared memory capabilities of MPI-3. Here is a snapshot of ensemble spread of sea surface temperature from the 96-member EnKF ocean forecasting system with MOM5 based OFAM3 model (51 x 1500 x 3600 grid), assimilating about 14M super-observations at each 3-day cycle.

![](sst-spread.png)

For more information see [README](https://github.com/sakov/enkf-c/blob/master/enkf/README) and [user guide](https://github.com/sakov/enkf-c/blob/master/enkf/doc/enkf-userguide.pdf). (An older version of the user guide is also available from [arXiv](http://arxiv.org/abs/1410.1233).) Have a feel for how the code works by running the included examples.

Checkout **EnKF-C** by running `git clone https://github.com/sakov/enkf-c`.

***

## License
EnKF-C is public software. See the `LICENSE` file for details.

## Dependencies
Ensure your environment has the following libraries and compilers installed before building:
* **Compilers:** `clang` or `gcc`
* **Parallel Computing:** `openmpi`
* **Data Formats:** `libnetcdf`
* **Linear Algebra:** `liblapack` (or for optimized performance, Intel's MKL implementation: `libmkl_rt`)

## Compiling
EnKF-C is developed natively for the GNU/Linux platform. To build the project:

1. Create a `make.inc` configuration file adapted to your architecture (refer to the templates available in the `arch/` directory).
2. Execute the compilation script from the command line:
   ```bash
   make
   ```
3. Upon successful compilation, the following binaries will be generated in the `bin/` directory:
   * `bin/enkf_prep`
   * `bin/enkf_calc`
   * `bin/enkf_update`
   * `bin/ens_diag`

## Code Indentation
To maintain the structural indentation style native to the EnKF-C code base:
1. Download the modified indent tool configuration: [indent-2.2.8a-mod.tar.gz](https://github.com "sakov/enkf-c Release Component")
2. Run the formatting rule from this directory (the style constraints are hardcoded into `.indent.pro`):
   ```bash
   make indent
   ```

## Starting Up
* **Run Example 1:** Follow the setup steps explicitly detailed inside `examples/1/README`.
* **User Manual:** For a thorough technical manual and operational guide, open the comprehensive PDF documentation located at `doc/enkf-userguide.pdf`.

## Notes
1. EnKF-C is known to perform reliably when compiled with GCC and linked with openmpi v4.1.0.
2. For best performance model dumps should be chunked by layers.

***

## Troubleshooting & Desktop Execution

### Suppressing Open MPI 5 Hardware Warnings
When executing local test examples or assimilation pipelines on a personal workstation or desktop computer using **Open MPI 5**, you may observe an aggressive stream of hardware-related warnings printed to the terminal console. These typically mention missing enterprise-scale cluster interconnects, fabric drivers, or fabric managers (such as UCX, InfiniBand, or RoCE components).

While these warnings are harmless on a local single-node architecture, they significantly clutter runtime diagnostics. To restrict parallel process communication to a quiet, local desktop sandbox environment, export the following environment configurations in your active shell profile before initiating standard binary runs:

```bash
# Force the robust legacy point-to-point management engine
export OMPI_MCA_pml=ob1

# Restrict network operations explicitly to self-loopback and basic TCP 
export OMPI_MCA_btl=self,tcp

# Prevent Open MPI from attempting to initiate UCX drivers for shared memory layers
export OMPI_MCA_osc="^ucx"
```

> **[WARNING] HPC Environment Notice:** Do **not** include these parameter constraints in batch job-submission scripts running EnKF-C on highly specialized distributed supercomputing systems or computing clusters. Forcing legacy engines and basic TCP protocols will bypass the high-speed network fabrics (like InfiniBand backbones) and introduce immense communication bottlenecks into large ensemble cycles.

***

## Reference
Sakov, P., 2014: *EnKF-C user guide*. Technical Report, Bureau of Meteorology. Available online from the [EnKF-C arXiv Paper Repository](https://arxiv.org).
