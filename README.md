# Capstone TOV Solver

This project solves the Tolman-Oppenheimer-Volkoff equations for tabulated neutron-star equations of state (EOS). The EOS tables used by the solver are generated from CompOSE data and stored in `data/`.

## Repository Layout

- `main.py`: entry point for generating a mass-radius sequence and plot.
- `src/`: TOV solver, EOS interpolation, constants, sequence generation, and plotting code.
- `data/`: EOS tables used by the Python solver and NICER posterior CSV files used for plotting.
- `outputs/`: generated star-sequence text files and the combined mass-radius PDF.
- `eos_dat_compilers/`: one CompOSE working directory per EOS. Each folder contains the downloaded CompOSE files, Fortran code, generated `eos.table`, and a README with source citations.
- `legacy/`: older single-file solver scripts kept for reference.

## Python Requirements

The Python solver requires:

- Python 3
- `numpy`
- `scipy`
- `matplotlib`

The code also uses Python standard-library modules including `os` and `pathlib`.

Install the dependencies with:

```bash
python3 -m pip install numpy scipy matplotlib
```

## Run the Solver

The active EOS is selected in `main.py`:

```python
EOS_file = 'eos_ds_cmf7.table'
```

The available project EOS tables are:

- `data/eos_abht_qmc_rmf1_unified_crust.table`
- `data/eos_apr_apr.table`
- `data/eos_bl_chiral.table`
- `data/eos_ds_cmf1.table`
- `data/eos_ds_cmf4.table`
- `data/eos_ds_cmf7.table`

To run the solver:

```bash
python3 main.py
```

This writes a sequence file to `outputs/star_sequence_<eos-name>.txt` and refreshes `outputs/all_eos.pdf`.

## EOS Table Format Expected by the Solver

`src/eos.py` reads the EOS table with:

```python
data = np.loadtxt(file_path)
P_MeVfm3 = data[:, 3]
eps_MeVfm3 = data[:, 4]
```

So each EOS table must be a whitespace-delimited ASCII table with at least five columns:

1. temperature `T`
2. baryon number density `n_b`
3. hadronic charge fraction `Y_q`
4. pressure `p` in `MeV fm^-3`
5. energy density `epsilon` in `MeV fm^-3`

The current CompOSE exports in `eos_dat_compilers/*/eos.table` use this format.

## Regenerate EOS Tables From CompOSE

Each subdirectory under `eos_dat_compilers/` is a separate CompOSE working directory. The current folders are:

- `eos_dat_compilers/eos_abht(qmc-rmf1)_unified_crust_compiler`
- `eos_dat_compilers/eos_apr(apr)_unified_crust_compiler`
- `eos_dat_compilers/eos_bl(chiral)_with_crust_compiler`
- `eos_dat_compilers/eos_ds(cmf)-1_with_crust_compiler`
- `eos_dat_compilers/eos_ds(cmf)-4_with_crust_compiler`
- `eos_dat_compilers/eos_ds(cmf)-7_with_crust_compiler`

These folders already contain the downloaded CompOSE data files and the CompOSE Fortran sources. To regenerate one table, work inside the corresponding folder.

### 1. Install Compiler Tools

You need GNU Make and a Fortran compiler.

For a simple ASCII export, HDF5 is not required. Build with `USE_HDF5=0` so the Makefile uses `gfortran` instead of `h5fc`:

```bash
cd "eos_dat_compilers/eos_ds(cmf)-7_with_crust_compiler"
make clean
make USE_HDF5=0
```

If you want HDF5 output, install HDF5 with Fortran support and build without `USE_HDF5=0`.

### 2. Select Output Quantities

Run the CompOSE program:

```bash
./compose
```

Choose Task 1, "Selection of Output Quantities".

For this project, select only the quantities needed by the Python solver:

- number of regular thermodynamic quantities: `2`
- regular quantity indices: `1 24`
- additional thermodynamic quantities: `0`
- composition quantities: `0`
- microscopic quantities: `0`
- error quantities: `0`
- output format: `1` for ASCII

This creates or updates `eos.quantities`. The important line should contain:

```text
1 24
```

where `1` is pressure and `24` is energy density.

### 3. Define the Tabulation Grid

Run `./compose` again and choose Task 2, "Definition of Tabulation Scheme and Parameter Values".

For the cold beta-equilibrium EOS tables used here, the only varying grid parameter is baryon number density `n_b`. Use:

- interpolation order for `n_b`: `1`
- tabulation scheme: `1` for loop form
- minimum `n_b`: use the minimum listed by CompOSE for that EOS
- maximum `n_b`: use the density range needed for the stellar sequence
- number of grid points: the current tables use `1700`
- grid scaling: `1` for logarithmic

The current generated tables use:

- ABHT(QMC-RMF1) unified crust: `n_b = 1.0e-11` to `1.28`, `1876` points, logarithmic spacing
- APR(APR): `n_b = 1.0e-11` to `1.34`, `1833` points, logarithmic spacing
- BL(Chiral) with crust: `n_b = 1.0e-11` to `1.29`, `283` points, logarithmic spacing
- DS(CMF)-1 with crust: `n_b = 1.0e-7` to `3.03e+00`, `1191` points, logarithmic spacing
- DS(CMF)-4 with crust: `n_b = 1.0e-7` to `3.03e+00`, `1129` points, logarithmic spacing
- DS(CMF)-7 with crust: `n_b = 1.0e-7` to `1.87e+00`, `1021` points, logarithmic spacing

This creates or updates `eos.parameters`.

### 4. Generate `eos.table`

Run `./compose` a third time and choose Task 3, "Generation of EoS Table".

The expected outputs are:

- `eos.table`: ASCII table consumed by this project
- `eos.info.json`: column metadata
- `eos.report`: CompOSE table report

Check that `eos.info.json` identifies column 4 as pressure and column 5 as energy density.

### 5. Copy the Generated Table Into `data/`

Rename the generated table according to the project convention:

```bash
cp eos.table ../../data/eos_ds_cmf7.table
```

Use the matching output name:

- ABHT(QMC-RMF1) unified crust: `data/eos_abht_qmc_rmf1_unified_crust.table`
- APR(APR): `data/eos_apr_apr.table`
- BL(Chiral) with crust: `data/eos_bl_chiral.table`
- DS(CMF)-1 with crust: `data/eos_ds_cmf1.table`
- DS(CMF)-4 with crust: `data/eos_ds_cmf4.table`
- DS(CMF)-7 with crust: `data/eos_ds_cmf7.table`

Then rerun:

```bash
python3 main.py
```

## Downloading a New EOS From CompOSE

To add a new EOS:

1. Download the EOS data files from CompOSE.
2. Create a new folder under `eos_dat_compilers/`.
3. Add all required CompOSE data files, including files such as `eos.t`, `eos.nb`, `eos.yq`, `eos.thermo`, `eos.init`, `eos.pdf`, and any composition or microphysics files included with the EOS.
4. Add the CompOSE Fortran sources and Makefile from `https://gitlab.obspm.fr/data_and_software_compose/code-compose`.
5. Build and run Tasks 1, 2, and 3 as described above.
6. Copy the resulting `eos.table` into `data/` with a clear lowercase filename.
7. Add or update a README in the EOS compiler folder with the CompOSE URL, access date, EOS name, and original references.

## Citation Notes

Each EOS compiler folder has its own README with the CompOSE source and citations. Keep those files with the data so the provenance remains clear.

For reports or papers, also cite CompOSE and the original EOS references in the bibliography, not only in the per-folder README.
