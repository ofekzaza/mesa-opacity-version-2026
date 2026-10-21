# Intro

Research code base for the porosity-reduction extension to MESA, by Ofek Eshet, under the guidance of Prof. Nir Shaviv.

In order to simulate massive and supermassive stars with MESA, we add an opacity (porosity) reduction algorithm, which is a 1D approximation of the real 3D effect inside massive stars.

With this code base we managed to simulate stars up to 100 solar masses, from creation up to Si burning (the 10 solar mass models stop at C burning instead, since they never ignite Si).

## Model logic

- **Target**: reduce the effective Eddington ratio (gamma = L / L_Edd) to about 0.8.
- **Trigger**: whenever gamma exceeds 0.8, the opacity is reduced by a few percent.
- **Rate**: the reduction is applied slowly, in small steps, until gamma reaches the target.

The reduction percentage is computed from the current gamma, the previously applied reduction factor, and a smoothing factor (`alpha = 0.1`). This makes the descent slower and the derivatives more continuous, which ensures a smooth MESA run.

In practice, the model is implemented as an opacity multiplier (`extra_opacity_factor`) that is written into the stellar structure and read back by MESA when the opacity is evaluated.

# Running it yourself

## Requirements

- MESA 26.4.1
- MESA SDK 26.3.2
- [uv](https://docs.astral.sh/uv/getting-started/installation/)
- Run `uv sync` to create the Python environment used by the plotting scripts.

## Folder naming convention

Folder names follow the format `{mass}m-{model}_{options}`.

Options:

- `_reduction` — the run includes mass loss (Dutch winds, `Dutch_scaling_factor = 0.66`).
- `_post` — the star is evolved starting from the main sequence instead of from core formation (`create_pre_main_sequence_model = .false.`).

Models:

- `ours` — the porosity-reduction model described above.
- `normal` — the default MESA model, with no opacity reduction.
- `mlt++` — the default MESA model with MESA's MLT++ (predictive convection) treatment enabled, i.e. `okay_to_reduce_gradT_excess = .true.` with `gradT_excess_max_change = 0.01`.
- `supereduction_a=2` / `supereduction_a=5` — MESA's built-in superadiabatic reduction, with a scale factor of 2 or 5.

Masses range from 10 to 100 solar masses, with both pre-main-sequence and post-main-sequence runs for every model. The 300 solar mass star is run with the `ours` model only.

## How to run

Run `sh run.sh`. This takes about a week, plus another week for the 300 solar mass `ours` reduction run.

To run a single model, use:

```sh
sh run_one.sh <folder_name>    # e.g. sh run_one.sh 10m-ours_reduction
```

Run logs are written to `run_log.txt`; read it to check the state of a run. Each folder also keeps its own `watchdog.log`.

### watchdog.sh

If a run fails for a reason that is not a MESA failure, the watchdog restarts it automatically using the `re` command, resuming from the last saved photo.

# Creating graphs

Once all simulations are finished, run `uv run main.py` to generate all the graphs in the `plotted_results` folder.

# Adding the porosity-reduction model to your own project

## Mixing length

First, we recommend using `MLT_option = 'TDC'` in `inlist_project`. This increases the stability of MESA significantly.

## Required changes

The opacity reduction is implemented as a MESA function extension, enabled with `use_other_opacity_factor = .true.` in the inlist file.

Then copy the `other_opacity_factor` subroutine and the `extra_opacity_factor_memory_target` variable declaration into your own `run_star_extras.f90`. If you do not have this file, copy ours into your project.

### Data

If you also want the extra columns that we added, you need to declare them in the history and profile output functions — see the `run_star_extras.f90` file for an example.
