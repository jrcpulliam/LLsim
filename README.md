# Main README file for LLsim package

> Line list simulator for use in testing methods and practicing outbreak analysis

## Flow for analysis

- Input simulation parameters; Create population and run epidemic simulation; Output population and perfectly observed cases
- Input perfectly observed cases and observation parameters; Determine observed cases and observation delays; Output line list

## Development notes

### Kinds of Inputs and Outputs

- Inputs: parameter values (for default simulation function minimally includes incubation period mean, infectious period mean, R0, population size, cfr; optional: shape parameters)
- Inputs: control parameters (seed, maximum simulation time)
- Other inputs: reporting parameters (should this be part of parameter values above? probably not, as it will be a separate function and could be used on data generated outside the package)
- Inputs could be tabluar, key value, or just function arguments... need to decide
- Outputs: line list data (and optionally full population data) - format specified by user (options: CSV, Rdata / Rds, incidence object??, other/s?); all output is tabluar
- Outputs: vignettes, which include examples (probably including example analysis, potentially using RECON packages `incidence` and `epicontacts`)

### Intermediate Results

- The full line list (perfect observation) is intermediate data, in some sense, falling between the simulation step and the observation step, but for testing methods will be useful data for comparison so should be available to the user

### Testing and validation

- Should test that the required input files are available (if that is how input works)
- (What other tests are necessary? Certainly some...)
- Should probably check dates or allow users to input various dates

