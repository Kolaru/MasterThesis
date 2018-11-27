Master Thesis of Benoît Richard, done in 2018 in the University of Fribourg under the supervision of Shi Guiyuan and Prof. Yi-Cheng Zhang.

It is made public here to allow easy access to the code and to share it as an example of the use the Julia packages [`IntervalArithmetic.jl`]() and [`IntervalRootFinding.jl`]().

All code is available under MIT license, while the report, the presentation and all plots are published under the [Creative Commons Attribution-ShareAlike 4.0 License](https://creativecommons.org/licenses/by-sa/4.0/), for license compatibility with the [Konect database](http://konect.uni-koblenz.de/) which was use for the plots of the example networks.

Note that despite being made available here, most of the code will probably not work out of the box, even if you have all of the (numerous) dependencies installed. Possible reasons of failure include:
  - My `startup.jl` file add the `src/` folder to the `LOADPATH`.
  - The source data for real networks are not included in this repo, but they can be retrieved in the [Konect database](http://konect.uni-koblenz.de/), though.
  - Some functions requires to run other functions first to generate the needed data (typically all most files in the `Plot generation` folder do).
  - The working directory must be the root directory of this repo.
  - The needed version of `IntervalRootFinding` include a not yet merged [pull request](https://github.com/JuliaIntervals/IntervalRootFinding.jl/pull/97)
  - The `MastersDoctoralThesis.cls` used as document class to the `reportmain.tex` LaTeX document is not distributed to avoid any licensing problem as I may have modified it at some point (who knows ?). The original class file can be downloaded from [latextemplates.com](https://www.latextemplates.com/template/masters-doctoral-thesis).
  - This is a gigantic mess and I may have inadvertently left behind some broken files.

Dependencies includes (some my be missing from the list):
  - Julia
    - `IntervalArithmetic.jl`
    - `IntervalRootFinding.jl`
    - `Distributions.jl`
    - `PyCall.jl`
    - `ProgressMeter.jl`
    - `Espresso.jl`
    - `NLsolve.jl`
    - `StaticArrays.jl`
    - `ForwardDiff.jl`
  - Python
    - `numpy`
    - `matplotlib`
    - `mpmath`
    - `networkx`
