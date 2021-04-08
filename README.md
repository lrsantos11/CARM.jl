# CARM - Circumcentered-Approximate-Reflection Method

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CARM

It is authored by Guilherme Araújo, Roger Behling, Yunier Bello-Cruz, Luiz-Rafael Santos.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
2. If one wants to run all the tests, do on Julia console
   ```
   julia> include(scriptsdir("runtests.jl"))
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.
