# Main Features / Design Overview

* Subworkflows

Subworkflows are used extensively to create reusable building blocks.

For example, a basic molecular dynamics workflow is composed of minimization, equilibration, and production steps, the equilibration step is composed of constant volume (NVT) and constant pressure (NPT) steps, and each of those are composed of primitive implementation-specific steps. If we then want to do a stability analysis, we should be able to incorporate the molecular dynamics workflow as a black box, and we should only have to append the stability analysis subworkflow. See [here](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs/tutorial.wic)!

* Multiple implementations

There are often several implementation engines that implement the same algorithm (e.g. amber / gromacs / namd). In principle, each implementation ought to be exactly interchangeable, but in practice implementations may randomly crash, etc. For this reason, we want the ability to arbitrarily switch implementations at any step. Moreover, different users may be familiar with different implementations, but we still want to compose together their workflows.

For example, we should be able to compose system setup using amber/tleap, equilibration using namd, and metadynamics using gromacs. File format conversions should be automatically inserted between steps. This is not possible with other 'implementation-independent' software packages which require a single implementation to be fixed at the beginning and used throughout.

See [static dispatch](advanced.md#static-dispatch)

* Automated Real-time Analysis & Plots

For quality control purposes, it is highly desirable to have fully automated analyses. This is particularly important for Virtual Screening, where the number of simulations is so large that a user cannot possibly manually inspect each simulation and intervene. For example, a stability analysis ought to be done before an expensive binding free energy calculation is performed.

We support iteratively running an arbitrary analysis workflow in real-time (i.e. while the simulation is still running) and plotting the results. Timeseries data can be automatically segmented and clustered into stastically separate probability distributions, which are beautifully reflected in the associated histograms.