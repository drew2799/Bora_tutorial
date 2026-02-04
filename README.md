This is a tutorial for the usage of Bora.jl. At the end of the tutorial, you will know how to:
- Load a set of trained Bora emulators
- Load data (measurements and covariances)
- Define a likelihoods and priors
- Run chains and compute MLE and MAP

# Preliminaries

The suggested way of installing julia is through [juliaup](https://github.com/JuliaLang/juliaup).
This can be done on UNix systems with the command (see also instructions at the above link!)

curl -fsSL https://install.julialang.org | sh

After you installed julia, there are two steps that are missing:

- Installing the jupyter kernel. This can be done by starting a julia session (simply type julia in the terminal), then in the Julia REPL type "using Pkg; Pkg.add("IJulia"); Pkg.build("IJulia")". After doing that, you succesfully have installed the Jupyter kernel.
- Now you have to install all the required Julia packages. This can be done by opening a Julia session *in this folder*, with the command "julia --project=.". This will create a new session, that will automatically read the Project and Manifest files, reading the info about packages. Then, you have to type the command "using Pkg; Pkg.instantiate()". Now, the installation process will start! This might take a few minutes.

After this has been done, you are ready! Just start a jupyter notebook session in this folder and execute the jupyter notebook present in this folder.


# Bora.jl: Getting Started

Welcome to the **Bora.jl** tutorial. By the end of this guide, you will be able to:

* Load pre-trained **Bora emulators**.
* Import datasets, including measurements and covariance matrices.
* Define cosmological **likelihoods and priors**.
* Execute MCMC chains and calculate **MLE** (Maximum Likelihood Estimation) and **MAP** (Maximum A Posteriori) values.

---

## Preliminaries

### 1. Install Julia

The recommended way to manage your Julia installation is via [juliaup](https://github.com/JuliaLang/juliaup). For Unix-based systems (Linux/macOS), run the following command in your terminal:

```bash
curl -fsSL https://install.julialang.org | sh

```

### 2. Configure the Jupyter Kernel

To use Julia within Jupyter notebooks, you must install the `IJulia` package. Open a Julia session by typing `julia` in your terminal, then execute:

```julia
using Pkg
Pkg.add("IJulia")
Pkg.build("IJulia")

```

### 3. Environment Setup & Dependency Installation

This tutorial uses a dedicated Julia environment to ensure version compatibility. Navigate to this project folder in your terminal and launch Julia using the local project flag:

```bash
julia --project=.

```

Inside the Julia REPL, instantiate the environment to download and precompile all required packages (this may take a few minutes):

```julia
using Pkg
Pkg.instantiate()

```

---

## Running the Tutorial

Once the installation is complete, you are ready to begin. Launch a Jupyter session from this directory:

```bash
jupyter lab

```

Open the `.ipynb` file provided in this folder to start the tutorial.
