# LEMPo - RMHD

*Linear Eigenvalue Modular solver for Pressure- and current-driven modes (Resistive MHD) in tokamaks.*

## Overview

__LEMPo-RMHD__ is a modular eigensolver written in _MATLAB_ for studying the linear stability of internal magnetohydrodynamics (MHD) modes, in a tokamak plasma.
It computes linear growth rates and perturbations associated with moderate- to long-wavelength internal MHD instabilities, including **tearing modes**, **resistive interchange** modes, and **infernal modes** (in a torus).
The solver is based on a unified set of global equations for long-wavelength, pressure- or current-driven MHD modes, incorporating resistivity and compressibility in an inverse aspect ratio expansion. See the [References](#references) section for details.

<br>

<p align="center">
  <img src="docs/images/LEMPo_cover.png" width="300"><br>
  <em>
  The solver is named after Lempo: both a Kuiper Belt object (discovered in 1999)  
  and a demon from Finnish mythology. This reflects its focus on <strong>infernal</strong> modes,  
  and its ultimate integration within the VENUS-MHD framework.  
  (See credits in <a href="#image-sources">Image Sources</a>.)
  </em>
</p>

<br>

---

## Features

- **Solvers included:**  
  - `LEMPo_ideal.m` — Full ideal MHD solver with infernal corrections and explicit sidebands.
  - `IM_LEMPo_ideal.m` — Interchange model: simplified model without infernal corrections and explicit sidebands. This model is enough to retrieve ideal interchange modes.

- **Example scripts:**  
  - `m1n1_reversed_q_example.m` — Runs LEMPo_ideal for a reversed q-profile with m=1, n=1.
  - `m8n10_piecewise_q_example.m` — Runs LEMPo_ideal for a piecewise q-profile with m=8, n=10.
  - `m9n10_IM_scan_alpha_example` - Runs IM_LEMPo_ideal for the default q profile with m=9, n=10. Scan on the pressure drive. 

- **Utilities:**  
  - All required helper functions are in `src/utils/`.

---

## Folder Structure

```
.
├── src/
│   ├── LEMPo_ideal.m
│   ├── IM_LEMPo_ideal.m
│   └── utils/
│       └── [helper MATLAB functions]
├── examples/
│   └── ideal_examples/
│       ├── m1n1_reversed_q_example.m
│       ├── m8n10_piecewise_q_example.m
│       └── m9n10_IM_scan_alpha_example
└── setup.m
```

---

## Requirements

- **MATLAB** (tested on MATLAB R2024b)
- **Signal Processing Toolbox** 
- Can be run from a laptop

---

## Getting Started

### 1. **MATLAB Setup**

Before running any scripts or solvers, set up the project paths:

1. Open MATLAB and navigate to the project root folder.
2. Run the setup script:
   ```matlab
   setup
   ```
   This adds all source code and utilities to your MATLAB path.

### 2. **Running Example Scripts**

- All example scripts are in the `examples/` folder.
- Each example demonstrates how to set up profiles and options, and how to run the solver for a specific case.

### 3. **Using the Solvers Directly**

- The main solvers (`LEMPo_ideal.m` and `IM_LEMPo_ideal.m`) are in `src/`.
- These expect appropriate input profiles and options; see the example scripts for templates.

### 4. **Utilities**

- All dependencies and helper functions required by the solvers are in `src/utils/`.

### 5. **Notes**

- ⚠️ Note: When running the solver, you will see a MATLAB warning for the eigs solver: `Warning: The first input matrix, shifted by sigma, is close to singular or badly scaled ...` This is expected and does not indicate an error.

---

## Documentation

- **Inputs and outputs for each solver** are explained in the header comments of the `.m` files.
- **Example scripts** are fully commented and show typical usage.

---

## References

- Coste-Sarguet, M. & Graves, J. P. (2024). *Fundamental properties of ideal and resistive infernal modes in tokamaks*. Plasma Physics and Controlled Fusion. [https://doi.org/10.1088/1361-6587/ad5ff2](https://doi.org/10.1088/1361-6587/ad5ff2) 
- Graves, J. P., Coste-Sarguet, M. & Wahlberg, C. (2022). *Pressure driven long wavelength MHD instabilities in an axisymmetric toroidal resistive plasma*. Plasma Physics and Controlled Fusion, 64(1), 014001. [https://doi.org/10.1088/1361-6587/ac3496](https://doi.org/10.1088/1361-6587/ac3496)

---

## License

**If you use this code, or a modified version of it, in published work, please cite:**

> Coste-Sarguet, Margot and Graves, Jonathan P (2024).  
> *Fundamental properties of ideal and resistive infernal modes in tokamaks*.  
> Plasma Physics and Controlled Fusion.  
> [https://doi.org/10.1088/1361-6587/ad5ff2](https://doi.org/10.1088/1361-6587/ad5ff2)

---

## Image Sources

The symbolic image of LEMPo-RMHD combines the following materials:
- Planet Venus: NASA image (public domain)
- Kuiper Belt object: from [[Solar System Wiki](https://thesolarsystem.fandom.com/wiki/47171_Lempo%E2%80%93Hiisi?file=Lempo-Hiisi_Celestia.jpg)]
- Demon artwork: from the game *Lempo*, by One Trick Entertainment

---

## Contact

Questions, suggestions, or contributions?  
Please open an issue or contact [margotcs](https://github.com/margotcs) on GitHub.
