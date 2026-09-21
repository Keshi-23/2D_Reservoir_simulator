# 2-D Reservoir Simulator

## Numerical Simulation of Immiscible Two-Phase Flow in Porous Media with Capillary, Gravity, and Heterogeneity Effects

This repository contains a **2-D reservoir simulator** for numerical modeling of **immiscible two-phase flow in porous media**. The simulator is based on the **IMPES (Implicit Pressure, Explicit Saturation)** formulation and is designed for black-oil reservoir applications.

The model accounts for important reservoir-flow mechanisms, including:

- Pressure-driven multiphase flow
- Capillary pressure
- Gravity segregation
- Spatially heterogeneous rock properties
- Relative permeability effects
- Fluid mobility variations
- Injection and production wells
- Pressure and saturation evolution
- Well pressure and flow-rate performance

The simulator provides spatial maps of reservoir pressure and phase saturation together with well-performance plots for evaluating pressure and production/injection behavior through time.

---

## 1. Project Objectives

The main objective of this project is to provide a numerical framework for studying two-phase flow in heterogeneous porous media.

The simulator can be used to investigate the effects of:

- Reservoir heterogeneity
- Permeability and porosity distributions
- Fluid viscosity and density
- Relative permeability
- Capillary pressure
- Gravity
- Initial pressure and saturation
- Well locations
- Injection and production conditions
- Reservoir geometry

The code is organized so that the **numerical solver**, **input data**, and **simulation results** remain separated.

---

## 2. Numerical Method

The simulator uses the **IMPES** approach:

> **Implicit Pressure, Explicit Saturation**

At every time step, the numerical procedure can be summarized as:

1. Read the current pressure and saturation fields.
2. Evaluate fluid and rock-dependent properties.
3. Calculate phase relative permeabilities and mobilities.
4. Evaluate capillary-pressure and gravity contributions.
5. Assemble the pressure equation.
6. Solve the pressure equation implicitly.
7. Calculate phase fluxes from the updated pressure field.
8. Update saturation explicitly.
9. Apply well source/sink terms.
10. Advance to the next time step.
11. Save pressure, saturation, and well-performance results.

## 3. Future Development

Possible extensions of the simulator include:

- Adaptive time stepping
- Fully implicit formulation
- Three-phase flow
- Gas-oil-water black-oil formulation
- Compressible rock and fluids
- More advanced well models
- Multiple injection and production wells
- Unstructured grids
- Local grid refinement
- Anisotropic permeability
- Improved nonlinear relative-permeability models
- Advanced capillary-pressure models
- Reservoir history matching
- Parallel computing
- Coupling with geomechanics
- Thermal reservoir simulation

---

## 4. Project Scope

This repository is intended for:

- Reservoir simulation research
- Multiphase-flow studies
- Numerical-method development
- Porous-media flow analysis
- Educational applications
- Testing the influence of capillary pressure, gravity, and heterogeneity on two-phase displacement

---

## 5. Citation

If this repository contributes to your research or publication, please cite the repository and any associated publication when available.

```text
Author: Promise Longe
Repository: https://github.com/Keshi-23/2D_Reservoir_simulator
Method: IMPES numerical simulation of immiscible two-phase flow in porous media
```

---

## 6. License

Please refer to the repository license for terms governing the use, modification, and distribution of this code.

---

## 7. Contact

For questions, suggestions, or research collaboration, please open an issue in this repository or contact the repository author.

---

**Keywords:** reservoir simulation, IMPES, black-oil model, two-phase flow, porous media, capillary pressure, gravity, reservoir heterogeneity, permeability, saturation, pressure, well performance, numerical simulation.
