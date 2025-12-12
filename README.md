# Adaptive Control and System Identification

This repository contains a collection of **simulation-based studies in adaptive control and system identification**, focusing on the design, analysis, and comparison of modern adaptive control strategies under varying system dynamics, noise conditions, delays, and modeling assumptions.

The projects emphasize **theoretical correctness**, **algorithmic clarity**, and **systematic performance evaluation**, with implementations primarily developed in **MATLAB/Simulink**.

---

## Scope and Focus

The work in this repository covers core topics in adaptive control and identification, including:

- Online and offline **system identification**
- Recursive parameter estimation (LS, RLS, Kalman Filter)
- **Model Reference Adaptive Control (MRAC / MRAS)**
- Self-Tuning Regulators (STR)
- Stochastic and minimum variance control
- Model Predictive Control (MPC), including adaptive and delay-aware formulations
- Stability, robustness, and performance trade-offs

Each project investigates how controller structure, estimation method, excitation, noise, and modeling assumptions affect **tracking performance**, **control effort**, and **robustness**.

---

## Project Overview

### Adaptive System Identification
Studies offline and online identification using:
- Least Squares (LS)
- Recursive LS (RLS), λ-RLS, covariance resetting
- Kalman Filter–based estimation  
with analysis of excitation richness, noise color, model order mismatch, and feedback effects.

### Self-Tuning Regulator (STR) Design
Covers direct and indirect STR implementations based on:
- Model-Based Desired Pole Placement (MDPP)
- RLS-based parameter estimation
- Pole–zero cancellation strategies  
including disturbance rejection, integrator augmentation, and long-term behavior.

### STR vs Stochastic Control
Compares adaptive STR controllers with:
- Minimum variance (MV) control
- Moving average (MA) control  
under colored noise and non-minimum-phase dynamics, highlighting variance reduction and robustness limits.

### Model Predictive Control (MPC)
Investigates non-adaptive and adaptive MPC formulations, including:
- One-step-ahead and weighted MPC
- Constant-future MPC
- J2 and J3 loss functions
- Delay mismatch and adaptive delay estimation  
with root-locus–based stability analysis.

### Model Reference Adaptive Control (MRAS)
Provides a detailed comparison of:
- MIT rule (gradient descent)
- Normalized MIT
- Lyapunov-based MRAS  
with emphasis on adaptation gains, parameter convergence, control effort behavior, and reference tracking.

---

## Methodology

Across all projects, the following methodology is consistently applied:

- Clear mathematical formulation and derivation
- Explicit modeling assumptions
- Controlled simulation experiments
- Comparative analysis using:
  - Output tracking
  - Parameter convergence
  - Control effort
  - Robustness to noise and uncertainty

The goal is not only to achieve good tracking, but to **understand why certain adaptive strategies succeed or fail** under specific conditions.

---

## Tools

- MATLAB
- Simulink
- Discrete- and continuous-time modeling
- Numerical simulation and signal analysis

---

## Intended Audience

This repository is intended for:
- Control systems researchers
- Graduate students in control and robotics
- Engineers working on adaptive and data-driven control
- PhD application reviewers seeking evidence of applied control expertise

The material assumes familiarity with linear systems, control theory, and adaptive algorithms.

---

## License

This project is released under the **MIT License**.
