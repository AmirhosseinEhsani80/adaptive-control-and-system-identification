## Adaptive and Non-Adaptive Model Predictive Control with Delay Estimation

This report presents a comprehensive simulation study on **Model Predictive Control (MPC)**
for discrete-time dynamic systems with **input delay**, considering both **non-adaptive**
and **adaptive** controller formulations.

Multiple MPC strategies are implemented and compared under different noise, disturbance,
delay-mismatch, and system phase conditions.
<img width="872" height="422" alt="image" src="https://github.com/user-attachments/assets/f2788d42-1f55-41a8-99f9-627a5601ae13" />

### Topics Covered
- Discrete-time system modeling with input delay
- Minimum Degree Pole Placement (MDPP) controller design
- Non-adaptive MPC methods:
  - One Step Ahead MPC
  - Weighted One Step Ahead MPC
  - MPC with J2 loss function (non-minimum phase)
  - MPC with J3 loss function (integrator-augmented)
  - Constant Future MPC
- Root locus–based stability analysis
- Effect of delay mismatch between plant and controller
- White noise and disturbance robustness analysis
- Adaptive delay estimation
- Fully adaptive MPC implementations of all methods

### Non-Adaptive MPC Analysis
The report first evaluates MPC methods without system identification. Performance is
analyzed under:
- Noise-free conditions
- White noise
- White noise with disturbance
- Changing and mismatched delays

Results show that weighted MPC reduces control effort but introduces steady-state errors,
while constant-future MPC provides improved tracking but limited noise rejection.

### Loss Function Design (J2 and J3)
For non-minimum-phase systems, J2 and J3 loss functions are investigated.
- J2 improves stability but amplifies disturbance sensitivity
- J3 introduces an integrator to remove steady-state error but significantly increases
noise amplification

Root locus analysis is used to select weighting parameters and assess stability margins.

### Delay Mismatch and Estimation
The effect of incorrect delay modeling is explicitly studied. Results demonstrate that
delay mismatch can severely degrade transient performance and destabilize the system.
An adaptive delay estimation strategy is then introduced, allowing the controller to
update its structure online.

### Adaptive MPC
All MPC methods are extended to **adaptive versions** using online parameter estimation.
Adaptive MPC shows improved tracking and robustness in the presence of unknown delay,
with constant-future MPC exhibiting the best overall performance.

### Key Observations
- Delay mismatch is a dominant source of performance degradation
- Weighted MPC reduces control effort but worsens tracking accuracy
- J2 and J3 loss functions trade noise robustness for steady-state performance
- Integrator augmentation amplifies noise sensitivity
- Adaptive delay estimation significantly improves MPC robustness
- Constant-future MPC achieves the best balance between tracking and control effort
