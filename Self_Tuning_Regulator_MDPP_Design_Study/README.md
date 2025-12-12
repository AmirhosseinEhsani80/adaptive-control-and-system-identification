## Self-Tuning Regulator Design Using MDPP

This report presents an extensive simulation study on **Self-Tuning Regulator (STR) design**
based on **Model-Based Desired Pole Placement (MDPP)** for discrete-time and continuous-time systems.

The work systematically analyzes **direct and indirect STR structures**, with and without
pole–zero cancellation, and evaluates their performance under different system properties,
model orders, and disturbance conditions.

<img width="1464" height="800" alt="image" src="https://github.com/user-attachments/assets/1e1dc603-08fa-4465-8be4-10d543cb0cf5" />


### Topics Covered
- MDPP-based controller design using dynamic feedback
- STR design for minimum-phase systems
- STR behavior with and without pole–zero cancellation
- Direct vs indirect STR formulations
- Recursive Least Squares (RLS) parameter estimation
- Over-parameterized and under-parameterized model effects
- STR design under step-shaped disturbances
- STR with and without integrator augmentation
- Long-time simulation and convergence analysis
- STR design for non-minimum-phase systems
- Continuous-time STR using continuous RLS

### Discrete-Time STR Analysis
The report begins with MDPP-based controller design and proceeds to STR implementation
for discrete-time systems. Both **direct** and **indirect** STR approaches are examined,
highlighting their differences in parameter estimation, control effort, and tracking accuracy.

Special emphasis is placed on the role of **pole–zero cancellation**, showing how proper
zero handling significantly improves tracking performance and reduces excessive control effort.

### Model Order Effects
The impact of **over-parameterization** and **under-parameterization** is analyzed in detail.
Results demonstrate that higher-order models may achieve acceptable tracking despite poor
parameter estimation, while under-parameterized models lead to unacceptable tracking
performance and excessive control effort.

### Disturbance Rejection and Integrator Design
The STR controller is evaluated in the presence of step disturbances. An integrator is then
introduced into the controller structure, showing improved tracking performance, reduced
overshoot, and better disturbance rejection, at the cost of slower parameter convergence.

### Long-Time Simulation Behavior
Extended simulations reveal important differences in long-term stability between STR
configurations. While some controllers maintain acceptable control effort over time,
others exhibit divergence or instability in the estimated controller polynomials.

### Non-Minimum Phase and Continuous Systems
The report also investigates STR design for **non-minimum-phase systems**, demonstrating
the limitations of direct STR approaches in such cases. Finally, a **continuous-time STR**
design using continuous RLS and low-pass filtering is presented and evaluated.

### Key Observations
- Indirect STR generally outperforms direct STR in tracking and stability
- Pole–zero cancellation plays a critical role in STR performance
- Over-parameterization is preferable to under-parameterization for control objectives
- Integrator augmentation improves disturbance rejection
- Non-minimum-phase systems severely limit achievable STR performance
- Continuous-time STR is feasible but sensitive to noise and control effort
