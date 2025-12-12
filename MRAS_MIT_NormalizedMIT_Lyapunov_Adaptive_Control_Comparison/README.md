## Model Reference Adaptive Control: MIT and Lyapunov Methods

This report presents a detailed simulation study of **Model Reference Adaptive Control (MRAS)**
applied to a second-order dynamic system. Three adaptive control strategies are derived,
implemented, and compared:

- MIT gradient descent method
- Normalized MIT gradient descent
- Lyapunov-based MRAS

All methods are implemented in continuous time and evaluated using extensive Simulink
simulations.

### Topics Covered
- Continuous-time plant and reference model formulation
- MIT rule derivation using gradient descent
- Implementation of MIT adaptation laws
- Normalized MIT method with α-regularization
- Lyapunov-based MRAS derivation using Lyapunov stability theory
- Stability analysis based on boundedness and convergence
- Simulink implementation of all adaptive controllers
- Parameter convergence analysis
- Control effort evaluation
- Reference tracking for sine and square wave inputs
<img width="875" height="573" alt="image" src="https://github.com/user-attachments/assets/1e8087f1-73b0-4ab7-9476-40507ad7c15e" />

### MIT and Normalized MIT Analysis
The MIT method is analyzed for different values of the adaptation gain λ.
Increasing λ improves parameter convergence, reduces overshoot, and enhances tracking
performance, while maintaining reasonable control effort.

The normalized MIT method introduces an additional normalization factor α.
Results show that increasing λ degrades performance for normalized MIT, while increasing α
improves stability, tracking accuracy, and control effort behavior.

<img width="850" height="521" alt="image" src="https://github.com/user-attachments/assets/000dc1e8-87bd-477f-9074-8af5ccb44764" />

### Lyapunov MRAS Analysis
The Lyapunov-based MRAS is derived using a quadratic Lyapunov function and guarantees
error convergence without requiring parameter convergence.

Simulation results demonstrate:
- Excellent tracking performance across all λ values
- Superior robustness compared to MIT-based methods
- Significantly higher-frequency control effort oscillations

While the magnitude of the control effort remains acceptable, its high-frequency nature may
pose limitations for real-world actuators.

### Input Signal Comparison
Both sine and square wave reference inputs are evaluated. Results indicate that the type of
reference signal has minimal impact on the relative performance of the adaptive methods,
with consistent trends observed across all simulations.

### Key Observations
- Increasing λ improves performance for the MIT method
- Normalized MIT is sensitive to λ but benefits from higher α values
- Lyapunov MRAS provides the best tracking performance overall
- Lyapunov MRAS control effort exhibits high-frequency oscillations
- Parameter convergence is not required for error convergence in Lyapunov MRAS
