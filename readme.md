# CertifiedMonodromyComputation.jl

**CertifiedMonodromyComputation.jl** is a Julia package for **certified numerical homotopy tracking** and **monodromy group computation**. It uses interval arithmetic (via [Nemo.jl](https://github.com/Nemocas/Nemo.jl) and Arb) to provide mathematically rigorous results and integrates with [GAP](https://www.gap-system.org/) for group-theoretic analysis.

## Features

* **Certified Tracking:** Uses interval arithmetic and the Krawczyk test to certify solution paths.
* **Automated Monodromy:** Automatically tracks the complete graph to generate the monodromy group (`solve_monodromy`).


## Quick Start

### 1. Basic certified tracking

If you want to track a single solution path from $t=0$ to $t=1$:

```julia
using Nemo, AbstractAlgebra
using CertifiedMonodromyComputation

# 1. Setup the System
@monodromy_setup begin
    vars = (x, y)
end
const CCi = _CCi # Alias for the coefficient ring (Complex Interval Field)

# 2. Define your system F(x, y)
f1 = x^2 + 3*y - 4
f2 = y^2 + 3

F = [f1 f2]
G = [x^2-1 y^2-1]

# 3. Define the start point at t=0 and the homotopy H
H = straight_line_homotopy(F, G, t)
point = [CCi(1), CCi(-1)]

# 4. Track!
track(H, point)
track(H, point; iterations_count=true) # print the number of iterations
track(H, point; show_display=false) # turn off the display

```
