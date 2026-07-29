# V7.10.0 model audit

## Condition

| Quantity | Value | Role |
|---|---:|---|
| Vesicle geometry | Fig. 7 contract | Defines the reference contour and areas |
| Mesh | 321 nodes smoke; 641 nodes production | Axisymmetric meridian discretization |
| Streptavidin concentration | 10-80 nM | Changes preferred inner-Ld leaflet area |
| Streptavidin saturation concentration | 80 nM | Occupancy mapping |
| Inner-Ld biotin mole fraction | 0.002 | Enters the bound-site area shift |
| Bound-site area coupling, gamma_A | 4.0 | Amplifies the occupancy-derived preferred-area shift |
| Lo bending rigidity | 80 kBT | Active energy parameter |
| Ld bending rigidity | 25 kBT | Active energy parameter |
| Lo/Ld line tension | 0.3 pN | Active boundary energy parameter |
| Gaussian-rigidity ratio | 0.9 | Active boundary energy parameter |
| Combined Ld spontaneous curvature at 10 nM | -0.015 um^-1 | Fixed through the concentration sequence |
| Post-CTB Lo spontaneous curvature | 0.3333333333 um^-1 | Fixed through the concentration sequence |
| CTB outer-Lo area-expansion input | 0.0916 | Produces the fixed post-CTB Lo area term |
| Experimental dome-increase target | 2.95 um | Validation only; not an energy term |

## Terms entering the minimized energy

| Term | Domain or interface | Active |
|---|---|---|
| Helfrich bending energy, `(2H-C0)^2` | Lo and Ld separately | Yes |
| Domain area-strain energy | Lo and Ld separately | Yes |
| Leaflet area-difference energy | Lo and Ld separately | Yes |
| Lo/Ld line energy | Shared boundary | Yes |
| Gaussian-rigidity boundary energy | Shared boundary | Yes |
| Fixed enclosed-volume constraint | Whole vesicle | Yes |
| Streptavidin-dependent preferred inner-Ld area | Ld | Yes |
| Fixed CTB-dependent Lo preferred area and mismatch | Lo | Yes |

The natural interfacial force and bending-moment balances are consequences of
minimizing this shared two-domain energy. They are evaluated as residual
diagnostics and are not added as duplicate fitted penalties.

## Inactive or diagnostic-only quantities

| Quantity | Status |
|---|---|
| Direct phenomenological Lo-Ld coupling | Inactive |
| Explicit protein binding or work energy | Inactive |
| Target-shape force | Inactive |
| Concentration-dependent spontaneous curvature | Inactive; both domain C0 values are fixed |
| Experimental dome target and tolerance | Validation only |
| Optimizer callback frames | Numerical diagnostics only |

## Success gates

The production validator checks numerical acceptance, geometry validity,
constraint errors, concentration ordering, monotonic dome suppression,
restoration at high concentration, interface residuals, and availability of
the packaged tables, meshes, images, JSON trajectory, and MP4 files.
