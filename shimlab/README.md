# ShimLab

A shim-stack valving simulator for motorcycle forks and shocks. Open `index.html` in a browser.

- Edit the damper, piston face, shim stack, bleed and oil inputs; force, shim lift, bleed share and shim stress update live.
- Set a baseline to overlay and compare a change.
- Save setups in the browser, or copy/paste them as text.

`engine.js` is the physics model (also inlined in `index.html`). The shim stack is solved as an
axisymmetric plate with radius-varying stiffness (clamped inside, free outside); flow uses the
orifice equation through the bleed in parallel with ports and the lift curtain in series.
The plate solver matches the closed-form solution for a uniform annular plate to within 0.001%.
