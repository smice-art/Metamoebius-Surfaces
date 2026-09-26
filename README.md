# Metamoebius Surfaces

A Blender add-on that generates **one-sided (non-orientable) spanning surfaces** from symmetric knot and link diagrams — Möbius bands, and a large family of related "checkerboard surfaces" with holes, twisted bands, and woven membranes.

<p align="center"><i>N-panel generator · live preview · 12 built-in presets · plane or spherical canvas</i></p>

---

## Table of Contents

- [How it works](#how-it-works)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Presets](#presets)
- [Panel reference](#panel-reference)
  - [Top bar](#top-bar)
  - [Diagram](#diagram)
  - [Surface](#surface)
  - [Output](#output)
- [Reading the shape](#reading-the-shape)
- [Tips & troubleshooting](#tips--troubleshooting)
- [Background](#background)

---

## How it works

The construction is the classic **checkerboard (Tait) surface** of a knot or link diagram, extended with a few artistic controls:

1. A closed, symmetric planar curve (or several) is generated — the **diagram**.
2. Its regions are 2-coloured like a chessboard, using the even-odd winding rule.
3. The regions of one colour are triangulated into **membranes**.
4. A disk is cut out around every crossing; inside it, a **half-twisted band** connects the two membrane corners — the over strand bulges up, the under strand dips down.
5. Membrane heights are solved as a smooth (harmonic) field, so everything joins without creases, and can optionally be relaxed towards a true **soap-film / minimal surface**.

The result is a genuine spanning surface of the diagram. For alternating diagrams, one of the two checkerboard colourings is almost always **non-orientable** — a single continuous surface with only one side, in the same family as the Möbius strip but with far richer topology (multiple holes, multiple crossings, "crosscaps" instead of a single half-twist).

This is a general, mathematically documented construction (Tait graphs, checkerboard colourings, minimal-surface relaxation) rather than a reproduction of any single artist's specific, unpublished process — it is built to explore the same *visual family* of woven one-sided surfaces from first principles.

---

## Installation

1. Open Blender (4.5 LTS or newer recommended).
2. `Edit ▸ Preferences ▸ Add-ons ▸ Install from Disk…`
3. Select `metamobius_surfaces.py` and enable the checkbox.
4. Open the **N-panel** in the 3D Viewport (press `N`) and select the **Metamobius** tab.

Alternatively, paste the script into Blender's **Text Editor** and click *Run Script* — this registers the add-on for the current session without installing it.

---

## Quick start

1. Pick a **Surface** from the preset dropdown (defaults to *Möbius Trefoil*).
2. Click **Generate Surface**.
3. Turn on **Live Update** and drag any slider — the mesh regenerates automatically.
4. Switch **Colouring** between *A* and *B* to see the two complementary surfaces of the same diagram — one is often one-sided, the other two-sided.

---

## Presets

| Preset | Diagram | What you get |
|---|---|---|
| **Möbius Trefoil** | Rosette (2,3) | The textbook case: a Möbius band with 3 half-twists |
| **Tri-Web** | Epitrochoid | 3-fold web with three holes and twisted corners |
| **Tri-Web Swirl** | Epitrochoid (3 terms) | Tri-Web with a chiral swirl from the Shape Phase |
| **Tri-Weave Bowl** | Epitrochoid, sphere | 3-fold interwoven bands on a spherical bowl |
| **Borromean Web** | Rings | Three alternating rings, one-sided spanning surface |
| **Hexa Star** | Epitrochoid | 6-fold star web with diamond twists |
| **Hexa Shell** | Epitrochoid, sphere | 6-fold lace wrapped on a sphere |
| **Quad Lace Shell** | Epitrochoid, sphere | 4-fold lace with nested holes |
| **Penta Weave Ball** | Rosette, sphere | 5-lobed Turk's-head weave |
| **Tetra Ribbon Ball** | Rosette, sphere | Ribbons braided around a ball |
| **Pierced Disc** | Rosette, sphere | Dense ring of holes on a dome (two-sided) |
| **Lissajous Lattice** | Lissajous, sphere | 3:4 rectangular weave, 11 crosscaps |
| **Ring Chain** | Chain | A chain of rings — a two-sided counter-example |
| **Custom** | — | Leaves every parameter exactly as set |

Picking a preset only *sets* the parameters below — feel free to keep turning knobs afterwards.

---

## Panel reference

### Top bar

| Control | Description |
|---|---|
| **Surface** | The preset dropdown. Applies a bundle of parameters below in one click. |
| **Generate Surface** | Builds (or rebuilds) the mesh with the current settings. |
| **Live Update** | When on, every slider change regenerates the surface automatically. Turn it off while scrubbing large parameters (e.g. *Lobes*, *Mesh Density*) if it feels sluggish. |

### Diagram

Controls the 2D curve the surface is built from.

| Control | Description |
|---|---|
| **Diagram** | The curve family — see below. |
| **Rotation** | Rotates the diagram in-plane. Purely cosmetic, does not change the shape or its topology. |

**Rosette / Turk's Head** — `r = 1 + A·cos(L·t)`, angle `= W·t`

| Control | Description |
|---|---|
| **Windings W** | How many times the curve winds around the centre. |
| **Lobes L** | Number of petals/lobes. Together with *W*, this sets the crossing count and symmetry. |
| **Amplitude** | Depth of the petals. Small values give near-circular loops; larger values give deep, narrow lobes. |

**Epitrochoid / Farris Wheel** — `z(t) = e^(i·f1·t) + a2·e^(i·(f2·t)) + a3·e^(i·(f3·t + Shape Phase))`

| Control | Description |
|---|---|
| **Freq 1 / 2 / 3** | Rotation speeds of the three summed circles. `Freq 2 − Freq 1` sets the base symmetry. |
| **Amp 2** | Size of the second circle relative to the first. |
| **Amp 3** | Size of the third circle. **Set this above 0 to unlock the Shape Phase.** |
| **Shape Phase** | Only visible/active once *Amp 3 > 0*. With three terms, exactly one phase actually changes the diagram's shape (the other is a pure rotation and has been removed from the panel). Breaks mirror symmetry and adds a chiral swirl. |

**Lissajous** — `x = sin(a·t + φ), y = aspect · sin(b·t)`

| Control | Description |
|---|---|
| **Freq X / Freq Y** | The two Lissajous frequencies — their ratio sets the pattern. |
| **Phase** | Shifts the crossing pattern; unlike the epitrochoid case, this genuinely changes the diagram here. |
| **Aspect Y** | Vertical stretch of the pattern. |

**Ring Link (Borromean …)**

| Control | Description |
|---|---|
| **Rings** | Number of rings arranged around a circle. |
| **Center Distance** | How far each ring's centre sits from the middle. |
| **Ring Radius** | Radius of each individual ring. |
| **Ring Squash** | Vertical squash of each ring (1 = perfect circle). |

**Ring Chain**

| Control | Description |
|---|---|
| **Rings** | Number of rings in the row. |
| **Spacing** | Distance between neighbouring ring centres — controls how deeply they interlock. |
| **Ring Squash** | Vertical squash of each ring. |

### Surface

Controls how the diagram becomes a 3D surface.

| Control | Description |
|---|---|
| **Colouring (A / B)** | Which checkerboard colour becomes the membranes. The two options give the diagram's two complementary spanning surfaces — usually one one-sided, one two-sided. |
| **Canvas: Plane / Sphere** | Draw the diagram on a flat plane, or wrap it onto a sphere (via stereographic projection) for ball/shell/bowl shapes. |
| **Sphere Coverage** | *(Sphere only)* How much of the sphere the diagram reaches around — small values give a shallow bowl, large values a nearly closed ball. |
| **Crossings** | How each crossing's over/under is decided: **Alternating** (standard weave, all bands twist the same way), **Torus Knot** (uses the real over/under sequence of a (W, L) torus knot — Rosette diagrams only), or **Random** (seeded). |
| **Mirror Twist** | Flips the handedness of every crossing at once. |
| **Seed** | *(Random crossings only)* Changes which random over/under pattern is used. |
| **Band Size** | Radius of the twisted disk cut around each crossing, relative to the spacing between crossings. |
| **Twist Height** | How far the over/under bands rise and dip at each crossing. Also the main control for how "wavy" a Soap Film relaxation will look. |
| **Band Resolution** | Number of segments along each twisted band arc — raise for smoother twists, especially at high Twist Height. |
| **Inflate Membranes** | Puffs the flat membranes up (positive) or down (negative), like a cushion. |
| **Edge Smoothing** | Smooths the knot boundary curve itself — higher values round off sharp diagram corners. |
| **Soap Film** | Blends the surface towards a true discrete **minimal surface** (soap film) spanning the same knot, using a cotangent-Laplacian solver. `0` = original surface, `1` = fully relaxed. Most visible on spherical presets or with a high Twist Height — flatter presets are already close to minimal and barely move. |
| **Soap Iterations** | *(Soap Film > 0 only)* Solver passes. The default of 4 is enough for nearly all shapes; raise it only if you push Soap Film to 1 on a very fine mesh. |
| **Mesh Density** | Overall resolution of the interior triangulation. Raise it for cleaner membranes on complex diagrams; lower it for fast iteration. |

### Output

| Control | Description |
|---|---|
| **Size** | Uniform scale of the final object. |
| **Solidify** | Builds a closed, seamless thick shell directly (an *orientation double cover*), rather than using Blender's Solidify modifier. This is what makes one-sided surfaces solidify cleanly with **no seam, no face-orientation errors, and no cracks** at the twisted bands. |
| **Thickness** | *(Solidify only)* Shell thickness. |
| **Subdivision** | Adds a Subdivision Surface modifier for extra smoothness. |
| **Shade Smooth** | Smooth-shades the mesh. On a thin (non-solidified) one-sided surface, custom split normals are used so the unavoidable orientation seam doesn't show up as a dark crease. |
| **Edge Outline** | Draws the knot's boundary curve as a beveled tube alongside the surface (0 disables it). |
| **Preview Material** | Applies a simple Fresnel-driven preview material (pale blue → white) so the surface reads well in the viewport. |

---

## Reading the shape

Every generated surface is either:

- **One-sided (non-orientable)** — the surface has a single continuous side, like a Möbius strip. It has one or more **crosscaps** instead of holes-with-two-sides.
- **Two-sided (orientable)** — an ordinary surface with an inside and an outside, described by its **genus** (number of handles).

Both cases also report:

- **χ (Euler characteristic)** — `vertices − edges + faces` of the underlying combinatorial surface.
- **Boundary loops** — how many separate closed knot curves bound the surface (usually 1, but link diagrams like *Borromean Web* or *Ring Chain* have several).
- **Crossings** — how many twisted bands the surface contains.

This classification is computed every time you click *Generate*. It isn't shown in the panel itself (to keep the panel's height from jumping around as the text changes length), but it's available in Blender's own **Info log** for errors, and programmatically via:

```python
bpy.context.scene.metamobius.info
```

in the Python Console, right after generating.

---

## Tips & troubleshooting

- **"Surface broke up (N boundary loops for M knot component(s))"** — the twisted bands overlapped and split the mesh. Lower **Band Size** or raise **Mesh Density**.
- **A crossing looks pinched or the bands overlap** — reduce **Band Size**, or lower **Twist Height**.
- **Soap Film seems to do nothing** — flat, low-twist presets (e.g. *Tri-Web*) are already close to a minimal surface. Try it on a *Sphere* canvas preset, or raise **Twist Height** first.
- **Shape Phase is greyed out / missing** — it only affects the Epitrochoid family, and only once **Amp 3 > 0**.
- **Live Update feels slow** — turn it off while dragging *Lobes*, *Windings*, *Mesh Density*, or *Soap Iterations*, then click **Generate Surface** once you've settled on a value.
- **Colouring B fails or looks wrong on a Plane canvas** — some diagrams have no bounded "outer" region to invert; the add-on automatically switches that colouring to the *Sphere* canvas when needed.

---

## Background

- **Tait / checkerboard surfaces** — the standard way to build a spanning surface from a knot diagram by 2-colouring its regions.
- **Crosscaps, genus, Euler characteristic** — standard invariants of a closed(-with-boundary) surface, used here to classify each generated shape.
- **Rotation-minimizing frames** *(used internally for the tubular knot outline)* — Wang, Jüttler, Zheng & Liu, *Computation of Rotation Minimizing Frames*, 2008.
- **Discrete minimal surfaces** *(Soap Film)* — Pinkall & Polthier's cotangent-Laplacian iteration for computing discrete minimal surfaces with a fixed boundary.

---

*Generated with the help of Claude.*
