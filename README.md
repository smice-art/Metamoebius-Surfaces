<p align="center">
  <img src="images/addon.jpg" alt="Metamoebius width="100%">
</p>

# Metamoebius Surfaces (Multi Möbius Generator)

A Blender add-on that generates **one-sided (non-orientable) spanning surfaces** from symmetric knot and link diagrams — Möbius bands, and a large family of related "checkerboard surfaces" with holes, twisted bands, and woven membranes.

# Screen Shot
![Banner Image](images/screen.jpg)

## How it works

The construction is the classic **checkerboard (Tait) surface** of a knot or link diagram, extended with a few artistic controls:

1. A closed, symmetric planar curve (or several) is generated — the **diagram**.
2. Its regions are 2-coloured like a chessboard, using the even-odd winding rule.
3. The regions of one colour are triangulated into **membranes**.
4. A disk is cut out around every crossing; inside it, a **half-twisted band** connects the two membrane corners — the over strand bulges up, the under strand dips down.
5. Membrane heights are solved as a smooth (harmonic) field, so everything joins without creases, and can optionally be relaxed towards a true **soap-film / minimal surface**.

The result is a genuine spanning surface of the diagram. For alternating diagrams, one of the two checkerboard colourings is almost always **non-orientable** — a single continuous surface with only one side, in the same family as the Möbius strip but with far richer topology (multiple holes, multiple crossings, "crosscaps" instead of a single half-twist).

This is a general, mathematically documented construction (Tait graphs, checkerboard colourings, minimal-surface relaxation) rather than a reproduction of any single artist's specific, unpublished process — it is built to explore the same *visual family* of woven one-sided surfaces from first principles.

## Installation

Note: Please download the file located under the Assets section or the latest Releases page, rather than the "Download ZIP" button on the main page.

1. Open Blender (4.5 LTS or newer recommended).
2. Download the latest release
3. In Blender, go to **Edit > Preferences > Add-ons**.
4. Click the dropdown in the top right and select **Install from Disk**.
5. Select your `.zip` file.
4. Open the **N-panel** in the 3D Viewport (press `N`) and select the **Metamobius** tab.

## Quick start

1. Pick a **Surface** from the preset dropdown (defaults to *Möbius Trefoil*).
2. Click **Generate Surface**.
3. Turn on **Live Update** and drag any slider — the mesh regenerates automatically.
4. Switch **Colouring** between *A* and *B* to see the two complementary surfaces of the same diagram — one is often one-sided, the other two-sided.

## Documentation
The Add-on is easy to understand. All adjustment are placed in a N-Panel Slider **Metamobius**. There is a detailed manual in the section DOCS: [MANUAL](docs/handbook.md)

## Views
| Text | Preview | Preview |
| :--- | :--- | :--- |
| Möbius Example | <img src="images/001.jpg" width="250"> |<img src="images/012.jpg" width="250"> |
| Möbius Example | <img src="images/003.jpg" width="250"> |<img src="images/004.jpg" width="250"> |
| Möbius Example | <img src="images/005.jpg" width="250"> |<img src="images/006.jpg" width="250"> |
| Möbius Example | <img src="images/007.jpg" width="250"> |<img src="images/008.jpg" width="250"> |
| Möbius Example | <img src="images/009.jpg" width="250"> |<img src="images/010.jpg" width="250"> |
| Möbius Example | <img src="images/011.jpg" width="250"> |<img src="images/012.jpg" width="250"> |

### Release Notes

### v1.0.0 (September 26, 2026)
- **Publishing**: First public upload of the Add-on.

### Blender
![Blender](https://img.shields.io/badge/Blender-4.3%2B-orange)
![Blender](https://img.shields.io/badge/Blender-4.58-greenorange)
![Blender](https://img.shields.io/badge/Blender-5.0-orange)

* I made this Addon and this small Tutorial with the Help of Claude.
