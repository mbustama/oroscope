# Logo design brief — Oroscope

Feed the paragraph below to an image-generation session. It is deliberately descriptive
about the physics and deliberately restrictive about the drawing: the failure mode for
this logo is busyness, not blandness.

Destination: `docs/source/_static/oroscope_logo.png`, and **PNG only** — PyPI does not
render SVG, and a project carrying both formats eventually ships two different logos.
Square and large (the current one is 1024×1024 RGBA) so it downscales cleanly, and
transparent outside its disc so it needs no matte. `docs/source/conf.py` points
`html_logo` at it. It also sits at the top of the GitHub README, so it has to survive
being 200 px wide on a white background *and* on a dark one.

---

## The paragraph

> A clean-line logo for **Oroscope**, a tool that reads mountains to find where cosmic
> neutrinos can be caught. The name fuses *oro* — gold, and the Spanish-speaking Andes
> where the search is set — with *-scope*, an instrument for looking, and it deliberately
> echoes *horoscope*: both read the sky, but this one reads it through rock. The image to
> convey is a single physical story. A tau neutrino arrives from below the horizon,
> passes through the body of a mountain where nothing else could, and converts to a tau
> lepton that escapes the far wall of a deep canyon; in the open air above the valley
> floor the tau decays into a shower of particles that fans upward and outward toward a
> detector waiting on the opposite slope. Render this as a minimal geometric mark, not an
> illustration: two angular canyon walls in outline suggesting a V of rock seen in
> cross-section, one straight ray entering low through the left massif and emerging from
> the right wall, and a narrow cone opening from that exit point into the empty space
> between the walls. The upward cone is the only element permitted a flourish — a few
> radiating strokes, or a soft gradient from gold to warm white, evoking both the shower
> and the *oro* of the name. Everything else is a confident single-weight line. Use at
> most two colours besides the line: a deep slate or indigo for the rock, and gold for
> the shower. No text inside the mark, no gridlines, no stars, no detector hardware, no
> ground texture, no drop shadows. It must remain legible as a favicon at 32 px, so the
> ray, the exit point and the cone need to read at a glance, and the whole mark should
> sit comfortably in a square. Flat vector aesthetic, generous negative space, the
> restraint of a good journal figure rather than the density of a mission patch.

---

## Notes for whoever iterates on it

- **The ray must enter low and exit high.** Neutrinos arrive from *below* the horizon;
  a ray drawn coming down from the sky depicts the wrong physics and the wrong tool.
- **The cone opens away from the exit point**, in air, not inside the rock. The decay
  happens in the valley, and the whole point of the canyon is the air gap.
- The gold is the pun. Keep it on the shower, not the rock — gold rock reads as a mining
  logo.
- If a variant is needed for dark backgrounds, invert the rock and keep the gold.
