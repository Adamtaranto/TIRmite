"""
Colours for the report.

One palette serves the whole report, so a pairing group keeps the same colour in
the track browser, the legend, the alignment panels and the static figures.
Colours are assigned by position in a fixed order, never cycled through
generated hues, so a given run always produces the same colouring.

Both modes are selected rather than derived: the dark column is the same eight
hues re-stepped for a dark surface, not an automatic inversion of the light one.

The ordering is the colour-vision-deficiency safety mechanism, not decoration.
The set was validated on the measurable checks -- OKLCH lightness band, chroma
floor, OKLab separation under simulated protanopia and deuteranopia, separation
under normal vision, and WCAG contrast against the chart surface. Worst adjacent
pair: CVD deltaE 9.1 light / 8.4 dark against a target of 8, normal-vision
deltaE 19.6 light / 19.3 dark against a floor of 15. Re-order or substitute a
hue only after re-running those checks.

Three light-mode slots (aqua, yellow, magenta) fall below 3:1 contrast on the
light surface. The relief rule applies and is satisfied here: every group is
named in the legend and in the statistics tables, so identity is never carried
by colour alone.
"""

from typing import Dict, List, Sequence, Tuple

__all__ = [
    'GROUP_COLOURS',
    'GROUP_COLOURS_DARK',
    'NEUTRAL_GREY',
    'NEUTRAL_GREY_DARK',
    'group_colours',
]

# Slot order: blue, orange, aqua, yellow, magenta, green, violet, red.
GROUP_COLOURS: Sequence[str] = (
    '#2a78d6',
    '#eb6834',
    '#1baf7a',
    '#eda100',
    '#e87ba4',
    '#008300',
    '#4a3aa7',
    '#e34948',
)

# The same eight hues stepped for a dark surface.
GROUP_COLOURS_DARK: Sequence[str] = (
    '#3987e5',
    '#d95926',
    '#199e70',
    '#c98500',
    '#d55181',
    '#008300',
    '#9085e9',
    '#e66767',
)

# Used for alignment padding where the model was not covered, and for
# contig-end caps. Deliberately outside the categorical palette so it never
# reads as a pairing group.
NEUTRAL_GREY = '#8a8a86'
NEUTRAL_GREY_DARK = '#6e6e6a'


def group_colours(group_ids: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """
    Assign a colour to each pairing group.

    Parameters
    ----------
    group_ids : sequence of str
        Group identifiers in the order they should be coloured, normally the
        order the pairing map lists them. Repeats are ignored.

    Returns
    -------
    dict
        Mapping of group id to ``(light_hex, dark_hex)``.

    Notes
    -----
    With more groups than palette slots the palette repeats rather than
    generating intermediate hues, which would stop being colour-vision-safe
    almost immediately. Repeated colours stay distinguishable because every
    group is also named in the legend and the statistics tables.
    """
    colours: Dict[str, Tuple[str, str]] = {}
    seen: List[str] = []
    for group_id in group_ids:
        if group_id in colours:
            continue
        slot = len(seen) % len(GROUP_COLOURS)
        colours[group_id] = (GROUP_COLOURS[slot], GROUP_COLOURS_DARK[slot])
        seen.append(group_id)
    return colours
