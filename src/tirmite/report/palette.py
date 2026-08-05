"""
Colours for the report.

One palette serves the whole report so that a pairing group keeps the same
colour in the track browser, the legend, the alignment panels and the static
figures. Colours are assigned by position, so a given run always produces the
same colouring.
"""

from typing import Dict, List, Sequence

__all__ = ['GROUP_COLOURS', 'NEUTRAL_GREY', 'group_colours']

# Paul Tol's bright qualitative scheme. Chosen because it stays distinguishable
# under deuteranopia, protanopia and tritanopia, and because the hues survive
# greyscale printing -- these reports get screenshotted into manuscripts.
GROUP_COLOURS: Sequence[str] = (
    '#4477AA',  # blue
    '#EE6677',  # red
    '#228833',  # green
    '#CCBB44',  # yellow
    '#66CCEE',  # cyan
    '#AA3377',  # purple
    '#BBBBBB',  # grey
    '#EE7733',  # orange
    '#0077BB',  # deep blue
    '#009988',  # teal
)

# Used for model-coverage padding in alignment panels and for contig-end caps:
# deliberately outside the categorical palette so it never reads as a group.
NEUTRAL_GREY = '#9A9A9A'


def group_colours(group_ids: Sequence[str]) -> Dict[str, str]:
    """
    Assign a colour to each pairing group.

    Parameters
    ----------
    group_ids : sequence of str
        Group identifiers in the order they should be coloured, normally the
        order the pairing map lists them.

    Returns
    -------
    dict
        Mapping of group id to a hex colour string.

    Notes
    -----
    With more groups than palette entries the palette repeats. That is
    preferable to generating arbitrary intermediate hues, which quickly stop
    being colour-blind safe; the report additionally distinguishes repeated
    colours by label.
    """
    colours: Dict[str, str] = {}
    seen: List[str] = []
    for group_id in group_ids:
        if group_id in colours:
            continue
        colours[group_id] = GROUP_COLOURS[len(seen) % len(GROUP_COLOURS)]
        seen.append(group_id)
    return colours
