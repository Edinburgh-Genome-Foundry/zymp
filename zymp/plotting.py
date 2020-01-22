from dna_features_viewer import BiopythonTranslator

from .tools import annotate_enzymes_sites


def plot_sequence_sites(
    sequence,
    enzymes_names,
    forbidden_enzymes=(),
    unique_sites=True,
    ax=None,
    figure_width=18,
    annotate_inline=True,
):
    """Plot the location of sites in the sequence.

    Non-unique and forbidden sites can be highlighted in red.

    Parameters
    ----------

    sequence
      The sequence of interest. ATGC string.

    enzymes_names
      List of names of the enzymes to plot.

    forbidden_enzymes
      The sites of these enzymes will also be plotted, but with a red
      background

    unique_sites
      If true, for each enzyme in enzyme_name with more than one site in
      the sequence, these will be plotted on a red background.

    ax
      Matplotlib ax on which to draw the figure. If none is provided a new
      figure is created and the ax is returned at the end.

    figure_width
      Width of the figure if no ax is provided and a new figure is returned.

    annotate_inline
      If True, the enzyme names will be written inside the annotations
      when possible, instead of above.

    """

    record = annotate_enzymes_sites(
        sequence,
        enzymes_names,
        forbidden_enzymes=forbidden_enzymes,
        unique_sites=unique_sites,
    )
    default_props = dict(
        thickness=10,
        box_color=None,
        fontdict=dict(
            family="Impact", size=7, color="black", weight="normal"
        ),
    )
    translator = BiopythonTranslator(
        features_properties=lambda f: default_props
    )
    graphic_record = translator.translate_record(record)
    graphic_record.labels_spacing = 1
    ax, _ = graphic_record.plot(
        figure_width=figure_width, annotate_inline=annotate_inline, ax=ax
    )
    return ax
