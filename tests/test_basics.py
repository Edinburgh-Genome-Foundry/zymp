import os
import matplotlib

matplotlib.use("Agg")

from zymp import (
    stacked_sites_array,
    plot_sequence_sites,
    annotate_enzymes_sites,
    write_record,
)

enzymes_names = [
    "AccI",
    "AclI",
    "AflII",
    "AflIII",
    "AgeI",
    "ApaLI",
    "AseI",
    "AvaI",
    "BamHI",
    "BanII",
    "BlnI",
    "BmtI",
    "BsmI",
    "BssHII",
    "DdeI",
    "DraI",
    "Eco47III",
    "EcoRI",
    "EcoRV",
    "HindII",
    "HindIII",
    "HinfI",
    "HpaI",
    "KpnI",
    "MfeI",
    "MluI",
    "MspA1I",
    "MunI",
    "NaeI",
    "NcoI",
    "NdeI",
    "NheI",
    "NotI",
    "NsiI",
    "NspI",
    "PstI",
    "PvuI",
    "PvuII",
    "SacI",
    "SacII",
    "SalI",
    "ScaI",
    "SfaNI",
    "SnaBI",
    "SpeI",
    "SphI",
    "SspI",
    "StyI",
    "VspI",
    "XhoI",
    "XmaI",
    "ZraI",
]

forbidden_enzymes = ["BsmBI", "BsaI"]


def test_basic_example(tmpdir):

    seq, sites_in_seq, leftover = stacked_sites_array(
        enzymes_names, forbidden_enzymes=forbidden_enzymes, tries=200
    )
    assert len(seq) < 170
    assert len(sites_in_seq) > 43

    # PLOT A SUMMARY
    ax = plot_sequence_sites(seq, enzymes_names)

    # WRITE THE SEQUENCE AND SITE ANNOTATIONS AS A RECORD
    record = annotate_enzymes_sites(
        seq, enzymes_names, forbidden_enzymes=forbidden_enzymes
    )
    write_record(record, os.path.join(str(tmpdir), "test.gb"))


def test_basic_example_with_condition(tmpdir):
    enzymes_names = [
        "AccI",
        "AclI",
        "AflII",
        "AflIII",
        "AgeI",
        "ApaLI",
        "AseI",
        "AvaI",
        "XhoI",
        "SacII",
        "XmaI",
    ]

    for i in range(3):
        for e in ["XhoI", "XmaI"]:

            def success_condition(seq, sites_in_seq, leftover):
                return {e}.issubset(sites_in_seq)

            seq, sites_in_seq, leftover = stacked_sites_array(
                enzymes_names,
                forbidden_enzymes=forbidden_enzymes,
                tries=500,
                success_condition=success_condition,
            )

            assert {e}.issubset(sites_in_seq)
