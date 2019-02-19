.. raw:: html

    <p align="center">
    <img alt="stacked array" title="stacked array" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/zymp/master/docs/_static/images/title.png" width="300">
    <br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/zymp.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/zymp
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/zymp/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/zymp?branch=master

Zymp is a Python library to produce small sequences of DNA packed with enzyme
restriction sites. You specify the enzymes you want, the ones you don't want,
whether you want the sites to be unique, or any other condition, and Zymp will
attempt to find a compact sequence verifying all of this (it really focuses on
sequence shortness).

**Warning:** Zymp is implemented with a "whatever works well enough"
philosophy. It has a lot of "whatever" but it generally works "well enough".
The algorithm is greedy with many simplifications so don't expect perfect solutions.

Examples
--------

Here is how you design a sequence

.. code:: python

    from zymp import (stacked_sites_array, plot_sequence_sites,
                  annotate_enzymes_sites, write_record)

    enzymes_names = [
        'AccI', 'AclI', 'AflII', 'AflIII', 'AgeI', 'ApaLI', 'AseI',
        'AvaI', 'BamHI', 'BanII', 'BlnI', 'BmtI', 'BsmI', 'BssHII',
        'DdeI', 'DraI', 'Eco47III', 'EcoRI', 'EcoRV', 'HindII',
        'HindIII', 'HinfI', 'HpaI', 'KpnI', 'MfeI', 'MluI',
        'MspA1I', 'MunI', 'NaeI', 'NcoI', 'NdeI', 'NheI', 'NotI',
        'NsiI', 'NspI', 'PstI', 'PvuI', 'PvuII', 'SacI', 'SacII',
        'SalI', 'ScaI', 'SfaNI', 'SnaBI', 'SpeI', 'SphI', 'SspI',
        'StyI', 'VspI', 'XhoI', 'XmaI', 'ZraI'
    ]

    forbidden_enzymes=['BsmBI', 'BsaI']

    # DESIGN AN OPTIMIZED SEQUENCE WITH ZYMP
    seq, sites_in_seq, leftover = stacked_sites_array(
            enzymes_names, forbidden_enzymes=forbidden_enzymes,
            unique_sites=True, tries=100)

    print ("Sequence length:", len(seq),
        "\nRestriction sites:", len(sites_in_seq),
        "\nSites not included: ", leftover)
                    
    # PLOT A SUMMARY
    ax = plot_sequence_sites(seq, enzymes_names)
    ax.figure.savefig("stacked_array.pdf", bbox_inches='tight')
                    
    # WRITE THE SEQUENCE AND SITE ANNOTATIONS AS A RECORD
    record = annotate_enzymes_sites(
        seq, enzymes_names, forbidden_enzymes=forbidden_enzymes)
    write_record(record, 'stacked_site_array.gb')

**Plot output:**

.. raw:: html

    <p align="center">
    <img alt="stacked array" title="stacked array" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/zymp/master/docs/_static/images/example_array.png" width="800">
    <br />
    </p>


**Console output:**

.. code:: bash

    Sequence length: 159
    Restriction sites: 49
    Sites not included:  {'NcoI', 'HpaI', 'SacII'}

Zymp has created a 159-nucleotide sequence with 49 of the 52 restriction sites
we specified, that's only ~3 nucleotides per site ! and the sequence is free
of BsaI or HpaI sites, so it is compatible with Golden Gate assembly.

If NcoI and HpaI are your favorite enzymes, you may be disappointed that they
are not in the final sequence. Zymp allows you to add validity conditions
for the result:

.. code:: python

    from zymp import stacked_sites_array

    def success_condition(seq, sites_in_seq, leftover):
        return {'NcoI', 'HpaI'}.issubset(sites_in_seq)

    seq, sites_in_seq, leftover = stacked_sites_array(
            enzymes_names, forbidden_enzymes=forbidden_enzymes,
            tries=100, success_condition=success_condition)

    print ("Sequence length:", len(seq),
        "\nRestriction sites:", len(sites_in_seq),
        "\nSites not included: ", leftover)

**New console output:**

.. code:: bash

    Sequence length: 158 
    Restriction sites: 47 
    Sites not included:  {'SacII', 'SacI', 'XhoI', 'BlnI', 'XmaI'}


Installation
-------------

You can install zymp through PIP

.. code::

    sudo pip install zymp

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

License = MIT
--------------

Zymp is an open-source software originally written at the
`Edinburgh Genome Foundry <http://genomefoundry.org>`_ by
`Zulko <https://github.com/Zulko>`_ and
`released on Github <https://github.com/Edinburgh-Genome-Foundry/zymp>`_
under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

Zymp is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
