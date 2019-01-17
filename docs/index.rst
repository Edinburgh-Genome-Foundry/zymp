.. raw:: html

    <p align="center">
    <img alt="zymp Logo" title="zymp Logo"
    src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/zymp/master/docs/_static/images/title.png" width="600">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/zymp.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/zymp
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/zymp/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/zymp?branch=master


Zymp is a Python library to use the Kappa biological modeling system.
It is built on top of ``kappy``, the official Python binding for Kappa, and
provides a more *pythonic* interface to define models and analyse results:
Python objects to define agents and rules, methods for visualizing complexes and
plotting time series, pretty error printing, etc.

Example
--------

Here is a standard Kappa script defining a situation with 3 agents A, B, C which
can irreversibly bind together to form A-B-C chains. We end the simulation after
10 (virtual) seconds and take a snapshot.

.. code::

    %agent: A(a, b)
    %agent: B(b, c)
    %agent: C(c, d)

    'a.b' A(b[.]),B(b[.]) -> A(b[1]),B(b[1]) @ 0.005
    'b.c' B(c[.]),C(c[.]) -> B(c[1]),C(c[1]) @ 0.021

    %init: 100 A()
    %init: 100 B()
    %init: 100 C()

    %mod: alarm 10.000 do $SNAPSHOT "end";

    %plot: |B(b[.])|
    %plot: |C(c[.])|

Now here is the same model defined in Python with Zymp:

.. code:: python

    from zymp import KappaModel, KappaAgent, KappaRule, KappaSiteState)

    model = KappaModel(
        agents = [
            KappaAgent('A', ('a', 'b')),
            KappaAgent('B', ('b', 'c')),
            KappaAgent('C', ('c', 'd'))
        ],
        rules = [
            KappaRule(
                'a.b',
                [
                    KappaSiteState('A', 'b', '.'),
                    KappaSiteState('B', 'b', '.')
                ],
                '->',
                [
                    KappaSiteState('A', 'b', '1'),
                    KappaSiteState('B', 'b', '1')
                ],
                rate=0.5e-2
            ),
            KappaRule(
                'b.c',
                [
                    KappaSiteState('B', 'c', '.'),
                    KappaSiteState('C', 'c', '.')
                ],
                '->',
                [
                    KappaSiteState('B', 'c', '1'),
                    KappaSiteState('C', 'c', '1')
                ],
                rate=2.1e-2
            )
        ],
        initial_quantities={
            'A': 100,
            'B': 100,
            'C': 100,
        },
        duration=10,
        snapshot_times={
            'end': 10
        },
        plots=[KappaSiteState('B', 'b', '.'),
               KappaSiteState('C', 'c', '.')]
    )

    # Now simulate the model !
    simulation_results = model.get_simulation_results()

The Python version is admittedly longer than the original Kappa script, but
where Python shines is in its ability to programmatically define much more
complex and dynamic models with dozens of agents and rules.

Zymp also makes it easy to simulate the result and get the final data, and
provide a few utilities to vizualize the results. For instance let us plot
the two time series records during the simulation:

.. code:: python
    
    from zymp import plot_simulation_time_series
    ax = plot_simulation_time_series(simulation_results['plots'])
    ax.figure.savefig('basic_example_time_series.png')

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/zymp/master/examples/basic_example_time_series.png" width="640">
    </p> 

And here is how you plot the products present at the end of the simulation:

.. code:: python

    end_agents = simulation_results['snapshots']['end']['snapshot_agents']
    fig, axes = plot_snapshot_agents(end_agents)
    fig.savefig('basic_example_agents_graphs.png')

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/zymp/master/examples/basic_example_agents_graphs.png" width="640">
    </p> 

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

Zymp is an open-source software originally written at the `Edinburgh Genome Foundry <http://genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_ and `released on Github <https://github.com/Edinburgh-Genome-Foundry/zymp>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

Zymp is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.


.. toctree::
    :hidden:
    :maxdepth: 3

    self


.. toctree::
   :hidden:
   :caption: Reference
   :maxdepth: 3

   ref


.. _Github: https://github.com/EdinburghGenomeFoundry/dnachisel
.. _PYPI: https://pypi.python.org/pypi/dnachisel

