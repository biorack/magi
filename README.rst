Metabolite Atlas Reaction Database
==================================

Metabolomics is the comprehensive profiling of the small molecule composition of a biological sample. This approach is being used to provide new insights into a variety of biological systems (clinical, bioenergy, etc.). A grand challenge for metabolomics is the complexity of the data, which often include many experimental artifacts. This is compounded by the tremendous chemical diversity of metabolites where identification of each uncharacterized metabolite is in many ways its own puzzle. The Metabolite Atlas project will provide easy to use tools that enable a scientist to quickly capture knowledge about what compounds they have observed and to propagate that knowledge to future experiments. Instead of having to sift through billions of data points, a scientist will use pre-existing method specific metabolite atlas results to suggest candidate identifications as soon as the data is generated.


Features
--------
- A database of all publicly available reactions
- Compounds in those reactions
- Atom mapping of products and reactants
- Comparison of all reactions
- API for accessing this information programmaticallly

Local Installation
------------------

.. code-block:: bash

    $ git clone https://github.com/biorack/metatlas_reactions.git
    $ cd metatlas_reactions
    $ pip install .
