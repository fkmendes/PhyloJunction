.. PhyloJunction documentation master file, created by
   sphinx-quickstart on Mon Aug  1 09:58:05 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

%%%%%%%%%%%%%
PhyloJunction
%%%%%%%%%%%%%

.. toctree::
   :maxdepth: 3
   :hidden:
   
   Home <self>
   Installation <installing.rst>

In a nutshell, |pj| (PJ) is a pure-|Python3| library whose goals are to:

    (i) facilitate and expedite model development and testing, and
    (ii) make model specification and simulation seamless and reproducible.

In order to meet these goals, |pj| implements an expandable ecosystem of models and parsing functionalities that it exposes to users through command-line (CLI) and graphical user interfaces (GUI).

Interacting with the CLI and GUI consists mainly in specifying a model by entering commands (or loading a `.pj` script) written in the `phylojunction` (lowercase) scripting language.
`phylojunction` is a lightweight language that finds its roots in STAN and JAGS, and is largely inspired by its older cousins from phylogenetic modeling, Rev (RevBayes_) and LPhy ().

Once a model is specified, |pj| then samples (i.e., simulates) values for each of the model variables, parses and summarizes those values, and outputs the result as text-formatted tables.
If interacting with the GUI, users can further inspect graphical summaries of those values.

On this website, you will find instructions and tutorials on how to install, use, and develop |pj|.
Please refer to the relevant sections for help with questions, bug reports, and feature requests.

% Subroutines for building models, simulating data and parsing results are nonetheless functional outside of the standalone user interfaces, and can be accessed via custom scripts or from within the Python interpreter.

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
