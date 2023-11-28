%%%%%%%%%%%%%
Documentation
%%%%%%%%%%%%%

.. contents::
    :local:
    :backlinks: none

--------
Overview
--------

Users may find it useful to read the original manuscript introducing |pj| (Mendes and Landis, Journal, 2024).
That paper summarizes the main features of |pj| while providing an account on the motivations underlying its design and development. 

Graphical models
=============================

|pj| follows the graphical model paradigm popularized in the last decade by programs like RevBayes, BEAST and BEAST 2.
In what follows, we shall offer just a brief refresher on what probabilistic graphical models are, and one simple example of their use in evolutionary modeling.
Those interested in a more detailed exposition will find it in Höhna et al. (2014) and references therein (also, take a look `here <https://revbayes.github.io/tutorials/intro/graph_models.html>`_).

Just as with the RevBayes and BEAST programs, specifying a model in |pj| amounts to writing down any number of conditional dependencies among random variables.
When looked at collectively, random variables comprise a type of probabilistic graphical model called a Bayesian network.
As the name suggest, such models can be described as graphs -- directed acyclic graphs (**DAG**), more specifically -- where each node (or vertex) represents a random variable, and each edge a conditional dependency between the two nodes (variables) it connects.

A Bayesian network representing a very simple model can be seen in figure 1a:

..  figure:: ../images/simple_graphical_model_both_manual.png
    :figwidth: 60%
    :align: center

    **Figure 1.** A simple probabilistic model with :math:`\theta=\{sd, m\}` as represented by a (a) Bayesian network, and a (b) factor graph.

If we were to pen down what this Bayesian network represents in probability terms, we would arrive at a joint density (assuming continuous variables) given by expression:

.. math::

    f_{\Theta,D}(\theta,d) & = f_{D|\Theta}(d|\Theta=\theta)f_\Theta(\theta)\\
    & = f_{\Theta|D}(\theta|D=d)f_D(d),

where :math:`\theta=\{sd, m\}`.

By comparing the Bayesian network graph and the expression above, the attentive reader will notice that functions are not represented in Bayesian networks.
We do not know exactly how random variables are related, other than :math:`\boldsymbol{d}` depending somehow on :math:`sd` and :math:`m`.
This is why it can be helpful to adopt a more general graphical representation: a factor graph (Fig. 1b).

The model underlying the factor graph in figure 1b is the same as that shown by the Bayesian network, the main difference being the additional factor nodes (the filled squares).
Factor nodes can make the relationships between two or more random variables more explicit, especially if their (factor) functions are annotated.

In the example above we can see a factor whose function :math:`f_{D|\Theta}(d|\Theta=\theta)` gives us the probability of random variable :math:`\boldsymbol{d} = \{d_i: 1 \leq i \leq n\}` given :math:`\theta`.
It is annotated as "Normal", so we know it is the probability density function of a normal distribution.
Random variables :math:`sd` and :math:`m` stand for the standard deviation and the mean of that normal distribution.
(For those curious to see just how complicated factor graphs can be, an example can be found in Zhang et al., 2023; see their Supplementary Fig. X.)

How should we read a DAG?
-------------------------

Random variables whose values we do not know (but would like to estimate) are referred to as the **parameters** of the model.
Parameter estimation requires that we know the value(s) of one or more random variables, which we then naturally refer to as **data**.

(transition to Bayes theorem)

.. math::

    f_{\Theta|D}(\theta|D=d) = \frac{f_{D|\Theta}(d|\Theta=\theta)f_\Theta(\theta)}{f_D(d)}

But because unlike those programs |pj| is first and foremost a simulator, 
Some of these variables may be assumed known, in which case they take constant values, while others are 



------------------------------
Graphical user interface (GUI)
------------------------------

----------------------------
Command-line interface (CLI)
----------------------------


-------
Lexicon
-------

.. toctree::
    :maxdepth: 2

    parametric.rst