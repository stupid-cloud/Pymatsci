Usage
=====

.. _installation:

Installation
------------

To use Lumache, first install it using pip:

.. code-block:: console

   (.venv) $ pip install pymatsci

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``pymatsci.get_random_ingredients()`` function:

.. autofunction:: pymatsci.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`pymatsci.get_random_ingredients`
will raise an exception.

.. autoexception:: pymatsci.InvalidKindError

For example:

>>> import pymatsci
>>> pymatsci.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

