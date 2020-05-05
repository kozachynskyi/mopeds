=======
par_est
=======

This is the documentation of **par_est**.

.. note::

    This is the main page of your project's `Sphinx`_ documentation.
    It is formatted in `reStructuredText`_. Add additional pages
    by creating rst-files in ``docs`` and adding them to the `toctree`_ below.
    Use then `references`_ in order to link them from this page, e.g.
    :ref:`authors` and :ref:`changes`.

    It is also possible to refer to the documentation of other Python packages
    with the `Python domain syntax`_. By default you can reference the
    documentation of `Sphinx`_, `Python`_, `NumPy`_, `SciPy`_, `matplotlib`_,
    `Pandas`_, `Scikit-Learn`_. You can add more by extending the
    ``intersphinx_mapping`` in your Sphinx's ``conf.py``.

    The pretty useful extension `autodoc`_ is activated by default and lets
    you include documentation from docstrings. Docstrings can be written in
    `Google style`_ (recommended!), `NumPy style`_ and `classical style`_.


Jacobian Matrix
===============

Jacobian Matrix from Simulator is provided in a given form (first parameter is tau):

.. math::

        \begin{bmatrix}
        \begin{matrix}
            \frac{\partial{Y_1}^{t=1}}{\partial{p_1}}  & \frac{\partial{Y_1}^{t=1}}{\partial{p_2}}  & ... & \frac{\partial{Y_1}^{t=1}}{\partial{p_{np+1}}}  \\
            \frac{\partial{Y_2}^{t=1}}{\partial{p_1}}  &  \frac{\partial{Y_1}^{t=1}}{\partial{p_2}} & ...  & \frac{\partial{Y_2}^{t=1}}{\partial{p_{np+1}}}  \\
            ...  & ... & ... & ...  \\
            \frac{\partial{Y_{ns}}^{t=1}}{\partial{p_1}} & \frac{\partial{Y_{ns}}^{t=1}}{\partial{p_2}} & ... & \frac{\partial{Y_{ns}}^{t=1}}{\partial{p_{np+1}}}         
        \end{matrix}
        \begin{vmatrix}
            \frac{\partial{Y_1}^{t=2}}{\partial{p_1}}  & \frac{\partial{Y_1}^{t=2}}{\partial{p_2}}  & ... & \frac{\partial{Y_1}^{t=2}}{\partial{p_{np+1}}} \\
            \frac{\partial{Y_2}^{t=2}}{\partial{p_1}}  &  \frac{\partial{Y_1}^{t=2}}{\partial{p_2}} & ...  & \frac{\partial{Y_2}^{t=2}}{\partial{p_{np+1}}} \\
            ...  & ... & ... & ... \\
            \frac{\partial{Y_{ns}}^{t=2}}{\partial{p_1}} & \frac{\partial{Y_{ns}}^{t=2}}{\partial{p_2}} & ... & \frac{\partial{Y_{ns}}^{t=2}}{\partial{p_{np+1}}}
        \end{vmatrix}
        \begin{matrix}
            \frac{\partial{Y_1}^{t=nt}}{\partial{p_1}}  & \frac{\partial{Y_1}^{t=nt}}{\partial{p_2}}  & ... & \frac{\partial{Y_1}^{t=nt}}{\partial{p_{np+1}}} \\
            \frac{\partial{Y_2}^{t=nt}}{\partial{p_1}}  &  \frac{\partial{Y_1}^{t=nt}}{\partial{p_2}} & ...  & \frac{\partial{Y_2}^{t=nt}}{\partial{p_{np+1}}} \\
            ...  & ... & ... & ... \\
            \frac{\partial{Y_{ns}}^{t=nt}}{\partial{p_1}} & \frac{\partial{Y_{ns}}^{t=nt}}{\partial{p_2}} & ... & \frac{\partial{Y_{ns}}^{t=nt}}{\partial{p_{np+1}}}
        \end{matrix}
        \end{bmatrix}

     
Covariance measurement matrix is diagonal matrix in form:

.. math::

        \begin{bmatrix}
            \sigma_1 & 0 & 0 & ... & 0 \\
            0 & \sigma_2 & 0 & ... & 0 \\
            0 & 0 & \sigma_3 & ... & 0 \\
            ... & ... & ... & ...  & 0 \\
            0 & 0 & 0 & 0 & \sigma_{ns} 
        \end{bmatrix}

Covariance matrix is calculated as following:

.. math::

    \left(\cfrac{\partial{Y}^{t}}{\partial{p}}\right)^T \times COV_Y \times \cfrac{\partial{Y}^{t}}{\partial{p}}

which outputs a matrix, where only diagonal sub-matrixes are of interest:

.. math::


      
    \begin{equation*}
    \begin{array}{|c c c c|c c c c|c c c c|}
    \hline
    \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}}\\
    \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}}\\
    \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots \\
    \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}}\\
    \hline
    \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}}\\
    \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}}\\
    \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots \\
    \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}}\\
    \hline
    \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{1}^{t=1}}{\partial p_{np+1}}\\
    \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{2}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{1}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{2}^{t=1}}{\partial p_{np+1}}\\
    \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots \\
    \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{1}} & \frac{\partial Y_{ns}^{t=1}}{\partial p_{2}} & \vdots  & \frac{\partial Y_{ns}^{t=1}}{\partial p_{np+1}}\\
    \hline
    \end{array}
    \end{equation*}

And

.. math::
    

    \begin{gather*}
    \begin{array}{|c|c|}
    & Vova \\
    Viva & \begin{array}{|c c c c|c c c c|c c c c|}
    \hline
    \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}}\\
    \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}}\\
    \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots \\
    \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}}\\
    \hline
    \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}}\\
    \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}}\\
    \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots \\
    \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}}\\
    \hline
    \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{1}}{\partial p_{np+1}}\\
    \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{2}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{1}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{2}}{\partial p_{np+1}}\\
    \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots  & \cdots  & \cdots  & \vdots  & \cdots \\
    \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{1}} & \frac{\partial Y^{t=1}_{ns}}{\partial p_{2}} & \vdots  & \frac{\partial Y^{t=1}_{ns}}{\partial p_{np+1}}\\
    \hline
    \end{array} \\
    \end{array}
    \end{gather*}
    



Contents
========

.. toctree::
   :maxdepth: 2

   License <license>
   Authors <authors>
   Changelog <changelog>
   Module Reference <api/modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _toctree: http://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html
.. _reStructuredText: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _references: http://www.sphinx-doc.org/en/stable/markup/inline.html
.. _Python domain syntax: http://sphinx-doc.org/domains.html#the-python-domain
.. _Sphinx: http://www.sphinx-doc.org/
.. _Python: http://docs.python.org/
.. _Numpy: http://docs.scipy.org/doc/numpy
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _matplotlib: https://matplotlib.org/contents.html#
.. _Pandas: http://pandas.pydata.org/pandas-docs/stable
.. _Scikit-Learn: http://scikit-learn.org/stable
.. _autodoc: http://www.sphinx-doc.org/en/stable/ext/autodoc.html
.. _Google style: https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings
.. _NumPy style: https://numpydoc.readthedocs.io/en/latest/format.html
.. _classical style: http://www.sphinx-doc.org/en/stable/domains.html#info-field-lists
