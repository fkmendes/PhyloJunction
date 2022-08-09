# Notes for developers

## Making packages and modules discoverable on VSCode

The first thing you must do before you can start coding within the `PhyloJunction` framework
is to set it up within an integrated development environment (IDE).
The IDE of choice here is VSCode, a very popular, clean and fast IDE.

The first step is to go under "File" > "Open Folder..." and choose `PhyloJunction/`, the git repository root folder.
If you click around for a bit, you will notice that the module files (`.py` files) are inside
the `src/` directory, nested within different packages (e.g., `phylojunction/`, `phylojunction/distribution/`).
We must point our python interpreter to the topmost package, which is called `phylojunction/` so that all modules can be found.
Here, we must create a `launch.json` file, which VSCode can help us with.

Click the button on the left that has a bug with a triangle ("Run and Debug"), and then "Create launch.json file", choosing module.
This should create a that looks like this:

```
{
    // Use IntelliSense to learn about possible attributes.
    // Hover to view descriptions of existing attributes.
    // For more information, visit: https://go.microsoft.com/fwlink/?linkid=830387
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Python: Current File",
            "type": "python",
            "request": "launch",
            "program": "${file}",
            "console": "integratedTerminal",
            "justMyCode": true
        }
    ]
}
```

In order to tell the interpreter to look inside `src/phylojunction/` for modules, add the following
line below `"justMyCode"`:

```
"env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
```

After doing this, you might notice if you choose one of the `.py` files that the import statements at the top of the file will be marked as problematic.
This is because VSCode itself cannot discover the imported packages and modules.
We must fix that, which we do by creating a `settings.json` file inside hidden folder `.vscode`.
Enter shortcut "CMD + Shift + P", type `settings.json` and choose "Open Workspace Settings (JSON)".

Now inside `settings.json`, add this:

```
{
    "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
}
```

## Configuring unit tests (using unittest) on VSCode

If you have an Apple machine, the shortcut on VSCode is "CMD + Shift + P", "Python: Configure Tests".
This will add to your `settings.json` file inside hidden folder `.vscode`.
Choose `.` as the directory containing tests, and `test_*` as the prefix for tests.

Your `settings.json` will look like this:

```
{
    "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ],
    "python.testing.unittestArgs": [
        "-v",
        "-s",
        ".",
        "-p",
        "test_*.py"
    ],
    "python.testing.pytestEnabled": false,
    "python.testing.unittestEnabled": true
}
```

If VSCode refuses to find the tests, one thing that might be set up wrong is the Python interpreter VSCode is using.
Unit tests usually have many dependencies -- if those dependencies are not found because the wrong interpreter is used (the native one that comes with the machine, say, instead of an interpreter installed with Homebrew), then the tests will not be found.

## Type checking

`mypy`'s documentation suggests calling mypy on entire packages, e.g.,

```
cd src/phylojunction/
mypy data/
mypy calculation/
```

But it's also possible to call it on individual files:

```
mypy data/tree.py
mypy utility/helper_functions.py --namespace-packages --explicit-package-bases
mypy distribution/dn_discrete_sse.py --namespace-packages --explicit-package-bases
```

`mypy` will verify type hinting in source files against mypy stubs (files ending with `.pyi`).
Stubs must be kept up-to-date by developers so `mypy` can find incongruencies.

Sometimes I have noticed that calling `mypy` on individual modules that depend on other modules within the same package fails (e.g., `helper_functions.py` is in the same package as, and depends on `exception_classes.py`).
Fixing it can be done by guiding `mypy`'s import discovery; 
This is what the flags `--namespace-packages` together with `--explicit-package-bases` do  above.
They are necessary even when an `__init__.py` is present in the same folder as the modules being type checked.

PJ has a script that calls `mypy` automatically on everything:

```
cd PhyloJunction/
python3 typecheck.py
```

This script should be called and be successful every time new code is written, before committing and pushing.

## Editable pip install for developers

Building and installing PhyloJunction is done using `setuptools`, which is specified in `pyproject.toml`.
Instructions for `setuptools` are provided in `setup.cfg`.

Start by cloning the PhyloJunction repository, and going to its root (i.e., the
`PhyloJunction` folder).

If you are on an Apple machine that obstrusively curtails file system access, run the following command:

```
python3 -m pip install --user -e .
```

The `--user` flag is the way around said restrictions.
This will tell `pip` to place the entry point executables (`pjcli` and `pjgui`) in  `/Users/user_name/Library/Python/3.9/bin/` if you are using Python 3.9, say.
Remember to add this path to the PATH system variable if you want to call the executables from anywhere on your file system.
Another option is to run the following command:

```
python3 -m pip install --prefix ~/.local/bin -e .
```

Which creates directories `bin/` (placing the executables therein) and `lib/python3.9/site-packages/` (placing the egg-link therein) inside `~/.local`.
This is an attractive option if you normally already have `~/.local/bin` as part of your PATH variable.

## Documentation

PhyloJunction is documented automatically with `sphinx`, so make sure you have this program installed:

```
pip3 install sphinx
```

`sphinx` is a Python program that converts docstrings into .rst files, and then from those and other .rst files it creates it can produce .html files for online documentation.

The first time documentation was done, a `docs/` directory was manually created in the git repo's root, `PhyloJunction/`.
Then:

```
cd docs/
sphinx-quickstart
```

The answer was "yes" for separate source and build folders and "PhyloJunction" for project name.

Then we needed to activate the Napoleon extension of sphinx, so that we can make use of the Google-style documentation strings PhyloJunction uses.
Inside `source/conf.py` we added:

```
extensions = [ 'sphinx.ext.napoleon' ]
```

With the following options (mostly the defaults):

```
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False # not using admonitions                                                                                                                          
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True
```

Because PhyloJunction's code is not in `docs/`, but in fact in `src/phylojunction/`, we must tell sphinx via `conf.py`.
At the top of `conf.py`, add:

```
import os
import sys
sys.path.insert(0, os.path.abspath("../src/phylojunction/"))
```

### Theme and other configurations

We will use the "Read the Docs" (RTD) theme to document PJ, as it is straightforward and familiar.
First, install it:

```
pip3 install sphinx-rtd-theme
```

And the following must be added to `conf.py`:

```
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
html_theme = "sphinx_rtd_theme"
```

We will also add a few .rst directives to `conf.py`.
These work like "macros" that will be available for all .html files rendered by sphinx, and help us replace strings and add hyperlinks to them automatically.

```
rst_prolog = """
.. |fkm| replace:: Fabio K. Mendes
.. |pj| replace:: PhyloJunction"""

# related to installation, requirements and executing 
rst_prolog += """\
.. |Python| replace:: Python
.. _Python: http://www.python.org
.. |Python3| replace:: Python 3
.. |dendropy| replace:: DendroPy
.. _DendroPy: https://dendropy.org/
.. |msprime| replace:: msprime
.. _msprime: https://tskit.dev/msprime
.. |numpy| replace:: Numpy
.. _NumPy: https://numpy.org/
.. |scipy| replace:: SciPy
.. _SciPy: https://scipy.org
.. |matplotlib| replace:: Matplotlib
.. _Matplotlib: https://matplotlib.org
.. |pandas| replace:: pandas
.. _pandas: https://pandas.pydata.org
.. |pip| replace:: pip
.. _pip: https://pypi.python.org/pypi/pip
.. |pysimplegui| replace:: PySimpleGUI
.. _PySImpleGUI: https://pysimplegui.org/en/latest/"""

# other
rst_prolog += """\
.. _R: https://www.r-project.org
.. |java| replace:: Java
.. _Java: https://www.java.com/en
.. |beast| replace:: BEAST
.. _BEAST: https://beast.community
.. |beast2| replace:: BEAST 2
.. _BEAST 2: https://www.beast2.org
.. |rb| replace:: RevBayes
.. _RevBayes: https://revbayes.github.io
.. |cpp| replace:: C++
.. _C++: https://isocpp.org/"""
```

Now that the configuration files are ready, we can call the program that automatically generate .rst files from all docstrings in PJ:

```
cd docs/
sphinx-apidoc -o ./source ../src/phylojunction
```

Do not forget to run this command again if docstrings are updated.
Finally, we call the command that produces the .html files:

```
make clean && make html
```

This command will need to be called every time we want to update the .html files once the .rst files are modified.

### Documentation website structure

The .rst file that will be converted into `index.html`, the welcome page, is `source/index.rst` (the master .rst file).
This file contains the information for citing and installing PJ.

The remaining .rst files are organized within manually created folders:

* user/index.rst
* developer/index.rst

These files are included in the website via the `.. toctree::` .rst directive in the master .rst file.