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

## 

## Type checking

`mypy`'s documentation suggests calling mypy on entire packages, e.g.,

```
mypy data/
mypy calculation/
```

But when doing the above, type checking succeeds when otherwise (if called on individual files separately) it would not.
It is unclear why this is the case.
So it is recommended to call `mypy` as follows:

```
mypy data/tree.py
mypy calculation/discrete_sse.py
mypy calculation/math_utils.py
mypy utility/exception_classes.py
mypy utility/helper_functions.py --namespace-packages --explicit-package-bases
mypy distribution/dn_parametric.py
mypy distribution/dn_discrete_sse.py --namespace-packages --explicit-package-bases
mypy interface/cmd_parse_utils.py
mypy interface/grammar/dn_grammar.py
mypy interface/grammar/det_fn_discrete_sse.py

TODO:
mypy interface/grammar/make_dn_discrete_sse.py
mypy interface/grammar/dn_grammar.py
```

Note that calling `mypy` on modules that depend on other modules within the same package (e.g., `helper_functions.py` is in the same package as, and depends on `exception_classes.py`) requires guiding `mypy`'s import discovery.
This is what the flags `--namespace-packages` together with `--explicit-package-bases` do (they are necessary even when an `__init__.py` is present in the same folder as the modules being type checked).