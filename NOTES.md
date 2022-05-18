# Notes for developers

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