"""Placeholders for OCP names pykentron does not provide."""


class Placeholder:
    """A name that imports but cannot be used."""

    _name = "?"

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(f"OCP.{self._name} is not provided by pykentron")

    def __class_getitem__(cls, item):
        return cls


def placeholders(module_name: str, names: str) -> dict:
    out = {}
    for name in names.split():
        out[name] = type(name, (Placeholder,), {"_name": f"{module_name}.{name}"})
    return out


def fallback(module_name: str):
    """A module `__getattr__` making every other name of the module a placeholder."""

    def __getattr__(name):
        if name.startswith("__"):
            raise AttributeError(name)
        return type(name, (Placeholder,), {"_name": f"{module_name}.{name}"})

    return __getattr__
