import warnings

# Dependencies taken from class_resolver: https://github.com/cthoyt/class-resolver/blob/5e51b41d08381e1754f3445c5ce4fcea4f6612d6/src/class_resolver/base.py#L276
from importlib.metadata import entry_points

# Dependency from importlib.metadata: https://github.com/cthoyt/class-resolver/blob/5e51b41d08381e1754f3445c5ce4fcea4f6612d6/src/class_resolver/utils.py#L45
from typing import TypeVar

from class_resolver import ClassResolver

X = TypeVar("X")


# Create a custom ClassResolver class to modify the _from_entrypoint method.
class PhevalClassResolver(ClassResolver):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    # Overload the _from_entrypoint method.
    @staticmethod
    def _from_entrypoint_custom(group: str) -> set[X]:
        elements: set[X] = set()
        for entry in entry_points(group=group):
            try:
                element = entry.load()
                # TODO: maybe display this with a logger?
                print(f"Loaded {entry.name} correctly")
            except (ImportError, AttributeError) as e:
                warnings.warn(f"could not load {entry.name}. See error message below.", Warning)
                raise
            else:
                elements.add(element)
        return elements

    def register_entrypoint(self, group: str) -> None:
        """Register additional entries from an entrypoint."""
        for element in self._from_entrypoint_custom(group).difference(self.lookup_dict.values()):
            self.register(element)
