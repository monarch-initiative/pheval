# Dependencies taken from class_resolver:
# https://github.com/cthoyt/class-resolver/blob/5e51b41d08381e1754f3445c5ce4fcea4f6612d6/src/class_resolver/base.py#L276
from importlib.metadata import entry_points

# Dependency from importlib.metadata:
# https://github.com/cthoyt/class-resolver/blob/5e51b41d08381e1754f3445c5ce4fcea4f6612d6/src/class_resolver/utils.py#L45
from typing import TypeVar

from class_resolver import ClassResolver

from pheval.utils.logger import get_logger

X = TypeVar("X")
logger = get_logger()


# Create a custom ClassResolver class to modify the _from_entrypoint method.
class PhevalClassResolver(ClassResolver):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    # Modified _from_entrypoint method to raise error instead of warning
    @staticmethod
    def _from_entrypoint_custom(group: str) -> set[X]:
        elements: set[X] = set()
        for entry in entry_points(group=group):
            try:
                element = entry.load()
                logger.info(f"Loaded {entry.name} correctly")
            except (ImportError, AttributeError):
                logger.warn(f"could not load {entry.name}. See error message below.")
                raise
            else:
                elements.add(element)
        return elements

    def register_entrypoint(self, group: str) -> None:
        """Register additional entries from an entrypoint."""
        for element in self._from_entrypoint_custom(group).difference(self.lookup_dict.values()):
            self.register(element)
