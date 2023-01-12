"""implementations"""
from functools import cache

from class_resolver import ClassResolver

from pheval.runners.runner import PhEvalRunner


@cache
def get_implementation_resolver() -> ClassResolver[PhEvalRunner]:
    """get_implementation_resolver

    Returns:
        ClassResolver[PhEvalRunner]: _description_
    """
    implementation_resolver: ClassResolver[PhEvalRunner] = ClassResolver.from_subclasses(
        PhEvalRunner,
        suffix="Implementation",
    )

    # implementation_resolver.synonyms.update(
    #     {
    #         "exomiser": DefaultPhEvalRunner,
    #     }
    # )

    # Plugins which want to register an implementation should use
    # the entrypoint group "pheval.plugins". The name of the entry
    # point will be used as a possible match against the input scheme
    # prefix. The value of the entry point should be an implementation
    # class.
    #
    # See also:
    # https://packaging.python.org/en/latest/specifications/entry-points/
    # https://class-resolver.readthedocs.io/en/latest/api/class_resolver.ClassResolver.html#class_resolver.ClassResolver.register_entrypoint
    implementation_resolver.register_entrypoint("pheval.plugins")
    return implementation_resolver
