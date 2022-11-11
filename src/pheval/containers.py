"""Containers module."""

from dependency_injector import containers, providers

from pheval import runner, runners


class Container(containers.DeclarativeContainer):
    config = providers.Configuration()

    selector = providers.Selector(
        config.runner.type,
        foo=providers.Singleton(runners.FooRunner),
        exomiser=providers.Singleton(runners.ExomiserRunner),
    )

    service = providers.Factory(
        runner.Runner,
        runner=selector,
    )
