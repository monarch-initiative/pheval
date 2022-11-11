"""
Monarch Initiative
"""

import click
from dependency_injector.wiring import Provide, inject

from .containers import Container
from .runner import Runner


@click.command()
@click.option("--input", "-i", metavar="INPUT", required=True, help="Input file")
@inject
def run(input, runner: Runner = Provide[Container.service]):
    print("running pheval::run command")
    runner.run()
