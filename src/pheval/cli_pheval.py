"""
Monarch Initiative
"""

import click


@click.command()
@click.option("--input", "-i", metavar="INPUT", required=True, help="Input file")
def run(input):
    print("running pheval::run command")
