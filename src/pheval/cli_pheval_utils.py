"""
"""

import click

@click.command()
@click.option("--input", "-i",
              required=True, metavar='FILE',
              help="Path to the semantic similarity profile to be scrambled.")
def scramble_semsim():
    print("running pheval_utils::scramble_semsim command")

@click.command()
@click.option("--input", "-i",
              required=True, metavar='FILE',
              help="Path to the phenopacket to be spiked.")
def scramble_phenopacket():
    print("running pheval_utils::scramble_phenopacket command")
