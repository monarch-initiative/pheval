import click
from pheval_benchmark.create_batch_commands import prepare_exomiser_batch
from pheval_benchmark.assess_prioritisation import assess_prioritisation


@click.group()
def pheval_benchmark():
    """ PhEval - A benchmarking CLI. """


pheval_benchmark.add_command(prepare_exomiser_batch)
pheval_benchmark.add_command(assess_prioritisation)


if __name__ == '__main__':
    pheval_benchmark()
