import click
from pheval_benchmark.create_batch_commands import prepare_exomiser_batch
from pheval_benchmark.assess_prioritisation import assess_prioritisation
from pheval_benchmark.create_spiked_vcf import spike_vcf
from pheval_benchmark.create_noisy_phenopackets import create_noisy_phenopackets
import warnings
warnings.filterwarnings("ignore")


@click.group()
def pheval_benchmark():
    """ PhEval - A benchmarking CLI. """

pheval_benchmark.add_command(spike_vcf)
pheval_benchmark.add_command(create_noisy_phenopackets)
pheval_benchmark.add_command(prepare_exomiser_batch)
pheval_benchmark.add_command(assess_prioritisation)


if __name__ == '__main__':
    pheval_benchmark()
