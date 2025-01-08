import click


@click.command(short_help="Create plots for a population of samples.")
@click.option(
    "-i",
    "--input_dir",
    type=click.Path(exists=True),
    required=True,
    help="Directory containing tapestry output folders."
)
def population(input_dir):
    """
    Plot Tapestry outputs for an individual sample

    """
    from .main import population
    population(input_dir)