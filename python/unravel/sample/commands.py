import click


@click.command(short_help="Create plots for an individual sample.")
@click.option(
    "-i",
    "--input_dir",
    type=click.Path(exists=True),
    required=True,
    help="Directory containing Tapestry outputs, for an individual sample."
)
def sample(input_dir):
    """
    Plot Tapestry outputs for an individual sample

    """
    from .main import sample
    sample(input_dir)

