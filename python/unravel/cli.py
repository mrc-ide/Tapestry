import click
# from unravel.coi.commands import coi
from unravel.sample.commands import sample
# from unravel.population.commads import population

@click.group()
def cli():
    """
    Plot output files from Tapestry
    
    """
    
    pass

# cli.add_command(coi)
cli.add_command(sample)
# cli.add_command(population)


if __name__ == "__main__":
    cli()