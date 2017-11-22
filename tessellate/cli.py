# -*- coding: utf-8 -*-

"""Console script for tessellate."""

import click


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# @click.option('--count', default=1, help='Number of greetings.')
# @click.option('--project',prompt='What is your project called?', help="Name of Project")
@click.argument('input-dir',
                type=click.Path(exists=True),
                nargs=-1,
                required=True
                # help = "Directory containing input files"
                )
@click.option('--input-format',
              type=click.Choice(['builtin', 'pdb', 'pdblist']),
              default="builtin",
              help="Input format: builtin, pdb, pdblist"
              )
@click.option('--output-format',
              type=click.Choice(['txt', 'json', 'bson']),
              default="txt",
              help="Output format: txt, json, bson"
              )
@click.option('--ignore/--no-ignore', default=True,
              help="Choose to ignore amino acids/water in PDB data when looking for cycles")
@click.option('--ignore-append', 'ignore_behaviour', flag_value="append", default=True,
              help="Append residues in ignore file to existing ignore dictionary")
@click.option('--ignore-overwrite', 'ignore_behaviour', flag_value="overwrite",
              help="Use residues in ignore file and exclude existing ignore dictionary")
@click.option('--ignore-file', default="ignorenames", help="Use this ignore file")
@click.option('--ligand/--no-ligand', default=True, help="Use builtin ligand dictionary when choosing cycle ordering")
@click.option('--ligand-append', 'ligand_behaviour', flag_value="append", default=True,
              help="Append ligands in ligand file to existing dictionary")
@click.option('--ligand-overwrite', 'ligand_behaviour', flag_value="overwrite",
              help="Use ligand file and exclude existing dictionary")
@click.option('--ligand-file', default="ligandnames", help="Use this ligand file")
@click.option('-v', '--verbose',
              count=True,
              default=0,
              help="Increase output verbosity."
              )
@click.option('-q', '--quiet',
              is_flag=True,
              help="Only show log warnings"
              )
# @click.option('-c', '--config', 'config_file',
#                     type = click.Path(exists=True, readable=True),
#                     multiple=True,
#                     help = "Specific config file to load, after those in MultiQC dir / home dir / working dir."
# )

def main(input_dir, input_format, output_format, ignore, ignore_behaviour, ignore_file, ligand, ligand_behaviour,
         ligand_file, verbose, quiet):
    """Console script for tessellate."""
    # click.echo('Logging %s !' % project)
    click.echo('Logging %s !' % input_dir)
    click.echo('Logging %s !' % input_format)
    click.echo('Logging %s !' % output_format)
    click.echo('Logging %s !' % ignore)
    click.echo('Logging %s !' % ignore_behaviour)
    click.echo('Logging %s !' % ignore_file)
    click.echo('Logging %s !' % ligand)
    click.echo('Logging %s !' % ligand_behaviour)
    click.echo('Logging %s !' % ligand_file)
    click.echo('Logging verbose %s !' % verbose)
    click.echo('Logging quiet %s !' % quiet)
    import tessellate.tessellate as tt
    import logging
    logger = logging.getLogger(__name__)

    # . making logging changes here only affects this function - all other modules read the init.
    if quiet:
        logger.critical("Changing to ERROR level, doesn't propagate")
        logger.setLevel(level=logging.ERROR)  # . only show errors and above

    if verbose > 0:  # https://docs.python.org/3/library/logging.html#levels
        logger.critical("Changing to DEBUG level, doesn't propagate")
        logger.setLevel(level=logging.DEBUG)  # show debug and worse

    logger.critical('Critical Logging on  %s !' % input_dir)
    logger.error('Error Logging on  %s !' % input_dir)
    logger.warning('Warning Logging on %s' % input_dir)
    logger.info('Info Logging on  %s !' % input_dir)
    logger.debug('Debug Logging on  %s !' % input_dir)

    # . Check input dir exists
    if input_dir == "" or input_dir == None:
        logger.critical("No input_dir defined")
        exit(1)
    # . If builtin format call the appropriate function
    if input_format == "builtin":
        logger.critical("Treating input_dir as a file: %s", input_dir)
        tt.analyse_pucker(input_dir[0], output_format=output_format)
    else:
        logger.critical("Code not yet developed for this function")
        exit(1)


if __name__ == "__main__":
    main()
