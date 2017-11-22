# -*- coding: utf-8 -*-

"""Main module."""

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def createnodejson_forcrd(id, conformer, puckercoords, nextconformer, ringsize):
    roundedpuckercoords = []
    try:
        for x_ind in puckercoords[1:-2].split(','):
            roundedpuckercoords.append("{0:.2f}".format(float(x_ind)))
    except Exception as e:
        logger.error('%s %s', e, puckercoords)
    jsonobject = {
        'pdbid': id,  # note this could be any id in this case
        'conformer': conformer,
        'puckercoords': roundedpuckercoords,
        'nextconformer': nextconformer,
        'ringsize': ringsize
    }
    return jsonobject, "%s %s %s %s\n" % (conformer, puckercoords, nextconformer, ringsize)


def analyse_pucker(inputfile, outputfile="output_of_analyse", output_format="json",
                   working_folder=""):
    pass
    try:
        import os
        import sys
        import time
        import tessellate.utils.pucker as puc
        import json
        import tessellate.utils.helperfunctions as helperfunctions
    except Exception as e:
        print("Error - Cannot import module %s", e)
        exit(1)

    outputfile = ".".join([outputfile, output_format])
    txt = False
    if output_format == "txt":
        txt = True
        outputfile = open(os.path.join(working_folder, outputfile), 'w')
        outputfile.write("tessellate 0.1 txt\n")
        outputfile.write(
            "PUCKER_COORDS CONFORMER ORIG_CONFORMER RING_SIZE\n")
        logger.critical("to be improved, does not match json output at present")
    nodejson = []

    # - does the list of coordinates exist
    inputstream = open(inputfile, 'r')
    for line in inputstream:
        logger.debug('pucker line item %s', line)
        chunked = line.split()
        id = chunked[0]
        ringsize = int((len(chunked) - 1) / 3)
        logger.debug('Ringsize %s', ringsize)
        if ringsize < 5 or ringsize > 8:  # not calculating for these sizes at the moment
            logger.error('Systems containing less than 5 or greater than 8 atoms are not supported. skipping')
            exit(1)  # FIXME!, rather ignore and calculate all others..
        allcoors1 = chunked[1:len(chunked)]
        allcoors = [float(i) for i in allcoors1]
        pobj = puc.Pucker(tuple(allcoors))

        if pobj and pobj.isvalid:
            try:
                thisframeTD = pobj.calculate_triangular_tessellation()
                logger.debug(thisframeTD)
                conformer = pobj.deduce_canonical_conformation()
                nextconformer = pobj.deduce_canonical_conformation(nextguess=True)
                logger.debug('Coordinates are %s %s %s %s : ', str(thisframeTD), str(conformer[0]),
                             str(nextconformer[0]), str(pobj.ringsize))
                node, log = createnodejson_forcrd(id, conformer[0], str(thisframeTD), nextconformer[0], pobj.ringsize)
                logger.info(log)
                if txt:
                    outputfile.write(log)
                nodejson.append(node)

            except Exception as e:
                logger.error('Pucker object valid, but calc, classify etc. failed  %s %s', str(pobj._coords), e)
                raise e
                pass

    if output_format == "json":
        d2 = [key for key in nodejson]
        helperfunctions.write_to_json(d2, os.path.join(working_folder, outputfile))
    elif output_format == "bson":
        d2 = [key for key in nodejson]
        helperfunctions.write_to_bson(d2, os.path.join(working_folder, outputfile))
