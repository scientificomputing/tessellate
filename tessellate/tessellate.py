# -*- coding: utf-8 -*-

"""Main module."""


def createnodejson_forcrd(id, conformer, puckercoords, nextconformer, ringsize):
    try:
        roundedpuckercoords=[]
        for x_ind in puckercoords[1:-2].split(','):
            roundedpuckercoords.append("{0:.2f}".format(float(x_ind)))
    except Exception as e:
        logging.error('%s %s', e, puckercoords)
    jsonobject = {
        'pdbid': id, # note this could be any id in this case
        'conformer': conformer,
        'puckercoords': roundedpuckercoords,
        'nextconformer': nextconformer,
        'ringsize': ringsize
    }
    return jsonobject, "%s %s %s %s\n" % (conformer, puckercoords, nextconformer, ringsize)

def analyse_pucker(inputfile, logfile = 'log_of_analyse.log', outputfile = "output_of_analyse.out", verbose = False, working_folder=""):
    pass
    try:
        import os
        import sys
        import time
        import logging as logging
        import logging
        # import pucker as puc
        import tessellate.utils.pucker as puc
        import json
        import tessellate.utils.helperfunctions as helperfunctions
    except Exception as e:
        print("Error - Cannot import module %s", e)
        exit(1)

    all_pucker_json = helperfunctions.init_all_pucker_dictionary()
    PDATA = {}

    logging.basicConfig(filename=os.path.join(working_folder, logfile), level=logging.DEBUG)

    outputfile = open(os.path.join(working_folder, outputfile), 'w')
    outputfile.write(
        "PUCKER_COORDS CONFORMER ORIG_CONFORMER RING_SIZE\n")

    nodejson = []

    # redirect all verbose info to verbose
    verbose_log = []

    # - does the list of coordinates exist
    inputstream = open(inputfile, 'r')
    for line in inputstream:
        logging.debug('pucker line item %s', line)
        chunked = line.split()
        id = chunked[0]
        print(chunked)
        ringsize = (len(chunked) - 1) / 3
        logging.critical('Ringsize %s',ringsize)
        print(ringsize)
        if ringsize < 5 or ringsize > 8:  # not calculating for these sizes at the moment
            logging.error('Systems containing less than 5 or greater than 8 atoms are not supported. skipping')
            verbose_log.append("Systems containing less than 5 or greater than 8 atoms are not supported. ")
            print("exiting")
            exit(1) # FIXME!, rather ignore and calculate all others..
        allcoors1 = chunked[1:len(chunked)]
        allcoors = [float(i) for i in allcoors1]
        pobj = puc.Pucker(tuple(allcoors))
        print(str(pobj))
        if pobj and pobj.isvalid:
            try:
                thisframeTD = pobj.calculate_triangular_tessellation()
                print(thisframeTD)
                conformer = pobj.deduce_canonical_conformation()
                nextconformer = pobj.deduce_canonical_conformation(nextguess=True)
                verbose_log.append('<div> Coordinates are : ' + str(thisframeTD) + " " + str(conformer[0]) + " " + str(
                    nextconformer[0]) + " " + str(pobj.ringsize) + "<br /> </div>")
                print(str(thisframeTD))

                node, log = createnodejson_forcrd(id, conformer[0], str(thisframeTD), nextconformer[0], pobj.ringsize)
                logging.info(log)
                outputfile.write(log)
                nodejson.append(node)
                try:
                    all_pucker_json[pobj.ringsize][conformer[0]] += 1
                except:
                    all_pucker_json[pobj.ringsize][conformer[0]] = 1

                # create data for visualisation
                try:
                    PDATA[id].append(
                        {"conformer": conformer[0],
                         "TD": thisframeTD})
                except:
                    PDATA[id] = [
                        {"conformer": conformer[0],
                         "TD": thisframeTD}]

            except Exception as e:
                logging.debug('Pucker object valid, but calc, classify etc. failed  %s %s', str(pobj._coords), e)
                raise e
                print(str(pobj._coords),e)
                pass

    pucsummary = open('pucker_summary.csv', 'w')
    pucsummary.write("size,conformer,count\n")
    for key in all_pucker_json[5]:
        pucsummary.write(",".join(("FIVE", key, str(all_pucker_json[5][key]), "\n")))
        try:
            all_pucker_json["five"].append({"label": key, "value": float(all_pucker_json[5][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(5,key), "color_aster":helperfunctions.getcolor(5,key)})
        except:
            all_pucker_json["five"] = [{"label": key, "value": float(all_pucker_json[5][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(5,key), "color_aster":helperfunctions.getcolor(5,key)}]

    for key in all_pucker_json[6]:
        pucsummary.write(",".join(("SIX", key, str(all_pucker_json[6][key]), "\n")))
        try:
            all_pucker_json["six"].append({"label": key, "value": float(all_pucker_json[6][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(6,key), "color_aster":helperfunctions.getcolor(6,key)})
        except:
            all_pucker_json["six"] = [{"label": key, "value": float(all_pucker_json[6][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(6,key), "color_aster":helperfunctions.getcolor(6,key)}]
    for key in all_pucker_json[7]:
        pucsummary.write(",".join(("SEVEN", key, str(all_pucker_json[7][key]), "\n")))
        try:
            #all_pucker_json["seven"].append({"label": key, "value": str(all_pucker_json[7][key])})
            all_pucker_json["seven"].append({"label": key, "value": float(all_pucker_json[7][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(7,key), "color_aster":helperfunctions.getcolor(7,key)})
        except:
            #all_pucker_json["seven"] = [{"label": key, "value": str(all_pucker_json[7][key])}]
            all_pucker_json["seven"] = [{"label": key, "value": float(all_pucker_json[7][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(7,key), "color_aster":helperfunctions.getcolor(7,key)}]
    for key in all_pucker_json[8]:
        pucsummary.write(",".join(("EIGHT", key, str(all_pucker_json[8][key]), "\n")))
        try:
            #all_pucker_json["eight"].append({"label": key, "value": str(all_pucker_json[8][key])})
            all_pucker_json["eight"].append({"label": key, "value": float(all_pucker_json[8][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(8,key), "color_aster":helperfunctions.getcolor(8,key)})
        except:
            #all_pucker_json["eight"] = [{"label": key, "value": str(all_pucker_json[8][key])}]
            all_pucker_json["eight"] = [{"label": key, "value": float(all_pucker_json[8][key]), "id": key, "weight": 1,"order": helperfunctions.getorder(8,key), "color_aster":helperfunctions.getcolor(8,key)}]
    pucsummary.close()
    # <>  dump all the pucker stats as json files
    d2 = [key for key in all_pucker_json["five"]]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'five.json'))
    d2 = [key for key in all_pucker_json["six"]]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'six.json'))
    d2 = [key for key in all_pucker_json["seven"]]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'seven.json'))
    d2 = [key for key in all_pucker_json["eight"]]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'eight.json'))

    # <> dump pucker stats out to be used in discrete bar charts, if value is not a float then graph doesn't work properly "1" -> 1.0
    d2 = [{'key': 'five cycles', "color": "#d67777", 'values': [key for key in all_pucker_json["five"]]}]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'five_bar.json'))
    d2 = [{'key': 'six cycles', "color": "#d67777", 'values': [key for key in all_pucker_json["six"]]}]
    #d2 = [{ 'key' : 'six cycles', "color": "#d67777", 'values' : [sorted([(key,value) for (key,value) in all_pucker_json["six"]])]} ]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'six_bar.json'))
    d2 = [{'key': 'seven cycles', "color": "#d67777", 'values': [key for key in all_pucker_json["seven"]]}]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'seven_bar.json'))
    d2 = [{'key': 'eight cycles', "color": "#d67777", 'values': [key for key in all_pucker_json["eight"]]}]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'eight_bar.json'))


    d = {"name": "PDB", "children": [{'name': key, "children": value} for key, value in PDATA.items()]}
    helperfunctions.write_to_json(d, os.path.join(working_folder, 'pucker_flare.json'))

    d2 = {"five": [key for key in sorted(all_pucker_json["five"], key=lambda dist: dist["order"] )],"six": [key for key in sorted(all_pucker_json["six"], key=lambda dist: dist["order"] )],"seven": [key for key in sorted(all_pucker_json["seven"], key=lambda dist: dist["order"] )],"eight": [key for key in sorted(all_pucker_json["eight"], key=lambda dist: dist["order"] )]}
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'aster.json'))

    d2 = [key for key in nodejson]
    helperfunctions.write_to_json(d2, os.path.join(working_folder, 'pucker_table.json'))
