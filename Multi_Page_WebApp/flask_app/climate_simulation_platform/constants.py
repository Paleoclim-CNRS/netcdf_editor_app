from collections import OrderedDict
import numpy

tasks = {
    "python": ["regrid", "routing", "heatflow", "ahmcoef", "pft", "invalidate"],
    "mosaicx": ["calculate_weights"],
    "panel": ["internal_oceans", "passage_problems", "subbasins"],
}

invalidates = OrderedDict(
    {
        "regrid": ["internal_oceans", "routing"],
        "internal_oceans": ["routing"],
        "routing": [
            "pft",
            "passage_problems",
            "subbasins",
            "heatflow",
            "ahmcoef",
            "calculate_weights",
        ],
        "passage_problems": ["subbasins", "heatflow", "ahmcoef", "calculate_weights"],
        "subbasins": ["calculate_weights"],
    }
)

file_dependancies = {
    "bathy": ["raw"],
    "pft": ["bathy"],
    "routing": ["raw"],
    "soils": ["raw"],
    "topo_high_res": ["raw"],
    "heatflow": ["bathy"],
    "ahmcoef": ["bathy"],
}

no_params = ["heatflow", "ahmcoef"]


def order_steps(steps):
    # Order the steps correctly
    blocking_steps = []
    dependant_steps = []
    for i in range(len(steps)):
        # Test to see if it is a blocking step -> it has steps that depend on it
        # Because invalidates is an OrderedDict we lookup the location in which it is seen
        # TODO We should probably use some sort of tree to run through this automatically and not depend on the order
        if steps[i] in invalidates.keys():
            blocking_steps.append([i, list(invalidates).index(steps[i])])
        else:
            dependant_steps.append(
                [
                    i,
                    len(invalidates.keys())
                    - numpy.argmax(
                        [steps[i] in val for val in invalidates.values()][::-1]
                    ),
                ]
            )
    dependant_steps = numpy.array(dependant_steps)
    blocking_steps = numpy.array(blocking_steps)
    out = []
    if len(blocking_steps):
        out = numpy.array(steps)[
            blocking_steps[numpy.argsort(blocking_steps[:, 1])][:, 0]
        ]
    if len(dependant_steps):
        out = list(
            numpy.insert(
                out, dependant_steps[:, 1], numpy.array(steps)[dependant_steps[:, 0]]
            )
        )
    return out


def invalidated(root):
    return walk_through(invalidates, root, [])


def dependant_files(root):
    results = walk_through(file_dependancies, root, [])
    results.append(root)
    return results


def walk_through(dictionnary, root, _list=[]):
    if root not in dictionnary.keys():
        return []
    for child in dictionnary[root]:
        _list.append(child)
        walk_through(dictionnary, child, _list)
    return list(set(_list))
