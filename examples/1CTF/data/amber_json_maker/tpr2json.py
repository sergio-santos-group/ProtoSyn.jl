"""
    tpr2json.py mol.txt > amber_top.json

Read a mol.txt file containing the dumped (using gmx dump) TPR information and parse into a amber topology JSON file.
"""

import re

class bondtypes:
    matcher         = re.compile('^\s+functype\[\d+\]=BOND')
    content_parser  = re.compile('(b0A|cbA)=([ \d\.\+\-e]+)')
    title_parser    = re.compile('^\s+functype\[(\d+)\]=BOND')

    @staticmethod
    def parse(line, *args):
        title       = int(bondtypes.title_parser.findall(line)[0])
        pre_content = dict(bondtypes.content_parser.findall(line))
        content     = dict((k, float(v)) for k, v in pre_content.iteritems())
        return title, content


class angletypes:
    matcher         = re.compile('^\s+functype\[\d+\]=ANGLES')
    content_parser  = re.compile('(thA|ctA)=([ \d\.\+\-e]+)')
    title_parser    = re.compile('^\s+functype\[(\d+)\]=ANGLES')

    @staticmethod
    def parse(line, *args):
        title       = int(angletypes.title_parser.findall(line)[0])
        pre_content = dict(angletypes.content_parser.findall(line))
        content     = dict((k, float(v)) for k, v in pre_content.iteritems())
        return title, content


class dihedraltypes:
    matcher         = re.compile('^\s+functype\[\d+\]=(PDIHS|PIDIHS)')
    content_parser  = re.compile('(phiA|cpA|mult)=([ \d\.\+\-e]+)')
    title_parser    = re.compile('^\s+functype\[(\d+)\]=(?:PDIHS|PIDIHS)')

    @staticmethod
    def parse(line, *args):
        title       = int(dihedraltypes.title_parser.findall(line)[0])
        pre_content = dict(dihedraltypes.content_parser.findall(line))
        content     = dict((k, float(v)) for k, v in pre_content.iteritems())
        return title, content


class nonbondedtypes:
    matcher         = re.compile('^\s+functype\[\d+\]=LJ_SR')
    content_parser  = re.compile('(c6|c12)=([ \d\.\+\-e]+)')
    title_parser    = re.compile('^\s+functype\[(\d+)\]=LJ_SR')

    @staticmethod
    def parse(line, *args):
        title       = int(nonbondedtypes.title_parser.findall(line)[0])
        pre_content = dict(nonbondedtypes.content_parser.findall(line))
        content     = dict((k, float(v)) for k, v in pre_content.iteritems())
        return title, content


class pairs:
    matcher         = re.compile('^\s+\d+ type=\d+ \(LJ14\)')
    content_parser  = re.compile('\d+')

    @staticmethod
    def parse(line, *args):
        content     = [int(k) for k in pairs.content_parser.findall(line)]
        #map(int, pairs.content_parser.findall(line))
        return content[0], {"a1": content[3], "a2": content[4]}


class proper_dihedrals:
    matcher         = re.compile('^\s+\d+ type=\d+ \(PDIHS')
    content_parser  = re.compile('\d+')

    @staticmethod
    def parse(line, *args):
        content     = [int(k) for k in proper_dihedrals.content_parser.findall(line)]
        return content[0], {"type": content[1], "a1": content[2],
            "a2": content[3], "a3": content[4], "a4": content[5]}


class improper_dihedrals:
    matcher         = re.compile('^\s+\d+ type=\d+ \(PIDIHS')
    content_parser  = re.compile('\d+')

    @staticmethod
    def parse(line, *args):
        content     = [int(k) for k in improper_dihedrals.content_parser.findall(line)]
        return content[0], {"type": content[1], "a1": content[2],
            "a2": content[3], "a3": content[4], "a4": content[5]}


class angles:
    matcher         = re.compile('^\s+\d+ type=\d+ \(ANGLES\)')
    content_parser  = re.compile('\d+')

    @staticmethod
    def parse(line, *args):
        content     = [int(k) for k in angles.content_parser.findall(line)]
        return content[0], {"type": content[1], "a1": content[2], "a2": content[3], "a3": content[4]}


class bonds:
    matcher         = re.compile('^\s+\d+ type=\d+ \(BONDS\)')
    content_parser  = re.compile('\d+')

    @staticmethod
    def parse(line, *args):
        content     = [int(k) for k in bonds.content_parser.findall(line)]
        return content[0], {"type": content[1], "a1": content[2], "a2": content[3]}


class atoms:
    matcher         = re.compile('^\s+atom\[\s+\d+]=')
    content_parser  = re.compile('(type|q)=( *[\d\.\+\-e]+)')
    title_parser    = re.compile('^\s+atom\[\s+(\d+)]=')

    @staticmethod
    def parse(line, *args):
        title       = int(atoms.title_parser.findall(line)[0])
        pre_content = dict(atoms.content_parser.findall(line))
        content     = dict((k, float(v)) for k, v in pre_content.iteritems())
        return title, content


class atomtypes:
    matcher = re.compile('^\s+atomtype\[\s+\d+]=')
    title_parser = re.compile('^\s+atomtype\[\s+(\d+)]=')

    @staticmethod
    def parse(line, *args):
        return int(atomtypes.title_parser.findall(line)[0]), {}

    
class atomtypenames:
    matcher        = re.compile('^\s+type\[\d+\]=')
    content_parser = re.compile('name="(.*)",')
    title_parser   = re.compile('^\s+type\[(\d+)\]=')

    @staticmethod
    def parse(line, *args):
        title       = int(atomtypenames.title_parser.findall(line)[0])
        content     = atomtypenames.content_parser.findall(line)[0]
        return title, content


class excls:
    matcher = re.compile('^\s+excls\[\d+\]')
    parser = re.compile('\d+')

    @staticmethod
    def parse(line, iterator, ffield):
        items = excls.parser.findall(line)
        while '}' not in line:
            line = next(iterator)
            items += excls.parser.findall(line)
        return items[0], [int(x) for x in items[3:]]


class Parser(object):
    fields = (
        bondtypes,
        angletypes,
        dihedraltypes,
        nonbondedtypes,
        atoms,
        excls,
        atomtypes,
        pairs,
        bonds,
        angles,
        proper_dihedrals,
        improper_dihedrals,
        atomtypenames
    )
    
    @staticmethod
    def parse(fhandle):
        iterator = iter(fhandle)
        ffield = {}

        while True:
            
            try:
                line = next(iterator)
            except StopIteration:
                break
            
            for field in Parser.fields:
                if field.matcher.match(line):
                    if field.__name__ not in ffield:
                        ffield[field.__name__] = {}
                    
                    (title, content) = field.parse(line, iterator, ffield)
                    ffield[field.__name__][title] = content
                    break
        
        return ffield


def clean(pre_ffield):
    #Convert dictionaries to arrays
    for key in ["bonds", "angles", "proper_dihedrals", "improper_dihedrals", "pairs", "atoms"]:    
        pre_ffield[key] = pre_ffield[key].values()

    #Remove 'A' characters from dictionary entries (ugly but it works)
    for key in ["bondtypes", "angletypes", "dihedraltypes"]:
        for key1, value1 in pre_ffield[key].items():
            for key2, value2 in value1.items():
                if key2[-1] == 'A':
                    nkey = key2.strip('A')
                    pre_ffield[key][key1][nkey] = value2
                    pre_ffield[key][key1].pop(key2)

    #Filter atomtypes to keep only the necessary combinations
    atomtype_count = len(pre_ffield["atomtypes"])
    nonbonded = {}
    for i in xrange(0, atomtype_count):
        idx = i * atomtype_count + 1 * i
        pre_ffield["atomtypes"][i] = pre_ffield["nonbondedtypes"][idx]
    
    #Cast atom types to int and add names to atomtypes
    for idx, atom in enumerate(pre_ffield["atoms"]):
        atom["type"] = int(atom["type"])
        pre_ffield["atomtypes"][atom["type"]]["name"] = pre_ffield["atomtypenames"][idx]

    #convert C6 and C12 to Sigma and Epsilon
    for atomtype in pre_ffield["atomtypes"].values():
        # sigma = (C12 / C6) ^ (1 / 6)
        # epsilon = (C6 ^ 2) / (4 * C12)
        if atomtype["c6"] == 0.0:
            atomtype["sigma"] = 0.0
        else:
            atomtype["sigma"] = pow((atomtype["c12"] / atomtype["c6"]), (1.0/6.0))
        if atomtype["c12"] == 0.0:
            atomtype["epsilon"] = 0.0
        else:
            atomtype["epsilon"] = pow(atomtype["c6"], 2) / (4.0 * atomtype["c12"])
        atomtype.pop("c6")
        atomtype.pop("c12")

    #Remove unnecessary entries
    pre_ffield.pop("nonbondedtypes")
    pre_ffield.pop("atomtypenames")

    return pre_ffield


if __name__ == '__main__':
    import json
    import sys

    fname = sys.argv[1]
    
    with open(fname) as fin:
        pre_ffield = Parser.parse(fin)
        ffield = clean(pre_ffield)
        print json.dumps(ffield, indent = 2)