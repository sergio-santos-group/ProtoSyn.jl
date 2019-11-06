import sys
import re


regex = re.compile(r"^\s*\w+\s+\w+\s+")
block = re.compile(r"^\s*\[\s*(?P<comp>\w+)types\s*\]")

fname = sys.argv[1]

ntypes = 0

components = []

class Component(object):
    def __init__(self, k, line):
        items = line.strip().split(";", 1)
        fields = items[0].split()
        self.key = tuple(fields[:k])

        i = 1 if k>1 else 0
        self.values = [[f] for f in fields[k+i:]]
        self.comments = items[1:]
    
    def __iadd__(self, other):
        if self.key != other.key:
            raise Exception("keys must be equal")
        
        self.values = [i+j for i,j in zip(self.values, other.values)]
        self.comments += other.comments
        return self

    def __str__(self):
        types = ["%4s"%repr(k).replace("\'","\"") for k in self.key]
        if len(self) == 1:
            values = ['%10s'%v[0] for v in self.values]
        else:
            values = ['[' + ', '.join(v) + ']' for v in self.values]

        s = ', '.join(types + values)
        return f'    ({s}),  # {";; ".join(self.comments)}'
    
    def __len__(self):
        return len(self.values[0])

with open(fname) as fin:
    for line in fin:
        m = block.match(line) 
        if m:
            
            if m["comp"] == "atom":
                ntypes=1
            elif m["comp"] == "bond":
                ntypes=2
            elif m["comp"] == "angle":
                ntypes=3
            elif m["comp"] == "dihedral":
                ntypes=4
            else:
                ntypes=0
                continue

            #print(f'{m["comp"]}s:')
            container = dict(ntypes=ntypes, lines={}, type=m["comp"])
            
            components.append(container)
            


        elif regex.match(line) and ntypes>0:
            fields = line.strip().split()
            key = ":".join(fields[:ntypes])

            lines = container["lines"]
            obj = Component(ntypes, line)

            if ntypes==4 and (fields[4]=='9'):
                if key in lines:
                    lines[key] += obj
                else:
                    lines[key] = obj
            else:
                lines[key] = obj
            # print(lines[key])



for component in components:
    print(f'@ffdef ffamber.components :{component["type"]}s TYPE begin')
    components = component["lines"]
    
    for comp in components.values():
        if len(comp) == 1:
            print(comp)
    for comp in components.values():
        if len(comp) > 1:
            print(comp)
#     ntypes = component["ntypes"]

#     for line in component["lines"]:
#         items = line.strip().split(";", 1)
#         fields = items[0].split()
        
#         # types = ["%3s"%("\"%s\""%f) for f in fields[:ntypes]]
#         types = ["%4s"%repr(f).replace("\'","\"") for f in fields[:ntypes]]
#         parms = ['%10s'%f for f in fields[ntypes+1:]]
#         s = ', '.join(types+parms)
#         print(f'    ({s}),  # {"".join(items[1:])}')
    
    print('end\n')