from collections import OrderedDict
from node_serializable import Serializable
from node_graphics_socket import QDMGraphicsSocket

LEFT_TOP     = 1
LEFT_CENTER  = 2
LEFT_BOTTOM  = 3
RIGHT_TOP    = 4
RIGHT_CENTER = 5
RIGHT_BOTTOM = 6



class Socket(Serializable):
    def __init__(self, node, index=0, position=LEFT_TOP, socket_type=1, multi_edges=True, count_on_this_node_side=1, is_input=False):
        super().__init__()

        self.node                    = node
        self.index                   = index
        self.position                = position
        self.socket_type             = socket_type
        self.is_multi_edges          = multi_edges
        self.count_on_this_node_side = count_on_this_node_side
        self.is_input                = is_input
        self.is_output               = not self.is_input


        self.grSocket = QDMGraphicsSocket(self, self.socket_type)

        self.setSocketPosition()

        self.edges = []


    def __str__(self):
        return "<Socket %s %s..%s>" % ("ME" if self.is_multi_edges else "SE", hex(id(self))[2:5], hex(id(self))[-3:])

    def setSocketPosition(self):
        self.grSocket.setPos(*self.node.getSocketPosition(self.index, self.position, self.count_on_this_node_side))

    def getSocketPosition(self):
        res = self.node.getSocketPosition(self.index, self.position, self.count_on_this_node_side)
        return res

    def addEdge(self, edge):
        self.edges.append(edge)

    def removeEdge(self, edge):
        if edge in self.edges: self.edges.remove(edge)

    def removeAllEdges(self):
        while self.edges:
            edge = self.edges.pop(0)
            edge.remove()

    def determineMultiEdges(self, data):
        if 'multi_edges' in data:
            return data['multi_edges']
        else:
            # probably older version of file, make RIGHT socket multiedged by default
            return data['position'] in (RIGHT_BOTTOM, RIGHT_TOP)

    def serialize(self):
        return OrderedDict([
            ('id', self.id),
            ('index', self.index),
            ('multi_edges', self.is_multi_edges),
            ('position', self.position),
            ('socket_type', self.socket_type),
        ])

    def deserialize(self, data, hashmap={}, restore_id=True):
        if restore_id: self.id = data['id']
        self.is_multi_edges = self.determineMultiEdges(data)
        hashmap[data['id']] = self
        return True