from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

from node_scene import Scene
from node_node import Node
from node_edge import Edge, EDGE_TYPE_BEZIER
from node_graphics_view import QDMGraphicsView


class NodeEditorWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.stylesheet_filename = 'qss/nodestyle.qss'
        self.loadStylesheet(self.stylesheet_filename)

        self.initUI()



    def initUI(self):
        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        # crate graphics scene
        self.scene = Scene()

        # Init Nodes
        self.addNodes()

        # create graphics view
        self.view = QDMGraphicsView(self.scene.grScene, self)
        self.layout.addWidget(self.view)

        # Save initial status for Redo
        self.scene.history.storeHistory("Start")


    def addNodes(self):
        node1 = Node(self.scene, "Monte Carlo Driver", inputs=[1,2,3], outputs=[4])
        node2 = Node(self.scene, "Dihedral Mutator", inputs=[3], outputs=[1])
        node3 = Node(self.scene, "Energy Function", inputs=[], outputs=[2])
        node1.setPos(200, -150)
        node2.setPos(-350, 0)
        node3.setPos(-75, -250)

        edge1 = Edge(self.scene, node3.outputs[0], node1.inputs[1], edge_type = EDGE_TYPE_BEZIER)
        edge2 = Edge(self.scene, node2.outputs[0], node1.inputs[0], edge_type = EDGE_TYPE_BEZIER)


    def loadStylesheet(self, filename):
        file = QFile(filename)
        file.open(QFile.ReadOnly | QFile.Text)
        stylesheet = file.readAll()
        QApplication.instance().setStyleSheet(str(stylesheet, encoding='utf-8'))
