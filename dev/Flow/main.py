import os
import sys
from PyQt5.QtWidgets import *

sys.path.append(os.path.abspath("/home/jpereira/project_c/ProtoSyn.jl/dev/Flow"))
from node_editor_window import NodeEditorWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)

    wnd = NodeEditorWindow()

    sys.exit(app.exec_())