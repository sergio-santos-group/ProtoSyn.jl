from node_graphics_edge import QDMGraphicsEdge

class SceneHistory():
    def __init__(self, scene):
        self.scene = scene

        self.history_stack = []
        self.history_current_step = -1
        self.history_limit = 32

    def undo(self):
        if self.history_current_step > 0:
            self.history_current_step -= 1
            self.restoreHistory()

    def redo(self):
        if self.history_current_step + 1 < len(self.history_stack):
            self.history_current_step += 1
            self.restoreHistory()


    def restoreHistory(self):
        self.restoreHistoryStamp(self.history_stack[self.history_current_step])


    def storeHistory(self, desc):
        # if the pointer (history_current_step) is not at the end of history_stack
        if self.history_current_step+1 < len(self.history_stack):
            self.history_stack = self.history_stack[0:self.history_current_step+1]

        # history is outside of the limits
        if self.history_current_step+1 >= self.history_limit:
            self.history_stack = self.history_stack[1:]
            self.history_current_step -= 1

        hs = self.createHistoryStamp(desc)

        self.history_stack.append(hs)
        self.history_current_step += 1


    def createHistoryStamp(self, desc):
        sel_obj = {
            'nodes': [],
            'edges': [],
        }
        for item in self.scene.grScene.selectedItems():
            if hasattr(item, 'node'):
                sel_obj['nodes'].append(item.node.id)
            elif isinstance(item, QDMGraphicsEdge):
                sel_obj['edges'].append(item.edge.id)

        history_stamp = {
            'desc': desc,
            'snapshot': self.scene.serialize(),
            'selection': sel_obj,
        }

        print("Returning history stamp:", history_stamp)
        return history_stamp


    def restoreHistoryStamp(self, history_stamp):
        self.scene.deserialize(history_stamp['snapshot'])

        # restore selection
        for edge_id in history_stamp['selection']['edges']:
            for edge in self.scene.edges:
                if edge.id == edge_id:
                    edge.grEdge.setSelected(True)
                    break

        for node_id in history_stamp['selection']['nodes']:
            for node in self.scene.nodes:
                if node.id == node_id:
                    node.grNode.setSelected(True)
                    break

