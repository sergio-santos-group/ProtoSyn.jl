#!/usr/bin/env python
# -*- coding: utf-8 -*-
from xmlrpc.server import SimpleXMLRPCServer
import torchani
import torch

print("Starting TorchANI XML-RPC server ...")
address = ('localhost', 50000)
server = SimpleXMLRPCServer(address)

class TorchANI(object):
    """Provide the TorchANI calculation functions"""

    def __init__(self):

        if torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")

        self.model = torchani.models.ANI2x(periodic_table_index = True).to(self.device)

    def calc(self, s, c, update_forces = False, model_index = 3):
        coordinates = torch.tensor([c], requires_grad = True, device = self.device).float()
        species     = torch.tensor([s], device = self.device)
        m1 = self.model.species_converter((species, coordinates))
        m2 = self.model.aev_computer(m1)
        m3 = self.model.neural_networks[model_index](m2)[1]
        if update_forces:
            f = torch.autograd.grad(m3.sum(), coordinates)[0]
            return m3.item(), f.tolist()
        else:
            return m3.item(), "none"

server.register_instance(TorchANI())
print("TorchANI XML-RPC Server is online!")

if __name__ == '__main__':
    try:
        print("TorchANI Python server running on %s:%s" % address)
        print("Use Ctrl-C to Exit")
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()
        print("Exiting")