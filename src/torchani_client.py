from xmlrpc import client
import numpy as np
import time

proxy = client.ServerProxy("http://127.0.0.1", verbose=False)
s = np.array([6, 1, 1, 1, 1])
c = np.array([[ 0.03192167,  0.00638559,  0.01301679],
              [-0.83140486,  0.39370209, -0.26395324],
              [-0.66518241, -0.84461308,  0.20759389],
              [ 0.45554739,  0.54289633,  0.81170881],
              [ 0.66091919, -0.16799635, -0.91037834]])
start = time.time()
a = proxy.calc(s.tolist(), c.tolist())
end = time.time()
print("RESULT XMLRPC:", a[0])
print(end-start, "s\n")

import torch
import torchani

device = torch.device("cuda")
model = torchani.models.ANI2x(periodic_table_index = True).to(device)
update_forces = False

start = time.time()
coordinates = torch.tensor([c], requires_grad = True, device = device).float()
species     = torch.tensor([s], device = device)
m1 = model.species_converter((species, coordinates))
m2 = model.aev_computer(m1)
m3 = model.neural_networks[3](m2)[1]
# return m3.item()
if update_forces:
    f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
    r = m3.item(), f.cpu().numpy()
else:
    r = m3.item(), "none"
end = time.time()
print("RESULT SERIAL:", r[0])
print(end-start, "s\n")

start = time.time()
coordinates = torch.tensor([c], requires_grad = True, device = device).float()
species     = torch.tensor([s], device = device)
m1 = model.species_converter((species, coordinates))
m2 = model.aev_computer(m1)
m3 = model.neural_networks[3](m2)[1]
# return m3.item()
if update_forces:
    f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
    r = m3.item(), f.cpu().numpy()
else:
    r = m3.item(), "none"
end = time.time()
print("RESULT SERIAL:", r[0])
print(end-start, "s\n")