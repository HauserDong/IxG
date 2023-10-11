import os
import numpy as np
from pydrake.geometry.optimization import HPolyhedron

def GcsDir():
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def FindModelFile(filename):
    assert filename[:6] == "models"
    return os.path.join(GcsDir(), filename)

def SerializeEdges(edges, file_path):
    with open(file_path, 'w') as file:
        file.write(len(edges).__str__())
        for edge in edges:
            file.write(',' + edge[0].__str__())
            file.write(',' + edge[1].__str__())

def SerializeRegions(regions, file_path):
    with open(file_path, 'w') as file:
        file.write(len(regions).__str__())
        for region in regions:
            file.write(',' + region.A().shape[0].__str__())
            file.write(',' + region.A().shape[1].__str__())
            for row in region.A():
                for element in row:
                    file.write(',' + element.__str__())
            for element in region.b():
                file.write(',' + element.__str__())


def DeserializeRegions(file_path):
    regions = []
    with open(file_path, 'r') as file:
        data = file.readline().split(',')
        num_regions = int(data[0])

        idx = 1
        for _ in range(num_regions):
            num_rows = int(data[idx])
            num_cols = int(data[idx + 1])
            idx += 2

            A = np.zeros((num_rows, num_cols))
            for i in range(num_rows):
                for j in range(num_cols):
                    A[i, j] = float(data[idx])
                    idx += 1

            b = np.zeros(num_rows)
            for i in range(num_rows):
                b[i] = float(data[idx])
                idx += 1

            regions.append(HPolyhedron(A, b))

    return regions

def DeserializeEdges(file_path):
    edges = []
    with open(file_path, 'r') as file:
        data = file.readline().split(',')
        num_edges = int(data[0])

        idx = 1
        for _ in range(num_edges):
            u = int(data[idx])
            v = int(data[idx + 1])
            idx += 2

            edges.append((u, v))

    return edges
