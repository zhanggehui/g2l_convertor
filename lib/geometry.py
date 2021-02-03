import numpy as np


# 解析标准格式的gro文件
class geometry:
    def __init__(self, in_file):
        self.file = in_file
        self.title = ''
        self.natoms = 0
        self.box = []
        self.gro_atoms = []
        self._extract_gro_file()

    def _extract_gro_file(self):
        count = 0
        with open(self.file, 'r') as file:
            for line in file:
                count += 1
                if count == 1:
                    self.title = line.lstrip().rstrip('\n').rstrip()
                elif count == 2:
                    self.natoms = int(line[0:5])
                    self.coordinates = np.empty([3, self.natoms])
                elif count == self.natoms + 3:
                    line = line.rstrip('\n')
                    # 斜方的盒子，暂时不做处理
                    if len(line) > 30:
                        pass
                    else:
                        for i in range(3):
                            self.box.append(float(line[i * 10: 10 * (i + 1)]))
                else:
                    ga_tmp = gro_atom(line)
                    self.gro_atoms.append(ga_tmp)
                    for i in range(3):
                        self.coordinates[i, ga_tmp.nr - 1] = ga_tmp.coordinates[i]


class gro_atom:
    def __init__(self, line):
        self.coordinates = []
        self.velocities = []
        self.resnr = int(line[0:5])
        self.res = line[5:10].rstrip()
        self.name = line[10:15].lstrip()
        self.nr = int(line[15:20])
        line = line.rstrip('\n')
        for i in range(3):
            self.coordinates.append(float(line[20 + i * 8: 28 + i * 8]))
            if len(line) > 44:
                self.velocities.append(float(line[44 + i * 8: 52 + i * 8]))
