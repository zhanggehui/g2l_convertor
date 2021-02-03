import os
import lib.geometry as gro
import lib.molecule as mol
import lib.handle_file as hf
from scipy.stats import mode
import numpy as np
import math
import sys

if len(sys.argv) > 1:
    GO_index = sys.argv[1]
else:
    GO_index = ''

# 以下为用到的全局变量
grp_names = ['GRA', 'EPO', 'HYD', 'CAR']
grp_anum = [1, 1, 2, 4]
first_EPO = 0
gra_y_range = []
CCAR_y_pos = []


def renameGO():
    print('---------------------------  Renaming GO  ---------------------------')
    grofile = 'reboxGO.gro'
    newgrofile = 'renameGO.gro'

    geo = gro.geometry(grofile)

    index_dict = {}
    for word in grp_names:
        index_dict[word] = []

    index_dict['GRA'].append(0)
    global first_EPO
    first_HYD = 0
    first_CAR = 0
    count = 0
    for atom in geo.gro_atoms:
        if atom.name == 'O' and first_EPO == 0:
            index_dict['GRA'].append(count)
            index_dict['EPO'].append(count)
            first_EPO = count
        elif atom.name == 'H' and first_HYD == 0:
            index_dict['EPO'].append(count - 1)
            index_dict['HYD'].append(count - 1)
            first_HYD = count
        elif atom.name == 'C' and first_HYD != 0 and first_CAR == 0:
            index_dict['HYD'].append(count)
            index_dict['CAR'].append(count)
            first_CAR = count
        count += 1
    index_dict['CAR'].append(len(geo.gro_atoms))

    gra_x_range = [np.min(geo.coordinates[0, index_dict['GRA'][0]:index_dict['GRA'][1]]),
                   np.max(geo.coordinates[0, index_dict['GRA'][0]:index_dict['GRA'][1]])]
    CCAR_coordy = []
    newnatoms = 0
    for word in grp_names:
        for i in range(index_dict[word][0], index_dict[word][1]):
            atom = geo.gro_atoms[i]
            if word == 'GRA':
                newnatoms += 1
                atom.name = 'CGRA'
            elif word == 'EPO':
                newnatoms += 1
                atom.name = 'OEPO'
            elif word == 'HYD':
                newnatoms += 1
                if (i - index_dict[word][0]) % 2 == 0:
                    atom.name = 'OHYD'
                elif (i - index_dict[word][0]) % 2 == 1:
                    atom.name = 'HHYD'
            elif word == 'CAR':
                if between(atom.coordinates[0], gra_x_range):
                    newnatoms += 1
                    if (i - index_dict[word][0]) % 4 == 0:
                        atom.name = 'CCAR'
                        CCAR_coordy.append(atom.coordinates[1])
                    elif (i - index_dict[word][0]) % 4 == 1:
                        atom.name = 'O1CAR'
                    elif (i - index_dict[word][0]) % 4 == 2:
                        atom.name = 'O2CAR'
                    elif (i - index_dict[word][0]) % 4 == 3:
                        atom.name = 'HCAR'
                else:
                    atom.name = 'Not'

    gra_y_range.extend([np.min(geo.coordinates[1, index_dict['GRA'][0]:index_dict['GRA'][1]]),
                        np.max(geo.coordinates[1, index_dict['GRA'][0]:index_dict['GRA'][1]])])
    CCAR_y_pos.extend([min(CCAR_coordy), max(CCAR_coordy)])

    with open(newgrofile, 'w') as f:
        f.write('GO\n')
        f.write('%5d\n' % newnatoms)
        count = 0
        for atom in geo.gro_atoms:
            if atom.name != 'Not':
                count += 1
                f.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %
                        (atom.resnr,
                         atom.res,
                         atom.name,
                         count,
                         atom.coordinates[0],
                         atom.coordinates[1],
                         atom.coordinates[2])
                        )
        f.write('%10.5f%10.5f%10.5f\n' % (geo.box[0], geo.box[1], geo.box[2]))
    if count != newnatoms:
        print('两次原子计数不相同,请检查!')
        sys.exit()
    else:
        print('一切顺利')


def between(x, x_list):
    return x_list[0] < x < x_list[1]


def resizeGO():
    print('---------------------------  Resizing GO  ---------------------------')
    grofile = 'renameGO.gro'
    newgrofile = 'resizeGO.gro'
    ratio = 1.42 / 3

    grpcount_dict = {}
    grpanum_dict = {}

    for i in range(len(grp_names)):
        grpcount_dict[grp_names[i]] = 0
        grpanum_dict[grp_names[i]] = grp_anum[i]

    count = 0
    with open(grofile, 'r') as grof:
        with open(newgrofile, 'w') as newgro:
            for line in grof:
                count += 1
                if count == 1:
                    newgro.write(line)
                elif count == 2:
                    natoms = float(line[0:5])
                    newgro.write(line)
                elif count == natoms + 3:
                    x = ratio * float(line[0:10])
                    y = ratio * float(line[10:20])
                    z = float(line[20:30])
                    newgro.write('%10.5f%10.5f%10.5f\n' % (x, y, z))
                else:
                    grpnm = line[12:15]
                    grpcount_dict[grpnm] += 1
                    if grpcount_dict[grpnm] % grpanum_dict[grpnm] == 1 or grpanum_dict[grpnm] == 1:
                        c_center = getcoordinates(line)
                        if grpnm != 'CAR':
                            A = np.array([ratio * c_center[0], ratio * c_center[1], c_center[2]])
                        else:
                            ind = -1
                            for i in range(len(gra_y_range)):
                                if c_center[1] == CCAR_y_pos[i]:
                                    ind = i
                            if ind == -1:
                                print('CCAR位置异常，此算法不对')
                                sys.exit()
                            else:
                                A = np.array([ratio * c_center[0],
                                              ratio * gra_y_range[ind] + (CCAR_y_pos[ind] - gra_y_range[ind]),
                                              c_center[2]])
                        C = A
                    else:
                        c_a = getcoordinates(line)
                        if grpnm != 'CAR':
                            B = A + np.array(c_a) - np.array(c_center)
                        else:
                            if grpcount_dict[grpnm] % grpanum_dict[grpnm] == 2:
                                tz = np.random.choice([math.pi / 2, -math.pi / 2])
                            rotate = np.array([[math.cos(tz), 0, -math.sin(tz)],
                                               [0, 1, 0],
                                               [-math.sin(tz), 0, math.cos(tz)]])
                            B = A + np.dot(rotate, np.array(c_a) - np.array(c_center))
                        C = B
                    left = line[0:20]
                    newgro.write('%s%8.3f%8.3f%8.3f\n' % (left, C[0], C[1], C[2]))


def getcoordinates(line):
    coordinates = []
    for i in range(3):
        coordinates.append(float(line[20 + i * 8: 28 + i * 8]))
    return coordinates


def make_itp():
    grofile = 'resizeGO.gro'
    moltop = 'n2tGO.top'
    outitp = 'GO.itp'
    outtypesitp = 'GOtypes.itp'

    # cg = np.random.randint(20, 30, 1)
    cg = 0
    avecharge = cg / first_EPO

    geo = gro.geometry(grofile)
    mole = mol.molecule('GO', moltop)

    k_dict = {'bond': 400000, 'angle': 400}
    func_dict = {'bond': 1, 'pair': 1, 'angle': 1}
    para_dict = {}
    for word in hf.mol_keywords[1:3]:
        para_dict[word] = {}
        for mp_tp in mole.mset_dict[word]:
            para_dict[word][mp_tp] = []

    for word in hf.mol_keywords[1:3]:
        for mp in mole.mlist_dict[word]:
            if word == 'bond':
                c1 = np.array(geo.gro_atoms[mp.a_index[0] - 1].coordinates)
                c2 = np.array(geo.gro_atoms[mp.a_index[1] - 1].coordinates)
                bondlength = np.linalg.norm(c1 - c2)
                bondlength = float('%.3f' % bondlength)
                if bondlength > 0.3:
                    print('键长太长!周期性边界出现问题!')
                    sys.exit(-1)
                para_dict[word][mp.id].append(bondlength)
                mp.value = bondlength
            elif word == 'angle':
                c1 = np.array(geo.gro_atoms[mp.a_index[0] - 1].coordinates)
                c2 = np.array(geo.gro_atoms[mp.a_index[1] - 1].coordinates)
                c3 = np.array(geo.gro_atoms[mp.a_index[2] - 1].coordinates)
                u = c1 - c2
                v = c3 - c2
                costheta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
                theta = float('%.1f' % math.degrees(math.acos(costheta)))
                para_dict[word][mp.id].append(theta)
                mp.value = theta

    with open(outtypesitp, 'w') as f:
        for word in hf.mol_keywords[1:3]:
            f.write('\n[ ' + word + 'types ]\n')
            for key, value in para_dict[word].items():
                keylist = key.split('_')
                for i in range(len(keylist)):
                    f.write('%s   ' % (keylist[i]))
                para_array = np.array(list(value))
                ave = np.mean(para_array)
                if np.max(para_array) - np.min(para_array) > 0.1 * ave:
                    tmp = para_array[para_array > ave]
                    avepara = mode(tmp)[0][0]
                else:
                    avepara = mode(para_array)[0][0]
                f.write('%d   %10.3f   %10.1f\n' % (func_dict[word], avepara, k_dict[word]))

                mobj_tmp = mole.mpobj_dict[word][key]
                mobj_tmp.para = avepara

    with open(outitp, 'w') as f:
        f.write('\n[ moleculetype ]\n')
        f.write('; Name            nrexcl     charge = %d\n' % cg)
        f.write('%-10s%7d\n' % ('GO' + GO_index, mole.nrexcl))

        f.write('\n[ atoms ]\n')
        f.write(';   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n')
        for atom in mole.mlist_dict['atom']:
            if atom.name != 'CGRA':
                avecharge = 0
            f.write('%6d%11s%7d%7s%7s%7d %10g%11.3f\n' %
                    (atom.nr,
                     atom.type,
                     atom.resnr,
                     atom.res,
                     atom.name,
                     atom.cgnr,
                     atom.charge - avecharge,
                     atom.mass))

        f.write('\n[ bonds ]\n')
        f.write(';  ai    aj funct            c0            c1            c2            c3\n')
        for bond in mole.mlist_dict['bond']:
            f.write('%5d %5d %5d\n' % (bond.a_index[0], bond.a_index[1], func_dict['bond']))

        f.write('\n[ pairs ]\n')
        f.write(';  ai    aj funct            c0            c1            c2            c3\n')
        for pr in mole.mlist_dict['pair']:
            f.write('%5d %5d %5d\n' % (pr.a_index[0], pr.a_index[1], func_dict['pair']))

        f.write('\n[ angles ]\n')
        f.write(';  ai    aj    ak funct            c0            c1            c2            c3\n')
        for ag in mole.mlist_dict['angle']:
            mobj_tmp = mole.mpobj_dict['angle'][ag.id]
            para_k = mobj_tmp.para
            if ag.value > 0.9 * para_k:
                f.write('%5d %5d %5d %5d\n' % (ag.a_index[0], ag.a_index[1], ag.a_index[2], func_dict['angle']))
            else:
                print('An unnecessary angel, type: ' + ag.id)


def main():
    renameGO()
    os.system('gmx x2top -f renameGO.gro -o n2tGO.top -ff oplsaaGO -name GO -nopbc -noparam -pairs -nexcl 3')
    resizeGO()
    make_itp()


if __name__ == "__main__":
    main()
