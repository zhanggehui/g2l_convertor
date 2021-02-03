from lib.handle_file import keywords
import lib.topology as top
import lib.geometry as gro
import lib.forcefield as ffield
import sys

"""
lammps需要的信息:
    atoms, bonds, angles, dihedrals, impropers
    atom types, bond types, angle types, dihedral types, improper types
    xlo xhi ylo yhi zlo zhi
    Masses
    Pair Coeffs
    Atoms
    Bond Coeff, Bonds
    Angle Coeffs, Angles
    Dihedral Coeffs, Dihedrals
    Improper Coeffs, Impropers
    
注意单位换算，距离单位要×10
原子用拓扑文件中的信息，键用力场中的信息
"""


class lammps:
    num_list = []  # atoms, bonds, angles, dihedrals, impropers
    name2obj = {}  # 当作一个有序的set使用
    name2lmptype = {}

    def __init__(self, in_gro, in_top, ff, outfile):
        self.gro = in_gro
        self.top = in_top
        self.ffield = ff
        self.out_file = outfile

        for word in keywords:
            self.name2obj[word] = {}
            self.name2lmptype[word] = {}

        for word in keywords:
            total_num = 0
            for mol in self.top.molecules_list:
                moltop = self.top.molname2molobj[mol.mol_name]
                total_num += mol.num * len(moltop.mlist_dict[word])
            self.num_list.append(total_num)

        for word in keywords:
            for moltop in self.top.molname2molobj.values():
                for key, value in moltop.mpobj_dict[word].items():
                    self.name2obj[word][key] = self.ffield.find_para[word][key]

    def write_data(self):
        e_ratio=4.18585182085
        with open(self.out_file, 'w') as f:
            f.write('%s\n\n' % self.gro.title)

            # 数目
            for i in range(len(keywords)):
                str_tmp = keywords[i] + 's'
                if self.num_list[i] > 0:
                    f.write('%18d  %-s\n' % (self.num_list[i], str_tmp))
            f.write('\n')

            # 各种类型的数目
            i = 0
            for word in keywords:
                str_tmp = word + ' types'
                if self.num_list[i] > 0:
                    f.write('%18d  %-s\n' % (len(self.name2obj[word]), str_tmp))
                i += 1
            f.write('\n')

            # 盒子尺寸
            xyz = ['x', 'y', 'z']
            for i in range(3):
                str_tmp = xyz[i] + 'lo ' + xyz[i] + 'hi'
                if i!=2:
                    f.write('%12.4f %12.4f %s\n' % (0., 10 * self.gro.box[i], str_tmp))
                else :
                    f.write('%12.4f %12.4f %s\n' % (0.-5, 10 * self.gro.box[i]+5, str_tmp))

            # 质量
            atomtys = self.ffield.find_para['atom']
            f.write('\nMasses\n\n')
            i = 0
            for key, value in self.name2obj['atom'].items():
                i += 1
                f.write('%8d%14.3f  # %-s, %s\n' % (i, value.mass, key, atomtys[key].name))
            f.write('\n')

            list_i = 0
            totalcg=0
            for word in keywords:
                list_i += 1
                if self.num_list[list_i - 1] > 0:
                    title_str = word.capitalize() + ' Coeffs'
                    if word == 'atom':
                        title_str = 'Pair Coeffs'
                    f.write('\n%s\n\n' % title_str)
                    count_tp = 0
                    for key, mp_tp in self.name2obj[word].items():
                        count_tp += 1
                        self.name2lmptype[word][key] = count_tp
                        if word == 'atom':
                            a_tmp = atomtys[key]
                            f.write('%8d%17.14f%17.14f # %-s, %s\n' % (
                                count_tp, a_tmp.epsilon/e_ratio, 10*a_tmp.sigma, key, a_tmp.name))
                        elif word == 'bond':
                            value_list = mp_tp.para
                            f.write('%12d%9.3f%9.3f #%s\n' % (count_tp, 0.5*value_list[1]/(e_ratio*100), 10*value_list[0], mp_tp.id))
                        elif word == 'angle':
                            value_list = mp_tp.para
                            f.write('%12d%11.5f%11.5f #%s\n' % (count_tp, 0.5*value_list[1]/e_ratio, value_list[0], mp_tp.id))

                    title_str = word.capitalize() + 's'
                    f.write('\n%s\n\n' % title_str)
                    count_p = 0
                    count_acc = 0
                    for mol in self.top.molecules_list:
                        moltop = self.top.molname2molobj[mol.mol_name]
                        molplist = moltop.mlist_dict[word]
                        mol_natoms = len(moltop.mlist_dict['atom'])
                        for i in range(mol.num):
                            for mp in molplist:
                                count_p += 1
                                if word == 'atom':
                                    gatom = self.gro.gro_atoms[count_p - 1]
                                    _resnr = gatom.resnr
                                    _lmptype = self.name2lmptype['atom'][mp.type]
                                    _charge = mp.charge
                                    _x = gatom.coordinates[0]
                                    _y = gatom.coordinates[1]
                                    _z = gatom.coordinates[2]
                                    f.write('%12d%9d%7d %10g%17.10f%17.10f%17.10f #%s\n' % (
                                        count_p, _resnr, _lmptype, _charge, 10 * _x, 10 * _y, 10 * _z, mp.type))
                                    totalcg+=_charge
                                elif word == 'bond':
                                    _lmptype = self.name2lmptype[word][mp.id]
                                    f.write('%12d%7d%10d%10d\n' % (
                                        count_p, _lmptype, mp.a_index[0] + count_acc, mp.a_index[1] + count_acc))
                                elif word == 'angle':
                                    _lmptype = self.name2lmptype[word][mp.id]
                                    f.write('%12d%7d%10d%10d%10d\n' % (
                                        count_p, _lmptype, mp.a_index[0] + count_acc, mp.a_index[1] + count_acc,
                                        mp.a_index[2] + count_acc))
                            count_acc += mol_natoms

            print('体系净电荷为: %g' %totalcg)

def gmx2lmp():
    if len(sys.argv) >= 4:
        grofile = sys.argv[1]
        pp_top = sys.argv[2]
        outfile = sys.argv[3]
    else:
        grofile = r'E:\Share_ubuntu\GO_MD\GO_build\NA_pool\0.7\GO_pool_0.7.gro'
        pp_top = r'E:\Share_ubuntu\GO_MD\GO_build\NA_pool\0.7\GO_pool_0.7_pp.top'
        outfile = r'E:\Share_ubuntu\GO_MD\GO_build\NA_pool\0.7\GO.data'
    sys_ff = ffield.forcefield(pp_top)
    sys_gmx_geo = gro.geometry(grofile)
    sys_gmx_top = top.toplogy(pp_top, sys_ff)

    sys_lmp = lammps(sys_gmx_geo, sys_gmx_top, sys_ff, outfile)
    sys_lmp.write_data()


if __name__ == '__main__':
    gmx2lmp()
