from functools import reduce

keywords = ['atom', 'bond', 'angle', 'dihedral', 'improper']  # 针对lammps系统
type_keywords = keywords + ['pair', 'constraint']  # 用于力场
mol_keywords = keywords + ['pair', 'settle', 'exclusion', 'constraint']  # 用于molecule的拓扑信息


def pureline(line):
    return line.replace('\t', ' ').split(';')[0].strip().rstrip('\n').rstrip()


def check_new(n_list):  # 如果不传参数，会导致静态
    if n_list[0] > n_list[-1]:
        n_list.reverse()
    return reduce(lambda x, y: x + '_' + y, n_list)


def extract_block_data(infile, bname):
    blockname = ''
    data = []
    red_flag = False
    with open(infile, 'r') as file:
        for line in file:
            line = pureline(line)
            if line:
                if line.find('[') != -1:
                    blockname = line[2:-2]
                    if blockname == bname:
                        red_flag = True
                elif blockname == bname:
                    linelist = line.split()
                    if bname == 'molecules':
                        data.append(molecules_block(linelist))
                elif red_flag:
                    break
    return data


class molecules_block:
    def __init__(self, linelist):
        self.mol_name = linelist[0]
        self.num = int(linelist[1])
