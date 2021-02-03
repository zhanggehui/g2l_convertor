import lib.handle_file as hf
import lib.forcefield as ffd


class molecule:

    def __init__(self, in_name, in_file, ff=None):
        self.name = in_name
        self.file = in_file
        self.ffield = ff
        self.nrexcl = 0
        # 初始化可变对象
        self.mlist_dict = {}
        self.mset_dict = {}
        self.mpobj_dict = {}  # 当作一个有序的set使用
        self.oldlen = {}

        for word in hf.mol_keywords:
            self.mlist_dict[word] = []
            self.mset_dict[word] = set()
            self.mpobj_dict[word] = {}
            self.oldlen[word] = 0

        self.nr2name_list = ['not_use']
        self.nr2type_list = ['not_use']

        self.__extract_moleculetype()

    def __extract_moleculetype(self):
        blockname = ''
        mol_flag = False
        red_flag = False
        with open(self.file, 'r') as file:
            for line in file:
                line = hf.pureline(line)
                if line:
                    if line.find('[') != -1:
                        blockname = line[2:-2]
                    elif blockname == 'moleculetype':
                        mol_name = line.split()[0]
                        if mol_name == self.name:
                            mol_flag = True
                            red_flag = True  # 当后面没有moleculetype时，会失效，但影响很小
                            self.nrexcl = int(line.split()[1])
                        else:
                            mol_flag = False
                    elif mol_flag:
                        linelist = line.split()
                        for word in hf.mol_keywords:
                            if blockname == word + 's':
                                tmp_mp = mol_part(word, linelist, self.nr2name_list, self.nr2type_list, self.ffield)
                                self.mlist_dict[word].append(tmp_mp)
                                self.mset_dict[word].add(tmp_mp.id)
                                if len(self.mset_dict[word]) > self.oldlen[word]:
                                    self.oldlen[word] = len(self.mset_dict[word])
                                    self.mpobj_dict[word][tmp_mp.id] = tmp_mp
                                if blockname == 'atoms':
                                    self.nr2name_list.append(tmp_mp.name)
                                    self.nr2type_list.append(tmp_mp.type)
                    elif red_flag:
                        break


class mol_part:
    def __init__(self, pname, linelist, n2nlist, n2tlist, ff):
        if pname == 'atom':
            self.nr = int(linelist[0])
            self.type = linelist[1]
            self.resnr = int(linelist[2])
            self.res = linelist[3]
            self.name = linelist[4]
            self.cgnr = int(linelist[5])
            self.charge = float(linelist[6])
            if len(linelist) == 8:
                self.mass = float(linelist[7])
            else:
                self.mass = ff.find_para['atom'][self.type].mass
            self.id = linelist[1]
        else:
            if pname == 'bond':
                self.__nums = 2
            elif pname == 'pair':
                self.__nums = 2
            elif pname == 'angle':
                self.__nums = 3
            elif pname == 'dihedral':
                self.__nums = 4
            elif pname == 'improper':
                self.__nums = 4
            elif pname == 'settle':
                pass
            elif pname == 'exclusion':
                pass
            elif pname == 'constraint':
                pass

            self.a_index = []
            self.a_name = []
            for i in range(self.__nums):
                self.a_index.append(int(linelist[i]))
                if ff:
                    an = ff.find_para['atom'][n2tlist[int(linelist[i])]].name
                else:
                    an = n2nlist[int(linelist[i])]
                self.a_name.append(an)
            self.id = hf.check_new(self.a_name)
            if len(linelist) > self.__nums + 1:
                tplist = self.a_name + linelist[self.__nums:]
                tmp_ftp = ffd.types_block(pname, tplist)
                if ff.find_para[pname].get(tmp_ftp.id):
                    if ff.find_para[pname][tmp_ftp.id].cmp(tmp_ftp):
                        print('拓扑文件中类型: '+tmp_ftp.id +'与力场已有中类型: '+tmp_ftp.id +'吻合, 无需处理')
                    else:
                        ff.find_para[pname][tmp_ftp.id] = tmp_ftp
                        print('根据拓扑文件对力场'+ pname + '类型: '+ self.id +'做出修改')
                else:
                    print('添加一个新的' + pname + '类型: ' + 'tmp_ftp.id')
                    ff.find_para[pname][tmp_ftp.id] = tmp_ftp
