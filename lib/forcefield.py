import lib.handle_file as hf
import lib.exp as myexp

class forcefield:
    def __init__(self, in_file):
        self.file = in_file
        self.find_para = {}
        self.__extract_ffield()

    def __extract_ffield(self):
        blockname = ''
        count = 0
        for word in hf.type_keywords:
            self.find_para[word] = {}
        with open(self.file, 'r') as file:
            for line in file:
                count += 1
                line = hf.pureline(line)
                if line:
                    if line.find('[') != -1:
                        blockname = line[2:-2]
                    elif blockname:
                        linelist = line.split()
                        if blockname == 'defaults':
                            self.ffdefaults = defaults_block(linelist)
                        else:
                            for word in hf.type_keywords:
                                if blockname == word + 'types':
                                    tmp_types_block = types_block(word, linelist)
                                    if self.find_para[word].get(tmp_types_block.id):
                                        try:
                                            if self.find_para[word][tmp_types_block.id].cmp(tmp_types_block):
                                                print('发现重复' + word + 'types: ' + tmp_types_block.id +
                                                      ' at line' + str(count) + ', 但是参数与相同')
                                            else:
                                                raise myexp.retype_Exception(
                                                    word + 'type: ' + tmp_types_block.id + ' at ' + str(count))
                                        except myexp.retype_Exception as e:
                                            print(e)
                                    else:
                                        self.find_para[word][tmp_types_block.id] = tmp_types_block





class defaults_block:
    def __init__(self, linelist):
        self.nbfunc = int(linelist[0])
        self.comb_rule = int(linelist[1])
        self.gen_pairs = linelist[2]
        self.fudgeLJ = float(linelist[3])
        self.fudgeQQ = float(linelist[4])


class types_block:
    def __init__(self, tname, linelist):
        if tname == 'atom':
            self.__nums = 1
            if len(linelist) == 8:
                self.type = linelist[0]
                self.name = linelist[1]
                self.a_nu = int(linelist[2])
                self.mass = float(linelist[3])
                self.charge = float(linelist[4])
                self.ptype = linelist[5]
                self.sigma = float(linelist[6])
                self.epsilon = float(linelist[7])
            else:
                self.type = linelist[0]
                self.name = linelist[0]
                self.a_nu = int(linelist[1])
                self.mass = float(linelist[2])
                self.charge = float(linelist[3])
                self.ptype = linelist[4]
                self.sigma = float(linelist[5])
                self.epsilon = float(linelist[6])
            self.id = self.type
        else:

            if tname == 'bond':
                self.__nums = 2
            elif tname == 'angle':
                self.__nums = 3
            elif tname == 'dihedral':
                self.__nums = 4
            elif tname == 'improper':
                self.__nums = 4
            elif tname == 'pair':
                self.__nums = 2
            elif tname == 'constraint':
                self.__nums = 2

            self.a_name = []
            self.para = []
            self.func = linelist[self.__nums]
            for i in range(self.__nums + 1, len(linelist)):
                self.para.append(float(linelist[i]))
            for i in range(self.__nums):
                self.a_name.append(linelist[i])
            self.id = hf.check_new(self.a_name)

    def cmp(self, v2):
        for i in range(len(v2.para)):
            if self.para[i] != v2.para[i]:
                return False
        return True
