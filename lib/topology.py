from lib.molecule import molecule
import lib.handle_file as hf


# 单例
class toplogy:
    molecules_list = []  # 分子及其个数
    molecules_set = set()  # 分子种类
    molname2molobj = {}  # 按出现的顺序, 分子名——分子拓扑

    def __init__(self, in_file, sys_ff):
        self.file = in_file
        self.ffield = sys_ff
        self._extract_molecules()

    def _extract_molecules(self):
        oldlen = 0
        self.molecules_list = hf.extract_block_data(self.file, 'molecules')
        for mol in self.molecules_list:
            self.molecules_set.add(mol.mol_name)
            if len(self.molecules_set) > oldlen:
                oldlen = len(self.molecules_set)
                tmp_moltop = molecule(mol.mol_name, self.file, self.ffield)
                self.molname2molobj[tmp_moltop.name] = tmp_moltop
