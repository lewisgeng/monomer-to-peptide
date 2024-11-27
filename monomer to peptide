from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import molzip
from rdkit.Chem import Draw, rdmolfiles
Draw.SetComicMode = True
Draw.rdDepictor.SetPreferCoordGen(True)
import copy, sys, django, os,re, rdkit, io, base64
from django.conf import settings
from PIL import Image

def generate_peptide(monomer_mol_list):
    CXSMILES_LIST = []
    SMILES_LIST = []
    def replace_inside_brackets(text, replacement):
        pattern = r'\[([^\]]+)\]'
        return re.sub(pattern, replacement, text)

    ps = Chem.SmilesWriteParams()
    for p in monomer_mol_list:
        mol = Chem.MolFromMolBlock(p)
        cxsmiles = Chem.MolToCXSmiles(mol, ps, Chem.rdmolfiles.CXSmilesFields.CX_ALL)
        cxsmiles=replace_inside_brackets(cxsmiles, "[*]")
        CXSMILES_LIST.append(cxsmiles)


    for cxsmiles in CXSMILES_LIST:
        comma_list = []
        cx_smiles_header = cxsmiles.split('|')[0]
        cx_smiles_body = cxsmiles.split('|')[1]
        coordinate = cx_smiles_body.split('(')[1].split(')')[0]
        r_info = cx_smiles_body.split('(')[1].split(')')[1].strip(',')
        f_r_index = r_info.split(':')[1].split('.')[0]
        f_r = r_info.split(':')[1].split('.')[2]
        if len(r_info.split(':'))>2:
            s_r_index = r_info.split(':')[2].split('.')[0]
            s_r = r_info.split(':')[2].split('.')[2]
        for n in range(coordinate.count(';')):
            comma_list.append(';')
        comma_list.insert(int(f_r_index), f_r)
        if len(r_info.split(':'))>2:
            comma_list.insert(int(s_r_index)+1, s_r)
        #SMILES_LIST = ['[*]CCCC(C)c1ccc([*])cc1 |$R1;;;;;;;;;;R2;;;$|','[*]C1CCC(C[*])CC1 |$R1;;;;;;R2;;$|']
        HELM = '%s|$%s$|' % (cx_smiles_header,''.join(comma_list))
        SMILES_LIST.append(HELM)
    mols = []
    for s in SMILES_LIST:
        m=Chem.MolFromSmiles(s)
        mols.append(m)

    def combine_fragments(m1, m2):
        m1 = Chem.Mol(m1)
        m2 = Chem.Mol(m2)
        for atm in m1.GetAtoms():
            if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == 'R2':
                atm.SetAtomMapNum(1)
        for atm in m2.GetAtoms():
            if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == 'R1':
                atm.SetAtomMapNum(1)
        return molzip(m1, m2)

    def make_peptide(monomerlist):
        monomerlist = copy.deepcopy(monomerlist)
        for idx, monomer in enumerate(monomerlist):
            if Chem.MolToSmiles(monomer).count("*") == 1:
                continue
            if idx==0:
                res = monomer
            else:
                res = combine_fragments(res, monomer)
        return res

    def cap_terminal(m):
        m = Chem.Mol(m)
        n_term = Chem.MolFromSmiles('CC(=O)[*:1]')
        c_term = Chem.MolFromSmiles('CO[*:2]')
        for atm in m.GetAtoms():
            if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == 'R1':
                atm.SetAtomMapNum(1)
            if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == 'R2':
                atm.SetAtomMapNum(2)
        res = molzip(m, n_term)
        res = molzip(res, c_term)
        return res

    #m1 = combine_fragments(mols[0], mols[1])
    #m2 = combine_fragments(mols[1], mols[0])
    #img=Draw.MolsToGridImage([mols[0], mols[1], m1, m2], molsPerRow=4)
    #img.save('1.png')

    m_f = make_peptide(mols)
    return_mol = Chem.MolToMolBlock(cap_terminal(m_f))

    img=Draw.MolsToGridImage([cap_terminal(m_f)], molsPerRow=1, subImgSize=(800,600))


    # 将PIL图像转换为PNG格式，并写入内存
    img_bytes = io.BytesIO()
    img.save(img_bytes, format='PNG')
    img_bytes.seek(0)

    # 读取内存中的图像数据，并转换为Base64字符串
    image_data = base64.b64encode(img_bytes.read()).decode('utf-8')

    return image_data,return_mol


