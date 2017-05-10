"""1. Creation of the tutorial project"""
from django.db import connection

with connection.cursor() as cursor:
    cursor.execute("SELECT mol_amw('C')")
    print(cursor.fetchone()[0])


"""2. Creation of a django application"""

from tutorial_application.models import Compound

Compound.objects.create(name='benzene', molecule='c1ccccc1')

from django_rdkit.models import *

for compound in Compound.objects.annotate(amw=AMW('molecule')):
    print(compound.name, compound.amw)


Compound.objects.all().delete()


"""3. Structures import and substructure queries"""

path = 'chembl16_30K.txt'

from rdkit import Chem

def chembl(path, limit=None):
    count = 0
    with open(path, 'rt') as f:
        for line in f:
            name, smiles = line.split()[:2]
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:
                yield name, molecule
                count += 1
                if limit and count == limit:
                    break

from tutorial_application.models import Compound

for name, molecule in chembl(path, limit=None):
    smiles = Chem.MolToSmiles(molecule)
    test_molecule = Chem.MolFromSmiles(smiles)
    if not test_molecule:
        print('smiles-mol-smiles roundtrip issue:', name)
    else:
        Compound.objects.create(name=name, molecule=molecule)

Compound.objects.count()

""""""

from django_rdkit.models import *

from tutorial_application.models import *

def smiles_substructure_query(substructure):
    query = Compound.objects.filter(molecule__hassubstruct=substructure)
    for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
        print(cmpd.name, cmpd.smiles)

smiles_substructure_query('c1cccc2c1nncc2')

from django.db.models import Value

def smiles_substructure_query(substructure):
    query = Compound.objects.filter(molecule__hassubstruct=Value(substructure))
    for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
        print(cmpd.name, cmpd.smiles)

smiles_substructure_query('c1cccc2c1nncc2')
smiles_substructure_query('c1ccccc1O')


"""4. SMARTS-based queries"""

def smarts_substructure_query(substructure):
    query = Compound.objects.filter(molecule__hassubstruct=QMOL(Value(substructure)))
    for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
        print(cmpd.name, cmpd.smiles)

smarts_substructure_query('c1[o,s]ncn1')


"""5. Using stereochemistry"""

smiles_substructure_query('NC(=O)[C@H]1CCCN1C=O')

from django_rdkit.config import config

config.do_chiral_sss = True

smiles_substructure_query('NC(=O)[C@H]1CCCN1C=O')


"""6. Similarity queries"""

from django_rdkit.models import *

from django_rdkit_tutorial.models import Compound

Compound.objects.update(
    torsionbv=TORSIONBV_FP('molecule'),
    mfp2=MORGANBV_FP('molecule'),
    ffp2=FEATMORGANBV_FP('molecule'),
)

""""""

from django_rdkit.models import *

from tutorial_application.models import Compound

smiles = 'Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1'

value = MORGANBV_FP(Value(smiles))

Compound.objects.filter(mfp2__tanimoto=value).count()

def get_mfp2_neighbors(smiles):
    value = MORGANBV_FP(Value(smiles))
    queryset = Compound.objects.filter(mfp2__tanimoto=value)
    queryset = queryset.annotate(smiles=MOL_TO_SMILES('molecule'))
    queryset = queryset.annotate(sml=TANIMOTO_SML('mfp2', value))
    queryset = queryset.order_by(TANIMOTO_DIST('mfp2', value))
    queryset = queryset.values_list('name', 'smiles', 'sml')
    return queryset


for name, smiles, sml in get_mfp2_neighbors('CN(C)c1c(-c2ccccc2)c2ccccc2[nH]c1=O')[:10]:
    print(name, smiles, sml)


"""7. Adjusting the similarity cutoff"""

print(get_mfp2_neighbors('CN(C)c1c(-c2ccccc2)c2ccccc2[nH]c1=O').count())

from django_rdkit.config import config

config.tanimoto_threshold = 0.7

print(get_mfp2_neighbors('CN(C)c1c(-c2ccccc2)c2ccccc2[nH]c1=O').count())

config.tanimoto_threshold = 0.6

print(get_mfp2_neighbors('CN(C)c1c(-c2ccccc2)c2ccccc2[nH]c1=O').count())

config.tanimoto_threshold = 0.5

print(get_mfp2_neighbors('CN(C)c1c(-c2ccccc2)c2ccccc2[nH]c1=O').count())