from __future__ import unicode_literals

from django_rdkit import models

from rdkit import Chem

from urllib import quote

class Compound(models.Model):

    name = models.CharField(max_length=256)
    molecule = models.MolField()

    #torsionbv = models.BfpField(null=True)
    #mfp2 = models.BfpField(null=True)
    #ffp2 = models.BfpField(null=True)

    def __str__(self):
        return '%s: %s' % (self.name, self.get_smiles(), )

    def get_smiles(self):
        return Chem.MolToSmiles(self.molecule)

    def get_image(self):
        encoded_smiles = quote(self.get_smiles())
        src = '/tutorial_application/structure_image/%s' % encoded_smiles
        return '<img src="%s" width=300 height=300>' % src
