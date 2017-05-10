from django_rdkit import models

from django.forms.models import ModelForm

from .models import Compound


class SubstructureSearchForm(ModelForm):
    class Meta:
        model = Compound
        fields = ('molecule', )
