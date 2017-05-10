#coding=utf-8
from django.shortcuts import render
from django.http import HttpResponse

from .models import Compound
from .forms import SubstructureSearchForm
from django.db.models import Value

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from urllib import unquote


def index(request):
    if request.method == 'POST':
        form = SubstructureSearchForm(request.POST)
        query = request.POST.get('molecule')
        try:
            round_trip_validation(query)

            compounds = Compound.objects.filter(molecule__hassubstruct=Value(query))
            compounds = compounds[:12]
            error_message = None
        except:
            compounds = None
            error_message = 'smilesキーが正しくありません'
    else:
        form = SubstructureSearchForm()
        compounds = None
        error_message = None
    context = {
                'form': form,
                'compounds': compounds,
                'error_message': error_message,
              }

    return render(request, 'tutorial_application/index.html', context)


def round_trip_validation(smiles):
    test_molecule = Chem.MolFromSmiles(smiles)
    test_smiles = Chem.MolToSmiles(test_molecule)
    Chem.MolFromSmiles(test_smiles)


def result(request):
    return HttpResponse('Result')


def hello(request):
    return HttpResponse('Hello World')


def random_smiles(request):
    compounds = Compound.objects.order_by('?')[:15]
    context = {'compounds': compounds, }

    return render(request, 'tutorial_application/random_smiles.html', context)


def random_image(request):
    compounds = Compound.objects.order_by('?')[:15]
    context = {'compounds': compounds, }

    return render(request, 'tutorial_application/random_image.html', context)


def structure_image(request, smiles):
    response = HttpResponse(content_type="image/png")
    smiles = unquote(smiles)
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    image = Draw.MolToImage(mol)
    image.save(response, "PNG")
    return response

