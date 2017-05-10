from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^result', views.result, name='result'),
    url(r'^hello', views.hello, name='hello'),
    url(r'^random_smiles', views.random_smiles, name='random_smiles'),
    url(r'^random_image', views.random_image, name='random_image'),
    url(r'^structure_image/(?P<smiles>.+)', views.structure_image, name='structure_image'),
]
