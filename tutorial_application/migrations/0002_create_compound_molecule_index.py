# -*- coding: utf-8 -*-
# Generated by Django 1.10.5 on 2017-03-09 01:55
from __future__ import unicode_literals

from django.db import migrations
from django_rdkit.operations import GiSTIndex


class Migration(migrations.Migration):

    dependencies = [
        ('tutorial_application', '0001_initial'),
    ]

    operations = [
        GiSTIndex('Compound', 'molecule')
    ]
