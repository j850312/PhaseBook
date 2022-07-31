from django.db import models

# Create your models here.
class Gene(models.Model):
    gene_name = models.CharField(max_length=20)
    tissue = models.CharField(max_length=15)
    age = models.CharField(max_length=50, blank=True)