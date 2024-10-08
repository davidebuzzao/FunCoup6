from django.contrib import admin

# Register your models here.
from .models import *

admin.site.register(Species)
admin.site.register(Proteome)
admin.site.register(Protein)
admin.site.register(GoldStandard)
admin.site.register(Evidence)