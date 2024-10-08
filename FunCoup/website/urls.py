from django.urls import path,include
from . import views

# Here, all website urls are defined.
# A url is linked to a function on the views.py script, and is given a name
# Where multiple urls are defined for the same view, this means that some of the query parameters are optional
 
urlpatterns=[
    path('', views.index, name="index"),
    path('search/', views.index, name="index"),
    path('contact/', views.contact, name="contact"),
    path('help/', views.help, name="help"),
    path('downloads/', views.downloads, name="downloads"),
    path('archive/', views.archive, name="archive"),
    path('statistics/', views.statistics, name="statistics"),
    path('download/<str:type>&<str:filename>', views.download, name="download"),
    path('linkdetails/<str:idA>&<str:idB>&<str:goldstandard>', views.linkDetails, name="linkDetails"),
    path('enrichment/<str:queryGenes>&<str:species>', views.enrichment, name="enrichment"),
    path('network/<str:geneQuery>&<str:genomeSelect>&<str:confidenceThreshold>&<str:grgLLRThreshold>&<str:depth>&<str:nodesPerStepOptions>&<str:expansionAlgorithm>&<str:prioritizeNeighbors>&<str:comparativeGenomes>&<str:individualEvidenceOnly>&<str:orthologsOnly>&<str:categoryID>&<str:constrainEvidence>&<str:evidenceSpecies>&<str:evidenceDataTypes>&<str:restriction>&<str:restrictPathway>&<str:restrictTissue>&<str:showAdvanced>&<str:missingMapping>&<str:multiMapping>/', views.network, name="network"),
    path('api/', views.api, name="api"),
    # path('api/geneQuery/<str:taxid>&<str:query>&<str:identifiersource>', views.apiGene, name="apiGene"),
    # path('api/geneQuery/<str:taxid>&<str:query>', views.apiGene, name="apiGene"),
    # path('api/species', views.apiSpecies, name="apiSpecies"),
    path('api/tsv/geneQuery/<str:taxid>&<str:query>/', views.apiGeneTSV, name="apiGeneTSV"),
    path('api/tsv/species/', views.apiSpeciesTSV, name="apiSpeciesTSV"),
    path('api/json/species/', views.apiSpeciesJSON, name="apiSpeciesJSON"),
    path('api/json/gene/<str:geneQuery>&<str:genomeSelect>/', views.apiGeneJSON, name="apiGeneJSON"),
    path('api/json/network/<str:geneQuery>&<str:genomeSelect>&<str:confidenceThreshold>&<str:grgLLRThreshold>&<str:depth>&<str:nodesPerStepOptions>&<str:expansionAlgorithm>&<str:prioritizeNeighbors>&<str:comparativeGenomes>&<str:individualEvidenceOnly>&<str:orthologsOnly>/', views.apiNetworkJSON, name="apiNetworkJSON"),
    path('api/json/linkdetails/<str:idA>&<str:idB>&<str:goldstandard>/', views.linkDetailsJSON, name="linkDetailsJSON")
]
