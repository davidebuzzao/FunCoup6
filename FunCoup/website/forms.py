from django import forms

# Here is the FunCoup query form defined
# All fields necessary to make a query are specified here, together with all static options in dropdowns as well as styling of the components
class Query(forms.Form):
    geneQuery = forms.CharField(widget=forms.Textarea(attrs={"rows":"1","id":"geneIdSearch","placeholder":"Enter one or multiple identifiers separated by space",'class': "form-control form-control-gene","oninput": "auto_grow(this)"}),required=True, strip=True, help_text="one or multiple space-separated identifiers")
    genomeSelect = forms.CharField(widget=forms.TextInput(attrs={"id":"genomeSelect","autocomplete":"off","placeholder":"Enter a species name",'class': 'form-control form-control-gene autocompleteInput'}),required=True, strip=True)

    showAdvanced = forms.CharField(widget = forms.HiddenInput(attrs={"id":"showAdvancedFlag"}), required = False, initial="false")
    confidenceThreshold = forms.CharField(widget=forms.TextInput(attrs={'type': 'range',"id":"confidenceThreshold","class":"form-range","min":0.50, "max":1, "step":0.01, "oninput":"this.nextElementSibling.value = this.value"}),initial=0.95)
    grgLLRThreshold = forms.CharField(widget=forms.TextInput(attrs={'type': 'range',"id":"grgLLRThreshold","class":"form-range","min":0, "max":2, "step":0.25, "oninput":"this.nextElementSibling.value = this.value"}),initial=1)
    depth = forms.CharField(widget=forms.TextInput(attrs={'type': 'range',"id":"depth","class":"form-range","min":0, "max":3, "step":1, "oninput":"this.nextElementSibling.value = this.value"}),initial=1)
    nodesPerStepOptions = forms.CharField(widget=forms.TextInput(attrs={'type': 'range',"id":"nodesPerStepOptions","class":"form-range","min":1, "max":50, "step":1, "oninput":"this.nextElementSibling.value = this.value"}),initial=15)
    
    expansionAlgorithmOptions=(('maxlink', 'MaxLink'),('topas', 'TOPAS'),('local', 'Each gene independently'), ('group', 'Genes as group'))
    expansionAlgorithm = forms.ChoiceField(widget=forms.RadioSelect(attrs={"id":"expansionAlgorithm"} ),choices=expansionAlgorithmOptions,initial=['group'])
    prioritizeNeighborsOptions=(('True', 'Prioritize common neighbors'))
    prioritizeNeighbors = forms.BooleanField(widget=forms.CheckboxInput(attrs={"style":"margin-left: 15px;"}),initial=True, required = False)

    comparativeGenomes = forms.ChoiceField(widget=forms.CheckboxSelectMultiple(attrs={"class":"form-check-input comparativeGenomesSelect"}),choices=[])
    individualEvidenceOnly = forms.BooleanField(widget=forms.CheckboxInput({"class":"form-check-input","type":"checkbox"}),initial=True, required = False)
    orthologsOnly = forms.BooleanField(widget=forms.CheckboxInput({"class":"form-check-input","type":"checkbox"}),initial=True, required = False)

    categoryID = forms.ChoiceField(widget=forms.RadioSelect(attrs={}),choices=[])
    constrainEvidenceOptions=(('all', 'Use all evidence'), ('species', 'Use evidence from species:'), ('evidence', 'Use evidence from data types:'))
    constrainEvidence = forms.ChoiceField(widget=forms.RadioSelect(attrs={"style":"margin-top: 1px;"} ),choices=constrainEvidenceOptions)
    evidenceSpecies = forms.ChoiceField(widget=forms.CheckboxSelectMultiple(attrs={"class":"form-check-input comparativeGenomesSelect","type":"checkbox","style":"margin-left:10px;"}),choices=[])
    evidenceDataTypes = forms.ChoiceField(widget=forms.CheckboxSelectMultiple(attrs={"class":"form-check-input comparativeGenomesSelect","type":"checkbox","style":"margin-left:10px;"}),choices=[])

    restrictionOptions=(('all', 'All genes (no restriction)'), ('pathway', 'Genes with the pathway:'),('tissue', 'Genes with the tissue:'))
    restriction = forms.ChoiceField(widget=forms.RadioSelect(attrs={"style":"margin-top: 1px;"} ),choices=restrictionOptions,initial=['all'])
    restrictPathway = forms.CharField(widget=forms.TextInput(attrs={"id":"restrictPathway","autocomplete":"off","placeholder":"Start typing a pathway name",'class': 'form-control form-control-gene autocompleteInput'}),required=False, strip=True)
    restrictTissue = forms.CharField(widget=forms.TextInput(attrs={"id":"restrictTissue","autocomplete":"off","placeholder":"Start typing a tissue name",'class': 'form-control form-control-gene autocompleteInput'}),required=False, strip=True)
