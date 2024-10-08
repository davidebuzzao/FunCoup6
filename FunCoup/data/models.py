from django.db import models
from postgres_copy import CopyManager

class Species(models.Model):
    tax_id = models.CharField(max_length=20)
    species_name = models.CharField(max_length=200)
    common_name = models.CharField(max_length=200)
    class Meta:
        indexes = [models.Index(fields=['tax_id']),
                    models.Index(fields=['species_name']),
                    models.Index(fields=['common_name'])]

class Proteome(models.Model):
    species = models.ForeignKey(Species, on_delete=models.CASCADE)
    version = models.CharField(max_length=100)
    class Meta:
        indexes = [models.Index(fields=['species','version'])]

class Protein(models.Model):
    uniprot_id = models.CharField(max_length=20)
    description = models.TextField()
    proteome = models.ForeignKey(Proteome, on_delete=models.CASCADE)
    class Meta:
        indexes = [models.Index(fields=['uniprot_id']),
                   models.Index(fields=['proteome'])]

class IdMapping(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    type = models.TextField()
    mapped_id = models.TextField()
    objects = CopyManager() 
    class Meta:
        indexes = [models.Index(fields=['protein'])]

class ProteinOrtholog(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_ortholog')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_ortholog')
    speciesA = models.CharField(max_length=20)
    speciesB = models.CharField(max_length=20)
    version = models.CharField(max_length=20) #InParanoid version
    objects = CopyManager() 
    class Meta:
        indexes = [models.Index(fields=['speciesA','speciesB'])]

class Pathway(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='protein_pathway')
    species = models.CharField(max_length=20)
    pathway_id=models.CharField(max_length=50)
    pathway_name=models.TextField()
    pathway_db=models.TextField()
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['protein']),
                   models.Index(fields=['species'])]

class Tissue(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='protein_tissue')
    species = models.CharField(max_length=20)
    tissue= models.CharField(max_length=50)
    tissue_id= models.IntegerField(default=1)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['protein']),
                   models.Index(fields=['species'])]

class GoldStandard(models.Model):
    type = models.CharField(max_length=20)
    species = models.ForeignKey(Species, on_delete=models.CASCADE)
    version = models.CharField(max_length=100)
    class Meta:
        indexes = [models.Index(fields=['type', 'version', 'species'])]

class GoldStandardLink(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_goldStandard')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_goldStandard')
    goldStandard = models.ForeignKey(GoldStandard, on_delete=models.CASCADE)
    direction = models.IntegerField(default=1) # 1:no direction 2: A->B; 3:A<-B; 5:A<->B
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['goldStandard'])]

class Evidence(models.Model):
    type = models.CharField(max_length=20)
    species = models.ForeignKey(Species, on_delete=models.CASCADE)
    version = models.CharField(max_length=100)
    scoreRange = models.TextField(null=True) #original min,max
    scoringMethod = models.CharField(max_length=200)
    class Meta:
        indexes = [models.Index(fields=['type', 'species', 'version', 'scoringMethod'])]

class DOM(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_DOM')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_DOM')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class MEX(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_MEX')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_MEX')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    sign = models.TextField(null=True)
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class MIR(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_MIR')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_MIR')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class GIN(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_GIN')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_GIN')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class GRG(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_GRG')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_GRG')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    direction = models.IntegerField(default=1) # 1:no direction 2: A->B; 3:A<-B; 5:A<->B
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class PEX(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_PEX')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_PEX')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class PHP(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_PHP')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_PHP')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class PIN(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_PIN')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_PIN')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class SCL(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_SCL')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_SCL')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class TFB(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_TFB')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_TFB')
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE)
    score = models.FloatField()
    metadata = models.TextField(null=True)
    objects = CopyManager()
    class Meta:
        indexes = [models.Index(fields=['evidence']),
                   models.Index(fields=['evidence','proteinA','proteinB'])]

class PolynomialLLR(models.Model):
    evidence = models.ForeignKey(Evidence, on_delete=models.CASCADE, related_name='evidence_polynomialLLR')
    goldStandard = models.ForeignKey(GoldStandard, on_delete=models.CASCADE, related_name='goldstandard_polynomialLLR')
    function = models.TextField(null=True)
    scoreRange = models.TextField(null=True) #min|max
    bw = models.TextField(null=True) # pos|neg 
    error = models.TextField(null=True)
    class Meta:
        indexes = [models.Index(fields=['evidence','goldStandard'])]

class Redundancy(models.Model):
    goldStandard = models.ForeignKey(GoldStandard, on_delete=models.CASCADE, related_name='goldstandard_redundancy')
    evidenceA = models.ForeignKey(Evidence, on_delete=models.CASCADE, related_name='evidenceA_redundancy')
    evidenceB = models.ForeignKey(Evidence, on_delete=models.CASCADE, related_name='evidenceB_redundancy')
    correlation = models.FloatField(default=0)
    objects = CopyManager() 
    # pvalue = models.FloatField(default=1)
    class Meta:
        indexes = [models.Index(fields=['evidenceA','evidenceB'])]
        
class Network(models.Model):
    proteinA = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinA_Network')
    proteinB = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='proteinB_Network')
    # species = models.ForeignKey(Species, on_delete=models.CASCADE)
    max_ppv = models.FloatField(default=0)
    ## Mirroring the website Interactions tab
    ppv = models.TextField(null=True) #Complex,Metabolic,Operon,PPI,Regulatory,Signaling
    fbs_goldstandard = models.TextField(null=True) #Complex,Metabolic,Operon,PPI,Regulatory,Signaling
    llr_evidence = models.TextField(null=True) #Complex:DOM,MEX,MIR,GIN,PEX,PHP,PIN,QMS,SCL,TFB|hsa,mmu,rno,cfa,...,osa;Metabolic:DOM,MEX,...||hsa,mmu,...;Operon:...
    GRG_direction = models.CharField(max_length=6, default='000000') # 0 undirected, 2 ->, 3 <- , 5 <->
    isGoldStandard = models.CharField(max_length=6, default='000000') # 0 no, 1 gs, 2 ->, 3 <- , 5 <->
    isPPVgoldstandard = models.CharField(max_length=6, default='000000') #0:PPV from evidence; 1: PPV from goldstandard
    objects = CopyManager() 
    class Meta:
        indexes = [models.Index(fields=['proteinA']),
                   models.Index(fields=['proteinB']),
                    models.Index(fields=['ppv'])]

class LogisticPPV(models.Model):
    goldStandard = models.ForeignKey(GoldStandard, on_delete=models.CASCADE, related_name='goldstandard_logisticPPV')
    function = models.TextField(null=True)
    error = models.TextField(null=True)
    class Meta:
        indexes = [models.Index(fields=['goldStandard'])]

class Stats(models.Model):
    instance_name = models.CharField(max_length=100)
    tax_id = models.CharField(max_length=20)
    species_name = models.CharField(max_length=200,default='A')
    common_name = models.CharField(max_length=200,default='A')
    num_links = models.IntegerField(default=0)
    num_nodes = models.IntegerField(default=0)
    genome_size = models.IntegerField(default=0)
    origin = models.CharField(max_length=20,default='A')
    class Meta:
        indexes = [models.Index(fields=['tax_id'])]
