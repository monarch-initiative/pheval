# Home

## Introduction

PhEval - Phenotypic Inference Evaluation Framework

### PhEval: Tool-specific processing (VP pipeline)

```mermaid
flowchart LR
    PC-->DP
    PC[(Phenopackets Corpus)]
    SSSOM[Semantic Similarity Profiles Mapping Commons]-->|OAK-SEMSIM|DP[Data Prepare]
    KG[Source data KG - Monarch KG]-->|KGX-BIOLINK|DP[Data Prepare]
    ONT[Ontologies - Phenio]-->|OAK-ONTO|DP[Data Prepare]
    DP-->RP[Run Prepare]
    RP-->PR[PhEval Runner]
    PR-->DP2[Data Process]
    ER[Exomiser Runner]-->PR
    EDP[Exomiser Data Prepare]-->DP
    ERP[Exomiser Run Prepare]-->RP
    PPP[Disease-profile similarity prediction Post-process]-->DP2
    PV[Phenotype/Variant]-->DP2
    GVP[Gene VP Post-process]-->DP2
    EPP[Exomiser Post Process]-->GVP
    GVP-->VPR[VP Report]
```

**Quick links:**

- [GitHub page](https://github.com/monarch-initiative/pheval/)
