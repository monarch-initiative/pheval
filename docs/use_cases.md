# Use Case - PhEval

## **Id:** UC-01

### **Use Case:** PhEval - Exomiser Scrambling process

### **Description:** PhEval runs the exomiser with the default configurations first. Next, PhEval scrambles the phenotypic mapping tables from the exomiser using a scalar factor. Scrambled values are stored in the phenotypic database. Finally, the exomiser runs using the updated scrambled values.

```mermaid
 sequenceDiagram
		autonumber
    PhEval->>Exomiser: run_exomiser()
    loop dump of mapping tables
        PhEval->>PhenotypeDB: dump_table(table_name)
    end
    loop Scramble (SIMJ,IC,SCORE) columns
        PhenotypeDB->>Scramble: scramble(column, scramble_factor)
    end
		
		PhEval-->Scramble: scrambled values
		PhEval->>PhenotypeDB: insert_scrambled_values()

		PhEval->>Exomiser: run_exomiser()
```
