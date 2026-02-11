# Plugins

This page lists **available PhEval plugins** and the phenotype-driven tools they integrate.

Plugins are separate Python packages that extend PhEval with tool-specific logic. Each plugin exposes one or more **runners** that can be invoked via the PhEval CLI.

For details on how plugins and runners fit into the execution model, see:

- [Using PhEval → Plugins and runners](../using_pheval/plugins_and_runners.md)

For tool-specific usage, configuration options, and examples, always refer to the **plugin README**.

---

## Available plugins

The table below lists currently implemented PhEval plugins, with links to both the PhEval runner and the original tool where applicable.

| Tool        | PhEval plugin | Original tool | Notes |
|-------------|---------------|---------------|-|
| Exomiser    | [pheval.exomiser](https://github.com/monarch-initiative/pheval.exomiser) | [Exomiser](https://github.com/exomiser/Exomiser) | |
| Phen2Gene   | [pheval.phen2gene](https://github.com/monarch-initiative/pheval.phen2gene) | [Phen2Gene](https://github.com/WGLab/Phen2Gene) | |
| LIRICAL     | [pheval.lirical](https://github.com/monarch-initiative/pheval.lirical) | [LIRICAL](https://github.com/TheJacksonLaboratory/LIRICAL) | |
| SvAnna      | [pheval.svanna](https://github.com/monarch-initiative/pheval.svanna) | [SvAnna](https://github.com/TheJacksonLaboratory/SvAnna) | |
| GADO        | [pheval.gado](https://github.com/monarch-initiative/pheval.gado) | [GADO](https://github.com/molgenis/systemsgenetics/wiki/GADO-Command-line) | |
| PhenoGenius | [pheval.phenogenius](https://github.com/monarch-initiative/pheval.phenogenius) | [PhenoGenius](https://github.com/kyauy/PhenoGeniusCli) | |
| AI MARRVEL  | [pheval.ai_marrvel](https://github.com/monarch-initiative/pheval.ai_marrvel.git) | [AI-MARRVEL](https://github.com/LiuzLab/AI_MARRVEL.git) | |
| OAK         | [pheval.oak](https://github.com/monarch-initiative/pheval.oak.git) | | |
| OntoGPT     | [pheval.ontogpt](https://github.com/monarch-initiative/pheval.ontogpt) | | |
| MALCO       | [pheval.llm](https://github.com/monarch-initiative/pheval.llm.git) | | |
| ELDER       | [ELDER](https://github.com/monarch-initiative/ELDER) | | |
| Template    | [pheval.template](https://github.com/monarch-initiative/pheval.template) | | Starter template |

---

## Notes

- Plugins may expose **multiple runners**, depending on supported modes or input types.
- Plugin availability does not imply identical functionality across tools.
- Experimental or research plugins may change more frequently than core plugins.

Users should always consult the plugin repository for the most up-to-date information.

