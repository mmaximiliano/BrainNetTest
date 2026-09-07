# Guía para compartir: ABIDE + BrainNetTest 1.0.0

## Para qué sirve esta carpeta

Esta carpeta permite que un co-autor pueda:

1. entender de dónde salen los datos clínicos de ABIDE;
2. revisar las condiciones de uso y la licencia;
3. reconstruir los conectomas binarios usados en el paper;
4. verificar, con números, cómo se armó la muestra balanceada;
5. correr BrainNetTest sobre la muestra completa y la balanceada; y
6. entender qué cambió entre BrainNetTest 0.2.1 y 1.0.0.

Los datos son observacionales. El análisis muestra asociación entre la etiqueta
diagnóstica y el patrón de conectividad para este pipeline. No demuestra
causalidad, no es un clasificador diagnóstico y no valida biomarcadores.

## Contenido

- `BrainNetTest_1.0.0.tar.gz`: source package nuevo.
- `00_install_brainnettest.R`: instala dependencias y el package local.
- `01_prepare_abide.R`: baja y ejecuta, con checksum, el pipeline clínico
  exacto fijado en el commit del paper.
- `02_check_balance.R`: audita el balance y genera tablas legibles.
- `03_run_analysis.R`: corre el test global y los edge-wise tests.
- `run_all.R`: ejecuta instalación, preparación, balance, análisis y
  verificación en orden.
- `99_verify_bundle.R`: verifica que los archivos y resultados estén completos.
- `LICENSE_ABIDE.md`: términos de uso, links y aclaración de licencias.
- `LICENSE_CODE.md`: licencia MIT completa para package y scripts.
- `PACKAGE_ARCHIVE.txt`: identidad y provenance del package tarball compartido.
- `CAMBIOS_BRAINNETTEST_1.0.md`: explicación de la nueva versión y migración.
- `LINKS.md`: links estables al paper, dataset, licencia y referencias.
- `EXPECTED_BALANCE.csv`: resumen esperado para comprobar la reproducción.
- `CHECKSUMS.txt`: checksums de los archivos compartidos.
- `MANIFEST.csv`: size+MD5 usado por `99_verify_bundle.R`.

Los manifests fijan inputs y outputs determinísticos. No fijan bytes de
`brainnet_results.rds` ni `edge_effects.pdf`, porque RDS guarda tiempos de
runtime y PDF incluye metadata de creación. Esos dos archivos se verifican por
estructura, contenido estadístico y tamaño no vacío.

No se incluyen las time series descargadas ni los conectomas subject-level en
esta carpeta. Sí se incluyen outputs estadísticos agregados y una auditoría de
pairs con IDs anonimizados, todos sujetos a los términos de ABIDE. Los datos de
entrada se generan localmente para que cada persona revise y respete esos
términos.

El package tarball es un fresh `R CMD build` del commit `695f534`.
`PACKAGE_ARCHIVE.txt` identifica este rebuild exacto y explica por qué su hash
puede diferir del archive de submission que contiene timestamps distintos.

## Links principales

### Paper de BrainNetTest

Manuscrito exacto usado para esta guía, fijado al commit `695f534`:

- Vista en GitHub:
  <https://github.com/mmaximiliano/BrainNetTest/blob/695f534/paper_source/article.pdf>
- Descarga directa:
  <https://raw.githubusercontent.com/mmaximiliano/BrainNetTest/695f534/paper_source/article.pdf>
- Rama de trabajo:
  <https://github.com/mmaximiliano/BrainNetTest/tree/feat/brainnettest-jss-rebuild>

El PDF es el manuscrito de resubmission/preprint, no una publicación JSS
aceptada y todavía no tiene DOI propio.

### Dataset y papers que hay que citar

- ABIDE I, página oficial y Usage Agreement:
  <http://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html>
- ABIDE Preprocessed / PCP:
  <http://preprocessed-connectomes-project.org/abide/>
- Instrucciones de descarga PCP:
  <http://preprocessed-connectomes-project.org/abide/download.html>
- Paper de ABIDE:
  <https://doi.org/10.1038/mp.2013.78>
- Paper/abstract del Preprocessed Connectomes Project:
  <https://doi.org/10.3389/conf.fninf.2013.09.00041>
- Atlas AAL:
  <https://doi.org/10.1006/nimg.2001.0978>
- Licencia CC BY-NC-SA 3.0:
  <https://creativecommons.org/licenses/by-nc-sa/3.0/>

## Qué datos se usan

No se descargan imágenes fMRI crudas. El script usa derivados públicos de
ABIDE I preprocesados por el Preprocessed Connectomes Project:

- site: `NYU`;
- pipeline: `C-PAC`;
- strategy: `filt_noglobal`;
- derivative: `rois_aal`;
- atlas: `AAL`;
- quality control: `mean framewise displacement < 0.2`;
- connectivity: Pearson correlation entre ROI time series;
- threshold: top 10% de las correlaciones Pearson más altas de cada persona
  (valor signed, no valor absoluto); y
- red final: binaria, simétrica, sin self-loops.

El source es el bucket HTTP público de PCP:

```text
https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/
```

La primera corrida baja aproximadamente 35 MB. Los archivos quedan en
`abide_cache/` y se reutilizan en corridas posteriores.

## Licencia y condiciones de uso

La página oficial de ABIDE I dice que el uso es irrestricto para investigación
no comercial bajo los protocolos INDI/FCP y lo identifica como
Creative Commons Attribution-NonCommercial-ShareAlike.

Para trabajar de forma conservadora:

1. usar los datos solamente para investigación no comercial;
2. citar ABIDE I y PCP;
3. identificar el site y el subset usados (`NYU`);
4. conservar la misma licencia al redistribuir derivados;
5. incluir los acknowledgments/funding del site cuando corresponda; y
6. no intentar reidentificar participantes.

El bucket PCP usado por el script es público y no requiere login. PCP pide
citar la iniciativa. El acceso a los datos ABIDE crudos por NITRC requiere
registración bajo los términos oficiales ABIDE/INDI. Que el derivative
preprocesado sea accesible por HTTP no elimina las condiciones de uso del
dataset original. `LICENSE_ABIDE.md` incluye también los funding
acknowledgements oficiales para NYU.

El código de BrainNetTest usa licencia MIT. Esa licencia de software no cambia
la licencia de los datos. Ver `LICENSE_CODE.md` para el código y
`LICENSE_ABIDE.md` para los datos.

## Cómo se selecciona la muestra completa

`01_prepare_abide.R` reproduce este filtro, sin mirar ningún network outcome:

1. lee el phenotypic file de ABIDE PCP;
2. elimina filas sin `FILE_ID`;
3. conserva solamente `SITE_ID == "NYU"`;
4. exige `func_mean_fd` disponible;
5. exige `func_mean_fd < 0.2`; y
6. ordena determinísticamente por diagnóstico y `FILE_ID`.

Resultado esperado:

- control: 98 personas;
- autism: 73 personas;
- total: 171 personas.

La muestra completa no está balanceada:

- sexo: autism 64 male / 9 female; control 72 male / 26 female;
- edad media: autism 14.922; control 15.674 años;
- mean FD: autism 0.07674; control 0.05460;
- standardized mean difference (SMD) de edad: `-0.113`; y
- SMD de mean FD: `0.604`.

El desbalance de movimiento es importante. Por eso el paper informa la muestra
completa y además una sensitivity sample balanceada.

## Cómo se construye la muestra balanceada

El procedimiento es deterministic coarsened-stratified 1:1 subsampling. Usa la
misma lógica de coarsened exact matching para formar strata, pero
`balance_id` se conserva solamente para auditar la selección: no convierte el
análisis posterior en un paired design. No es propensity score matching y no
usa ninguna arista ni resultado de red.

Para cada persona se arma:

```text
balance_stratum =
  SEX
  x floor(AGE_AT_SCAN / 5)
  x floor(func_mean_fd / 0.05)
```

Interpretación:

- sexo: matching exacto;
- edad: bins de 5 años;
- movimiento: bins de mean FD de 0.05.

Dentro de cada stratum:

1. se cuentan autism y control;
2. se toma `min(n_autism, n_control)` de cada grupo;
3. se eligen los primeros IDs anonimizados en orden de `FILE_ID`; y
4. se asigna el mismo `balance_id` a cada par.

Resultado:

- 55 autism + 55 control = 110 personas;
- 49 male + 6 female en cada grupo;
- edad media: control 14.249, autism 14.746;
- SMD de edad: `0.089`;
- mean FD: control 0.061748, autism 0.061826;
- SMD de mean FD: `0.0026`;
- cada `balance_id` tiene exactamente un control y un participante autism; y
- cada par comparte el mismo `balance_stratum`.

Retención:

- autism: 55/73 = 75.3%;
- control: 55/98 = 56.1%;
- total: 110/171 = 64.3%.

### Cómo entra esta muestra al test

`03_run_analysis.R` trata las 55 redes control y 55 autism como dos grupos de
observaciones independientes. El randomization test mezcla etiquetas entre las
110 personas manteniendo 55/55; no restringe las permutaciones dentro de
`balance_id`.

Los IDs de par sirven para verificar que el subsampling balanceó los strata,
no como blocking variable inferencial. BrainNetTest 1.0.0 todavía no implementa
restricted/block permutations. Si se quisiera una inferencia condicionada a
pairs, haría falta otro procedimiento y nuevos resultados.

### Qué balancea y qué no

Balancea exactamente sexo y de forma coarsened edad y movimiento. No balancea
IQ, medication, comorbidities, handedness, socioeconomic variables ni
confounders no medidos. Dentro de un mismo bin pueden quedar pequeñas
diferencias continuas. Por eso se llama sensitivity analysis y no prueba que
desapareció todo confounding.

`02_check_balance.R` vuelve a calcular estos números, verifica cada par y frena
con error si no coinciden las reglas.

## Balance de densidad de red

Hay un segundo balance, distinto al balance clínico: cada conectoma conserva
las 644 correlaciones más fuertes de 6.441 posibles aristas, aproximadamente
el 10%.

Así, las 171 redes tienen exactamente 644 aristas. Esto evita que el test se
limite a detectar que un grupo tiene redes globalmente más densas. El análisis
pregunta dónde aparecen las aristas, condicionado al threshold elegido.

Esta decisión también es una limitación: un threshold diferente puede cambiar
el resultado.

## Ejecución paso a paso

### Requisitos

- R 4.0 o superior;
- internet para la primera descarga;
- aproximadamente 35 MB para downloads, más espacio para cache/output;
- toolchain de compilación de R si `igraph` no tiene binary para el sistema.

### Opción rápida: todo junto

Desde esta carpeta:

```sh
Rscript run_all.R
```

### Opción recomendada: revisar cada paso

```sh
Rscript 00_install_brainnettest.R
Rscript 01_prepare_abide.R
Rscript 02_check_balance.R
Rscript 03_run_analysis.R
Rscript 99_verify_bundle.R
```

Outputs:

- `application/abide_connectomes.rds`;
- `application/abide_subjects.csv`;
- `application/abide_provenance.csv`;
- `application/input_manifest.csv` con size+MD5 de pipeline, phenotypic file y
  ROI time series;
- `results/balance_summary.csv`;
- `results/balance_diagnostics.txt`;
- `results/brainnet_global_results.csv`;
- `results/brainnet_balanced_edges.csv`;
- `results/brainnet_top_edges.csv`; y
- `results/brainnet_results.rds`.

## Resultados esperados del paper

Con 999 random assignments y seed 42:

- full cohort: `T` aproximadamente `-44.7`, `p = 0.001`;
- balanced cohort: `T` aproximadamente `-27.4`, `p = 0.01`;
- edge hypotheses en la muestra balanceada: 6.441;
- raw p-values `< 0.001`: 3;
- Holm-selected edges: 0; y
- BY-selected edges: 0.

Interpretación correcta: hay evidencia global de una diferencia distribuida
en el patrón de aristas para este site/pipeline/threshold, pero no hay evidencia
suficiente para afirmar una conexión puntual después de multiplicity control.

## Qué cambia en BrainNetTest 1.0.0

La versión nueva es breaking:

- ahora hay que construir `brainnet_data`;
- `brainnet_test()` reemplaza a `identify_critical_links()`;
- el resultado es un objeto S3 `brainnet_result`;
- `selected_edges()` reemplaza `result$critical_edges`;
- `selected_nodes()` reemplaza `get_critical_nodes()`;
- `plot(result)` reemplaza `plot_critical_edges(populations, result)`;
- el test global usa `|T|`, inclusive ties y Monte Carlo plus-one correction;
- cada arista usa Fisher two-sided y Holm por default; y
- la ablación adaptativa queda explícitamente como descriptive diagnostic.

Antes:

```r
result <- identify_critical_links(populations)
result$critical_edges
```

Ahora:

```r
networks <- brainnet_data(populations, node_labels = node_labels)
result <- brainnet_test(
  networks,
  n_permutations = 999,
  adjust = "holm",
  seed = 42
)
selected_edges(result)
```

La explicación completa está en `CAMBIOS_BRAINNETTEST_1.0.md`.

## Qué conviene decir y qué conviene evitar

Sí:

- "La etiqueta diagnóstica está asociada con un patrón global distribuido para
  este pipeline."
- "El resultado persiste en una sensitivity sample coarsened-balanced."
- "Ninguna arista individual sobrevive Holm o BY."

No:

- "El autismo causa este patrón."
- "Encontramos biomarcadores."
- "La muestra balanceada elimina todo confounding."
- "El 10% es el único threshold correcto."
- "Un p global significativo prueba una arista específica."

## Troubleshooting

- Si falla la descarga, volver a correr: el cache evita bajar de nuevo los
  archivos completos.
- Para cambiar el cache:

  ```sh
  ABIDE_CACHE=/ruta/al/cache Rscript 01_prepare_abide.R
  ```

- Si `01_prepare_abide.R` informa checksum incorrecto, no ejecutar el script
  remoto: verificar el commit/link.
- Si cambia el número de participantes, revisar disponibilidad del bucket,
  phenotypic file, `mean FD`, site y versionado del script.
- Si cambia el balance, revisar que el input esté ordenado por `FILE_ID` y que
  los bins sigan siendo 5 años y 0.05 mean FD.

