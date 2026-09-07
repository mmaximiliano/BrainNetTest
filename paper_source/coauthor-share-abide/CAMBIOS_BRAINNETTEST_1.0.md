# BrainNetTest 1.0.0: qué cambia

## Resumen

La versión 1.0.0 es una re-fundación del package. El objetivo principal es
separar tres cosas que antes estaban mezcladas:

1. global inference;
2. edge-wise inference con multiplicity control; y
3. exploratory ablation.

## API nueva

### Input validado

```r
networks <- brainnet_data(
  list(control = control_graphs, autism = autism_graphs),
  node_labels = aal_labels
)
```

`brainnet_data()` rechaza:

- matrices no cuadradas;
- valores distintos de 0/1;
- asimetría;
- diagonal no cero;
- distinto número u orden de nodos;
- missing/non-finite values; y
- grupos con menos de dos redes.

### Análisis

```r
result <- brainnet_test(
  networks,
  n_permutations = 999,
  edge_method = "fisher",
  adjust = "holm",
  seed = 42
)
```

### Output S3

```r
print(result)
summary(result)
as.data.frame(result)
selected_edges(result)
selected_nodes(result)
plot(result)
```

El resultado es `brainnet_result` y guarda el test global, todos los edge
tests, central graphs, parámetros, finite p-value resolution y metadata.

## Cambios estadísticos

### Global test

- cada undirected edge se cuenta una sola vez, usando upper triangle;
- se usa two-sided extremeness `abs(T)`;
- exact enumeration incluye el observed assignment;
- Monte Carlo sampling usa el complete fixed-size orbit;
- el p-value usa inclusive ties; y
- Monte Carlo usa `(1 + extremos) / (B + 1)`.

### Edge-wise inference

- se testean todas las posibles aristas;
- default: Fisher exact two-sided;
- default adjustment: Holm;
- BY está disponible para FDR bajo general dependence;
- BH es opt-in y requiere assumptions adicionales; y
- una arista se llama `selected` solamente por su adjusted p-value.

El global p-value y los edge p-values responden preguntas distintas.

### Ablation

El método viejo ordenaba aristas y paraba cuando el test dejaba de rechazar.
Esa parada usaba los mismos datos para seleccionar y evaluar, por lo que no se
presenta más como valid post-selection inference.

Se conserva:

```r
result <- brainnet_test(networks, ablation = TRUE, seed = 42)
ablation_path(result)
```

`descriptive_tail_fraction` es descriptivo; no es un p-value.

## Migración rápida

| Antes (0.2.1) | Ahora (1.0.0) |
|---|---|
| raw nested list | `brainnet_data()` |
| `identify_critical_links()` | `brainnet_test()` |
| `result$critical_edges` | `selected_edges(result)` |
| `get_critical_nodes()` | `selected_nodes(result)` |
| `plot_critical_edges()` | `plot(result)` |
| `compute_central_graph()` | `central_graph()` |
| `compute_distance()` | `graph_distance()` |

## Helpers y argumentos removidos

Estos helpers dejaron de ser public API:

- `compute_test_statistic()`;
- `compute_edge_frequencies()`;
- `compute_edge_pvalues()`; y
- `rank_edges()`.

La información equivalente se obtiene desde `brainnet_test()` y
`as.data.frame(result)`, que conserva proportions, effects, raw p-values,
adjusted p-values y decisions en un objeto validado.

También se removieron:

- argumento `a`: un positive scale factor no cambia el randomization order;
- argumento `adjust_method`: antes ajustaba valores usados para ordenar, pero
  no daba una regla válida de selección; ahora `adjust = "holm"` controla la
  decisión por arista;
- `modified_populations`: el análisis ya no devuelve copias de redes con
  aristas borradas; y
- el stopping set adaptativo del resultado principal.

## RNG behavior

Cuando se pasa `seed`, la versión 1.0 restaura el RNG state del caller al
terminar, incluso si ocurre un error. Esto evita que una llamada interna cambie
silenciosamente la secuencia aleatoria del análisis que la rodea.

## Consecuencia para ABIDE

El análisis ABIDE muestra:

- global association en la muestra completa;
- global association en la sensitivity sample balanceada;
- ningún edge individual seleccionado por Holm/BY.

La versión nueva evita llamar "critical" a los tres raw-smallest p-values. Se
informan como contexto anatómico, no como biomarcadores.
