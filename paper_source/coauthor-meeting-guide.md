# Paper BrainNetTest: guía para la reunión con el co-autor

## Tres puntos para acordar primero

1. El test global y los tests por arista responden preguntas distintas.
2. Un resultado global no prueba que una arista puntual sea diferente.
3. El resultado de ABIDE es una asociación para un solo pipeline de datos. No
   es causal y no es un resultado diagnóstico.

## Cómo se mueven los datos dentro del paquete

1. Cada persona tiene una matriz de adyacencia cuadrada.
2. `1` significa que la arista (conexión) está presente; `0` que está ausente.
3. La matriz debe ser simétrica, tener diagonal cero y usar el mismo orden de
   nodos para todas las personas.
4. `brainnet_data()` revisa estas reglas. No modifica las matrices de entrada.
5. El paquete lee solo el triángulo superior, así cada arista no dirigida se
   cuenta una sola vez.
6. Para cada grupo y arista, cuenta cuántas personas tienen esa arista.
7. El test global usa todas las frecuencias de aristas juntas.
8. El paso de aleatorización mueve personas enteras entre las etiquetas de
   grupo, manteniendo fijos los tamaños de grupo. **No** mezcla nodos ni
   aristas.
9. Cada arista además recibe su propio test de Fisher. La corrección de Holm
   tiene en cuenta la gran cantidad de aristas testeadas.
10. Los resultados se guardan en un objeto `brainnet_result`. Las redes
    originales no se editan.

## Secciones del paper

### Resumen (Abstract)

**Qué dice, simple:** Presenta el problema, el nuevo diseño del software, las
dos capas estadísticas, los resultados principales de simulación y el hallazgo
en ABIDE.

**Por qué está:** Quien lee tiene que entender el paquete y sus límites sin
leer todo el paper.

### 1. Introducción

**Qué dice, simple:** Un estudio clínico puede preguntar:

- ¿Los grupos difieren en toda la red?
- ¿Qué aristas puntuales difieren después de corregir por muchos tests?

También explica en qué se diferencia BrainNetTest de NBS, NBR, FIAR,
brainGraph, NetworkComparisonTest e igraph.

**Por qué está:** Muestra el hueco de software y hace que el paper sea
software-first, como pidió JSS.

### 2. Contrato de análisis y contexto clínico

**Qué dice, simple:** El paquete soporta grupos independientes de redes
alineadas, binarias y no dirigidas.

El paquete no arma las redes a partir de imágenes crudas. La elección del
atlas, el control de movimiento, la estimación de correlaciones, el umbral y
el control de calidad ocurren antes de BrainNetTest.

**Por qué está:** Le dice al usuario cuándo el método es válido y evita
afirmaciones fuera del tipo de dato soportado.

### 3. Diseño del software e interfaz de R

**Qué dice, simple:**

- `brainnet_data()` valida la entrada.
- `brainnet_test()` corre el análisis.
- `print()`, `summary()`, `plot()` y `as.data.frame()` dan el comportamiento
  normal de R.
- `selected_edges()` y `selected_nodes()` dan vistas enfocadas.

**Por qué está:** El paquete viejo devolvía listas simples. Las clases nuevas
hacen que los errores se detecten antes y que los resultados sean más fáciles
de usar.

### 4. Motor estadístico y computacional

#### Grafo central y score global

**Qué dice, simple:** Un grafo central es la matriz de adyacencia promedio de
un grupo. Un valor de `0.70` significa que el 70% del grupo tiene esa arista.
No es el grafo de una persona.

El score global combina todas las diferencias de frecuencia de aristas entre
grupos en un solo número.

**Por qué está:** Da un único test para una diferencia distribuida en la red.

#### P-valor por aleatorización

**Qué dice, simple:** El paquete mueve muchas veces redes enteras entre las
etiquetas y recalcula el score. Compara `|T|`, así cualquiera de los dos lados
puede ser extremo.

Los espacios chicos de asignaciones se listan por completo. Los espacios
grandes usan asignaciones al azar y la regla `(extremos + 1) / (B + 1)`.

**Por qué está:** Da un p-valor válido en muestra finita y evita p-valores
Monte Carlo iguales a cero.

#### Tests por arista

**Qué dice, simple:** Cada arista arma una tabla de presente/ausente y usa un
test de Fisher a dos colas. Holm es la corrección por defecto.

**Por qué está:** Testear miles de aristas sin corrección generaría muchos
falsos hallazgos.

#### Ablación

**Qué dice, simple:** Las aristas se ordenan y sus contribuciones al score se
van sacando con sumas de prefijos rápidas.

**Por qué está:** Conserva la mejora de velocidad útil del paquete viejo.

**Límite importante:** La fracción de cola de la ablación es descriptiva. No es
un p-valor post-selección y no define un conjunto "crítico" probado.

### 5. Flujos de trabajo del software

**Qué dice, simple:** El paper muestra el mismo flujo para dos y tres grupos.
El gráfico usa datos que ya están guardados en el resultado.

**Por qué está:** Se ve la API pensada sin listados de código largos.

### 6. Aplicación clínica: datos de autismo ABIDE

**Preparación de los datos:**

- Un solo sitio: NYU.
- Filtro de movimiento: desplazamiento medio (framewise) menor a `0.2`.
- Pipeline C-PAC, datos filtrados, sin regresión de señal global.
- Series temporales del atlas AAL.
- Correlaciones de Pearson entre regiones.
- Se conserva el 10% de aristas más fuertes para cada persona.
- Se quitan las regiones sin señal en alguna persona, quedan 114 regiones.

La muestra completa tiene 98 controles y 73 participantes con autismo. Una
segunda muestra balanceada, determinística, tiene 55 personas por grupo,
emparejadas por estrato de sexo, bins de edad de cinco años y bins de
movimiento. Los resultados de la red no se usaron para elegir la muestra
emparejada.

**Resultado principal:**

- Muestra completa: `T = -44.7`, `p = 0.001`.
- Muestra balanceada: `T = -27.4`, `p = 0.01`.
- Análisis balanceado: 6.441 tests por arista.
- Tres aristas tienen `p` crudo `< 0.001`.
- Ninguna arista sobrevive a la corrección de Holm ni de BY.

**Cómo leer la figura de la Sección 6:**

#### Panel A: mapa de calor de diferencia de frecuencia de aristas

- Cada fila y columna es una de las 114 regiones del atlas AAL.
- Las regiones están ordenadas por sistemas anatómicos amplios.
- Cada celda es la frecuencia de la arista en el grupo con autismo menos la
  frecuencia en el grupo control.
- Rojo: la arista es más común en el grupo con autismo.
- Azul: la arista es menos común en el grupo con autismo.
- Un color cercano al blanco: poca diferencia observada.
- Los bloques ayudan a ver si las diferencias están dentro o entre sistemas
  cerebrales.

Este panel muestra el patrón anatómico de los efectos observados. No muestra
qué aristas pasan una corrección estadística.

#### Panel B: efecto por arista contra evidencia cruda

- Cada punto es una de las 6.441 aristas testeadas.
- El eje horizontal es la diferencia de frecuencia autismo-menos-control.
- Los puntos a la derecha son más comunes en autismo; los de la izquierda,
  menos comunes en autismo.
- El eje vertical es `-log10(p-valor crudo)`. Un punto más alto tiene un
  p-valor más chico.
- Los diez puntos resaltados tienen los p-valores crudos más chicos.

Los puntos resaltados no son biomarcadores confirmados. Se resaltan solo para
mostrar los resultados crudos más fuertes. Ninguno pasa Holm ni BY.

**Por qué se usan los dos paneles:** El Panel A muestra dónde aparecen las
diferencias en el cerebro. El Panel B muestra la dirección del efecto y la
evidencia estadística cruda. Apoyan el resultado global de patrón distribuido,
pero no una afirmación sobre una arista puntual.

**Qué significa, simple:** Los grupos difieren en el patrón distribuido de
aristas, pero los datos no apoyan una afirmación sobre una arista específica.

**Por qué está:** Muestra un uso real y también por qué el resultado global y
el de aristas deben mantenerse separados.

### 7. Validación por simulación

**Tests de nulo completo:** Cuando los grupos se generan del mismo modelo, las
tasas de falsos positivos quedan dentro de los límites planeados.

- Mayor tasa de rechazo global bajo el nulo: `0.051`.
- Mayor tasa de error por familia de Holm: `0.015`.

**Tests de nulo parcial:** Holm es estricto.

- Tasa media de verdaderos positivos de Holm: `0.191`.
- Precisión media de Holm: `0.716`.

**Qué significa, simple:** Holm se pierde muchas aristas verdaderas en muestras
chicas, pero la mayoría de las que selecciona son verdaderas.

**Por qué está:** Chequea tanto el control de error como el costo de ser
conservador.

### 8. Rendimiento computacional

**Motor global:** Los cálculos matriciales por bloques son hasta cerca de
`19.1x` más rápidos en los casos medidos.

**Actualización de ablación:** Las sumas de prefijos son hasta cerca de `953x`
más rápidas para el paso de actualización aislado.

**Límite importante:** `953x` no es una aceleración del paquete de punta a
punta. El ranking y el cálculo de contribuciones quedan fuera de ese tiempo.

**Por qué está:** Prueba igualdad numérica y velocidad medida, sin usar
estimaciones de varios días.

### 9. Discusión

**Qué dice, simple:** El paquete ahora tiene un flujo clínico claro,
conclusiones globales y por arista válidas, y optimizaciones exactas de
software.

Límites que quedan:

- solo redes binarias y no dirigidas;
- solo personas independientes;
- sin covariables;
- sin diseño de familia, sitio, pareado o medidas repetidas;
- sin afirmación causal ni diagnóstica;
- los resultados dependen del preprocesamiento, el atlas y el umbral.

**Por qué está:** Deja claro qué prueba el paper y qué no prueba.

### Detalles computacionales y replicación

**Qué dice, simple:** Lista las dependencias de R y los scripts que
reconstruyen las simulaciones, los benchmarks, los conectomas de ABIDE, el
manuscrito y la información de sesión.

**Por qué está:** Quienes revisan en JSS tienen que poder reproducir el
trabajo.

## ¿El paquete todavía identifica las aristas críticas?

**Respuesta corta:** Sí, pero cambió cómo se hace. La identificación válida de
aristas ahora se hace con los tests por arista con control de multiplicidad
(Holm por defecto), y se expone con `selected_edges()` (y `selected_nodes()`
para el resumen por nodo).

**Qué pasó con el método viejo:** El `identify_critical_links()` original
sacaba aristas de forma adaptativa y devolvía un "conjunto crítico" con una
regla de parada. Eso se quitó como inferencia, porque usaba los mismos datos
para ordenar las aristas y para decidir dónde parar; no era un p-valor válido
post-selección. Ese nombre ya no se exporta.

**Qué se conservó:** El cálculo rápido por prefijos sigue vivo, pero solo como
diagnóstico descriptivo, en `ablation_path()`. Devuelve una columna
`descriptive_tail_fraction`, que no es un p-valor, y no devuelve un conjunto de
parada.

**Cómo se identifican las aristas ahora (pasos):**

1. `brainnet_test()` corre un test de Fisher por cada arista.
2. Aplica Holm por defecto para controlar los falsos positivos.
3. Una arista queda "seleccionada" si su p-valor ajustado es `<= alpha`.
4. `selected_edges()` devuelve esas aristas; `selected_nodes()` resume los
   nodos que participan.

**Dónde está esto en el paper:**

- Sección 3: lista `selected_edges()` y `selected_nodes()` como la interfaz.
- Sección 4.2 (Inferencia por arista con control de multiplicidad): define el
  test de Fisher por arista y las correcciones Holm/BY/BH.
- Sección 4.3 (Cómputo separable por arista y ablación): describe
  `ablation_path()` y aclara su límite (es descriptivo, no inferencia).
- Sección 6 (ABIDE): es la identificación real sobre datos reales. Tres aristas
  con `p` crudo `< 0.001`, pero ninguna sobrevive Holm ni BY. Se ve en la
  figura, Panel B.
- Sección 7: mide qué tan bien identifica aristas verdaderas (tasa de
  verdaderos positivos y precisión de Holm y BY en las simulaciones).

**Punto clave para la reunión:** "Identificar aristas críticas" sigue siendo
posible, pero ahora es honesto: solo se marcan aristas que pasan una corrección
estricta. En ABIDE ninguna pasa, así que no afirmamos aristas individuales; sí
afirmamos una diferencia global distribuida.

## Tests y términos, en palabras simples

### Test exacto de Fisher

Para una arista, armamos una tabla chica: cuántas personas tienen la arista y
cuántas no, en cada grupo. El test de Fisher pregunta qué tan probable es ese
reparto si los grupos tuvieran la misma probabilidad de arista. Un p-valor
chico significa que la frecuencia de la arista se ve distinta entre grupos. Se
llama "exacto" porque usa conteos exactos, no una aproximación de muestra
grande, así que es seguro para aristas raras.

### Test de aleatorización (permutación)

Es el test global. Dejamos cada red como está, pero mezclamos qué etiqueta de
grupo tiene cada red. Recalculamos el score muchas veces. Si el score real es
más extremo que la mayoría de los mezclados, las etiquetas de grupo importan.
El p-valor es la proporción de mezclas al menos tan extremas como el dato real.

### Tasa de error por familia (FWER)

La probabilidad de cometer al menos una selección falsa en todo el conjunto de
tests por arista. Controlar la FWER mantiene los falsos positivos muy poco
probables, incluso con miles de aristas.

### Corrección de Holm

Holm controla la tasa de error por familia. Ordena los p-valores por arista de
menor a mayor y usa un corte más estricto para los más chicos. Es nuestro
default porque es válido incluso cuando los tests por arista están
correlacionados, como pasa acá. En criollo: Holm obliga a cada arista a
"ganarse" la significancia frente al hecho de que testeamos miles de aristas.

### Tasa de falsos descubrimientos (FDR)

La proporción esperada de falsos positivos entre las aristas que sí
seleccionamos. La FDR es más relajada que la FWER: acepta unas pocas aristas
falsas para encontrar más verdaderas.

### Corrección BH (Benjamini-Hochberg)

Un método de FDR. Tiene más potencia que Holm, pero su garantía habitual
necesita que los tests tengan una estructura de dependencia positiva. Por eso
es opcional, no default.

### Corrección BY (Benjamini-Yekutieli)

Un método de FDR que sigue siendo válido bajo cualquier dependencia entre
tests. Es más estricto que BH. Lo ofrecemos cuando se quiere control de FDR sin
el supuesto de BH.

### Test a dos colas y `|T|`

Dos colas significa que nos importa una diferencia en cualquier dirección.
Comparamos el tamaño del score, `|T|`, así un efecto fuerte cuenta como extremo
tanto si el score es muy positivo como si es muy negativo.

### P-valor Monte Carlo y la regla del "más uno"

Cuando no podemos listar todas las mezclas posibles, sorteamos muchas al azar.
Usamos `(extremos + 1) / (B + 1)`, donde `B` es la cantidad de mezclas. El "más
uno" cuenta el dato real como un resultado válido, así el p-valor nunca puede
ser cero.

### Límite superior de confianza de Wilson

Se usa solo en la sección de validación. Después de muchas corridas de
simulación, da una cota superior cautelosa para una tasa de error. Pedimos que
esa cota se mantenga por debajo de `0.065`, para tener confianza de que la tasa
real está cerca del objetivo `0.05`.

### `-log10(p-valor)`

Es solo un truco de visualización en el Panel B. Convierte p-valores diminutos
en puntos altos, así los resultados más fuertes quedan más arriba y son más
fáciles de ver.

## Por qué se agregó Holm (y qué error corremos sin él)

**El problema:** comparaciones múltiples. En el análisis balanceado de ABIDE
testeamos 6.441 aristas. Si usamos un corte de `0.05` sin corregir, y no
hubiera ninguna diferencia real, igual esperaríamos alrededor de
`6.441 x 0.05 ≈ 322` aristas "significativas" que son falsas. Con miles de
tests, encontrar algo por azar deja de ser raro y pasa a ser casi seguro.

**Qué pasaba antes:** en el paquete viejo la corrección se usaba solo para
ordenar las aristas, no para decidir cuáles entraban en el conjunto "crítico".
O sea, no había garantía a nivel de familia. Con eso, el conjunto crítico podía
llenarse de falsos positivos y aun así presentarse como aristas importantes.

**Qué hace Holm:** controla la probabilidad de que haya aunque sea una arista
falsa en todo el conjunto (la FWER). Vale bajo cualquier dependencia entre los
tests, que es justo nuestro caso, porque las aristas comparten sujetos.

**Cuánto importa, con nuestros propios números de simulación** (comparador sin
corregir contra Holm, en las cuatro alternativas):

- Sin corregir: detecta cerca del 60–66% de las aristas verdaderas, **pero** en
  el 37%–53% de los estudios reporta al menos una arista falsa, y entre el 9% y
  el 15% de las aristas reportadas son falsas.
- Con Holm: la chance de reportar alguna arista falsa baja a alrededor del 1%, y
  la proporción de falsas entre las seleccionadas baja a menos del 1%.

**En criollo:** sin Holm "encontrás más", pero una buena parte es mentira. Con
Holm "encontrás menos", pero lo que marcás es confiable.

**Conexión con ABIDE:** las 3 aristas con `p` crudo `< 0.001` no sobreviven a
Holm. Sin Holm las habríamos llamado críticas; con Holm entendemos que no hay
evidencia suficiente para una conexión puntual.

## Cómo interpretar cada resultado

### "La mayor tasa de rechazo global bajo el nulo completo fue 0.051"

- **Qué es:** cuando NO hay diferencia real, cuántas veces el test global igual
  rechaza.
- **Objetivo:** `0.05`.
- **Lectura:** `0.051 ≈ 0.05`. El test está bien calibrado; no infla los falsos
  positivos. **Bueno.** Si este número fuera, por ejemplo, `0.20`, el test
  estaría roto.

### "La mayor tasa de error por familia de Holm fue 0.015"

- **Qué es:** cuando no hay ninguna arista distinta, la chance de marcar aunque
  sea una arista falsa.
- **Objetivo:** `<= 0.05`.
- **Lectura:** `0.015` está muy por debajo de `0.05`. Holm es seguro y un poco
  conservador. **Bueno.**

### "Tasa media de verdaderos positivos (TPR) de Holm = 0.191"

- **Qué es:** de las aristas que sí difieren de verdad, qué fracción
  detectamos.
- **Lectura:** cerca del 19%. Es **poca potencia**. No es un error: es el costo
  de ser estricto con muestras de 15–20 por grupo. Con más sujetos sube; también
  sube si se usa BY (control de FDR, más permisivo).
- **¿Malo?** No es "malo", es honesto. Preferimos no inventar aristas. Lo
  declaramos como limitación en el paper.

### "Precisión media de Holm = 0.716"

- **Qué es:** entre las aristas seleccionadas, qué fracción son verdaderas
  (promediando además las corridas que no seleccionan nada, que cuentan como 0).
- **Cuidado con la lectura fácil:** NO significa "el 28% de las seleccionadas
  son falsas". La proporción real de falsos descubrimientos es de alrededor de
  `0.004`, y la FWER de alrededor de `0.01`.
- **Qué pasa en realidad:** en aproximadamente el 72% de las corridas Holm
  selecciona al menos una arista, y cuando selecciona, casi todas son
  verdaderas. En cerca del 28% no selecciona nada (se pierde todo). Ese `0.716`
  mezcla las dos cosas.
- **Lectura:** cuando Holm marca algo, es confiable. **Bueno.**

### Las dos juntas (TPR + precisión)

El perfil es: "marca poco, pero lo que marca es casi siempre correcto". Es
justo lo que queremos en un contexto clínico, para no sobre-afirmar.

### Los números de ABIDE (global contra aristas)

- **Global significativo** (completo `T = -44.7`, `p = 0.001`; balanceado
  `T = -27.4`, `p = 0.01`): hay un patrón distribuido asociado al diagnóstico,
  en esta muestra y este pipeline.
- **Ninguna arista sobrevive Holm/BY:** no hay evidencia para una conexión
  puntual.
- **Por qué es coherente:** si el efecto está repartido en muchas aristas
  chiquitas, el test global lo ve sumando todo, pero cada arista sola no llega
  al umbral estricto. La TPR baja de las simulaciones anticipa exactamente este
  escenario.

## Preguntas que tu co-autor puede hacer

### ¿Por qué se sacó el viejo resultado de `identify_critical_links()`?

Su orden de aristas y su punto de parada usaban los mismos datos. El resultado
viejo podía servir como descripción, pero no era inferencia válida por arista.
La versión 1.0 usa tests por arista corregidos para la selección y conserva la
ablación solo como descripción.

### ¿Qué se mezcla exactamente?

Las redes de participantes enteras reciben nuevas etiquetas de grupo. El orden
de nodos, los valores de las aristas y los tamaños de grupo quedan fijos.

### ¿El paquete borra aristas de los datos de entrada?

No. Las matrices originales no se editan. La ablación resta contribuciones de
aristas guardadas al score.

### ¿Por qué se usa solo el triángulo superior?

Las redes son no dirigidas. El triángulo inferior repite las mismas aristas.
Contar los dos duplicaría la distancia.

### ¿Por qué se usa `|T|`?

Con tamaños desiguales o más de dos grupos, el score no siempre se mueve en una
dirección fija. `|T|` da un test a dos colas.

### ¿Por qué se sacó la constante de normalización `a`?

Multiplicar todos los scores por el mismo número positivo no cambia su orden de
aleatorización ni el p-valor.

### ¿Por qué el test global puede ser significativo cuando ninguna arista queda seleccionada?

Muchas diferencias chicas pueden combinarse en un patrón global fuerte. Una
arista sola tiene que pasar un umbral mucho más estricto porque se testean
6.441 aristas.

### ¿Por qué Holm en lugar de BH?

Holm controla la probabilidad de que haya cualquier arista falsa seleccionada
bajo cualquier dependencia. BH puede tener más potencia, pero necesita
supuestos extra de dependencia. BY está disponible cuando se prefiere control
de FDR bajo dependencia general.

### ¿Por qué la tasa de verdaderos positivos en simulación es baja?

Las muestras tienen solo 15 a 20 redes por grupo y Holm es estricto. La tasa
baja es el costo de un control fuerte de falsos positivos, no una falla oculta.

### ¿Por qué se fuerza a cada red a tener 10% de densidad en ABIDE?

Esto saca la cantidad total de aristas como principal diferencia entre grupos.
El test se enfoca entonces en dónde ocurren las aristas. El resultado igual
depende de la elección del 10%.

### ¿Por qué se usa solo el sitio NYU de ABIDE?

Mezclar sitios agregaría diferencias de escáner y de sitio. Un solo sitio da un
ejemplo más limpio, pero también limita cuánto se puede generalizar el
resultado.

### ¿Por qué se crea una muestra ABIDE balanceada?

Los grupos originales difieren en sexo y movimiento. La muestra balanceada
chequea si el resultado global se mantiene después de igualar más el sexo, la
edad y el movimiento medidos.

### ¿El balanceo saca toda la confusión?

No. Solo mejora el balance de las variables medidas y los bins usados. Pueden
quedar diferencias no medidas.

### ¿Podemos decir que el autismo causa la diferencia de red?

No. Son datos observacionales. Solo podemos decir que la etiqueta de
diagnóstico está asociada con el patrón de aristas para esta muestra y este
pipeline.

### ¿Podemos llamar biomarcadores a las aristas crudas más fuertes?

No. Ninguna sobrevive a Holm ni a BY. Se muestran solo por contexto anatómico.

### ¿El paquete maneja tres o más grupos?

Sí. El test global soporta varios grupos. Los tests por arista son ómnibus:
dicen que al menos un grupo difiere, no qué par difiere.

### ¿Los grupos pueden tener tamaños distintos?

Sí. El test de aleatorización a dos colas se eligió en parte para soportar
tamaños de grupo desiguales.

### ¿Ya podemos usar redes con peso o dirigidas?

No. Eso necesitaría un nuevo contrato de datos, fórmulas, tests y validación.

### ¿Cuál es la principal contribución de software?

Clases de R validadas, conclusiones globales y por arista separadas, salidas
reproducibles, aleatorización por bloques y actualizaciones exactas por sumas
de prefijos.

### ¿Qué conviene evitar decir en la reunión?

Evitar:

- "Las aristas seleccionadas causan autismo."
- "El test global encuentra todo tipo de diferencia de red."
- "Ninguna arista significativa quiere decir que no hay diferencia de aristas."
- "El resultado de `953x` es la aceleración total del paquete."
- "El resultado balanceado de ABIDE saca toda la confusión."

## Decisiones para confirmar juntos

- ¿Estamos los dos cómodos con el umbral proporcional del 10%?
- ¿La redacción de ABIDE es lo bastante cuidadosa sobre asociación y confusión?
- ¿Acordamos que no se afirma ninguna arista puntual como biomarcador?
- ¿Acordamos Holm como default y BY/BH como opciones?
- ¿Acordamos que la ablación es solo descriptiva?
- ¿Los límites del paquete están enunciados con suficiente claridad?
