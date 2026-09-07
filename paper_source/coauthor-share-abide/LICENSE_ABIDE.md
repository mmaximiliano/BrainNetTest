# Licencia y términos de uso (Usage Agreement) de ABIDE

## Datos clínicos

La fuente primaria es ABIDE I. La página oficial incluye una sección
`Usage Agreement`:

<http://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html>

Allí se establece:

- uso irrestricto para investigación no comercial;
- registración para acceder a los datos raw por NITRC/INDI;
- obligación de identificar los datasets/sites usados;
- pedido de reconocer las fuentes de funding correspondientes; y
- Creative Commons Attribution-NonCommercial-ShareAlike.

Licencia indicada por ABIDE:

<https://creativecommons.org/licenses/by-nc-sa/3.0/>

## Derivados preprocesados

Este pipeline usa ABIDE Preprocessed del Preprocessed Connectomes Project:

- proyecto: <http://preprocessed-connectomes-project.org/abide/>
- downloads: <http://preprocessed-connectomes-project.org/abide/download.html>
- publicación a citar:
  <https://doi.org/10.3389/conf.fninf.2013.09.00041>

Los archivos usados por el script están en un bucket S3 público. Eso permite
descargarlos sin login. PCP pide citar su iniciativa/publicación cuando se usan
los derivatives. La licencia BSD que puede verse en la ficha NITRC de PCP
describe el software/recurso y no debe interpretarse como una relicencia
comercial de los datos clínicos ABIDE.

El acceso a los datos ABIDE crudos por NITRC requiere registración bajo los
términos oficiales ABIDE/INDI. Que el derivative preprocesado sea accesible por
HTTP no elimina los términos del dataset original.

Por prudencia, los conectomas binarios derivados deben tratarse con las mismas
restricciones: uso no comercial, attribution y ShareAlike.

## Código

El package BrainNetTest y los scripts escritos para el análisis usan licencia
MIT. La licencia MIT cubre el software, no relicencia los datos clínicos.

## Qué se puede compartir

Se pueden compartir:

- scripts;
- parámetros del pipeline;
- tablas agregadas;
- resultados estadísticos no identificables; y
- derivados, solamente respetando los términos ABIDE/CC BY-NC-SA.

No se debe:

- usar los datos con fines comerciales;
- intentar reidentificar participantes;
- eliminar las citas/attribution;
- presentar los IDs anonimizados como identidad real; ni
- asumir que la licencia MIT del package reemplaza la licencia del dataset.

## Citaciones mínimas

1. Di Martino A et al. (2014), ABIDE:
   <https://doi.org/10.1038/mp.2013.78>
2. Craddock C et al. (2013), PCP:
   <https://doi.org/10.3389/conf.fninf.2013.09.00041>
3. Tzourio-Mazoyer N et al. (2002), AAL:
   <https://doi.org/10.1006/nimg.2001.0978>

## Funding acknowledgements para el site NYU

La página oficial de ABIDE I pide reconocer las fuentes de funding de los
datasets usados. Para NYU lista:

- NIH: `K23MH087770`, `R21MH084126`, `R01MH081218`, `R01HD065282`;
- Autism Speaks;
- The Stavros Niarchos Foundation;
- The Leon Levy Foundation; y
- el endowment de Phyllis Green and Randolph Cowen.

El acknowledgement general de ABIDE I también reconoce
`NIMH K23MH087770`, Leon Levy Foundation, Joseph P. Healy,
Stavros Niarchos Foundation y `NIMH R03MH096321`.

Esta nota resume los términos técnicos encontrados en las fuentes oficiales;
no constituye asesoramiento legal.
