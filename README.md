# Diseño de nuevas enzimas mediante scaffolding del centro activo de enzimas de estructura conocida
El objetivo principal de este trabajo es, partiendo del modelo [RF*diffusion*](https://github.com/RosettaCommons/RFdiffusion), diseñar nuevas enzimas partiendo de estructuras ya conocidas.
Para ello, hemos seguido el siguiente procedimiento:
- [Generar el backbone de la nueva enzima](#generar-el-backbone-de-la-nueva-enzima):
  - [Definir el motivo a partir del sustrato](#definir-el-motivo-a-partir-del-sustrato).
  - [Aplicar *motif-scaffolding* al motivo para generar la estructura de la proteína](#aplicar-motif-scaffolding-al-motivo-para-generar-la-estructura-de-la-proteína).
  - [Aplicar *inverse folding* para obtener la secuencia de aminoácidos](#aplicar-inverse-folding-para-obtener-la-secuencia-de-aminoácidos).
- [Medir la calidad de las enzimas generadas](#medir-la-calidad-de-las-enzimas-generadas):
  - [Generar usando *AlphaFold* la estructura de la proteína](#generar-usando-alphafold-la-estructura-de-la-proteína).
  - [Analizar la calidad del *bolsillo* generado en torno al sustrato](#analizar-la-calidad-del-bolsillo-generado-en-torno-al-sustrato).
    
En primera instancia definiremos y aplicaremos el flujo de trabajo al monómero [6VDZ](https://www.rcsb.org/structure/6vdz), para después tratar de extender el procedimiento al díımero [7KQU](https://www.rcsb.org/structure/7kqu). En ambos casos se trata de oxidoreductasas con sustrato un grupo HEMO.

## Generar el backbone de la nueva enzima 
Para realizar el trabajo hemos utilizado el modelo de difusión descrito en [RF*diffusion*](https://github.com/RosettaCommons/RFdiffusion), por lo que  que hemos construido un entorno local tal como se describe en el GitHub.
Sobre este entorno hemos añadido el código ```mytools/utils.py``` y hemos modificado ```rfdiffusion/inference/utils.py```.

### Definir el motivo a partir del sustrato.

### Aplicar *motif-scaffolding* al motivo para generar la estructura de la proteína.

### Aplicar *inverse folding* para obtener la secuencia de aminoácidos

## Medir la calidad de las enzimas generadas

### Generar usando *AlphaFold* la estructura de la proteína

### Analizar la calidad del *bolsillo* generado en torno al sustrato
