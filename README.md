# Diseño de nuevas enzimas mediante scaffolding del centro activo de enzimas de estructura conocida
El objetivo principal de este trabajo es, utilizando el modelo [RF*diffusion*](https://github.com/RosettaCommons/RFdiffusion), diseñar nuevas enzimas a partir de estructuras ya conocidas.
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

Para definir el motivo en torno al que queremos generar la nueva proteína, consideraremos los residuos cuya distancia al sustrato sea inferior a una longitud propuesta. Definimos la distancia entre residuo y sustrato como la mínima distancia de los átomos del residuo (capturados en el
PDB) y los átomos del sustrato.

Para ello usamos el código:
```python
from rfdiffusion.inference import utils as iu
from mytools import utils as myu

# Seleccionar motivo
# dist < 5A
pdb = iu. parse_pdb ("../TFM/RFdiffusion/inputs/6vdz.pdb", parse_hetatom = True )
substrateName = ' HEC'
chain = 'A'
distMotif = 5.0
CA = True

myu.motif_substr(pdb,substrateName,chain,distMotif,CA)
```
Para una distancia de 5&#x212b;, el motivo estaría definido por los residuos 
[191, 194, 195, 198, 231, 234, 238, 266, 268, 269, 271, 274, 277, 278, 281, 282, 306, 309, 310, 313, 316, 317, 318, 319, 320].

### Aplicar *motif-scaffolding* al motivo para generar la estructura de la proteína.

Para realizar la difusión inversa sobre el pdb ```6vdz``` aplicamos los comandos:

- Variante 1
```
scripts/run_inference.py \
inference.output_prefix=../TFM/RFdiffusion/outputs/6vdz_M/6vdz_M0_5.0A \
inference.input_pdb=../TFM/RFdiffusion/inputs/6vdz.pdb \
'contigmap.contigs=[173/A191-191/2/A194-195/2/A198-198/32/A231-231/2/A234-234/3/A238-238/27/A266-266/1/A268-269/1/A271-271/2/A274-274/2/A277-278/2/A281-282/23/A306-306/2/A309-310/2/A313-313/2/A316-320]' \
potentials.guide_scale=1 \
'potentials.guiding_potentials=["type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1","type:monomer_ROG,weight:1,min_dist:5"]' \
potentials.substrate=HEC \
inference.num_designs=10
```

- Variante 2:
```
scripts/run_inference.py \
inference.output_prefix=../TFM/RFdiffusion/outputs/6vdz_M/6vdz_M1_5.0A \
inference.input_pdb=../TFM/RFdiffusion/inputs/6vdz.pdb \
'contigmap.contigs=[173/A191-191/2/A194-195/2/A198-198/32/A231-231/2/A234-234/3/A238-238/27/A266-266/1/A268-269/1/A271-271/2/A274-274/2/A277-278/2/A281-282/23/A306-306/2/A309-310/2/A313-313/2/A316-320]' \
potentials.guide_scale=2 \
potentials.guide_decay="quadratic" \
'potentials.guiding_potentials=["type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1","type:monomer_ROG,weight:1,min_dist:5"]' \
potentials.substrate=HEC \
inference.num_designs=10
```
- Variante 3
```
scripts/run_inference.py \
inference.output_prefix=../TFM/RFdiffusion/outputs/6vdz_M/6vdz_M2_5.0A \
inference.input_pdb=../TFM/RFdiffusion/inputs/6vdz.pdb \
'contigmap.contigs=[173/A191-191/2/A194-195/2/A198-198/32/A231-231/2/A234-234/3/A238-238/27/A266-266/1/A268-269/1/A271-271/2/A274-274/2/A277-278/2/A281-282/23/A306-306/2/A309-310/2/A313-313/2/A316-320]' \
potentials.guide_scale=1 \
potentials.guide_decay="quadratic" \
'potentials.guiding_potentials=["type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1","type:monomer_ROG,weight:1,min_dist:5"]' \
potentials.substrate=HEC \
inference.num_designs=10
```
En las tres variantes usamos los potenciales:
- ```substrate_contacts```
- ```monomer_ROG```

### Aplicar *inverse folding* para obtener la secuencia de aminoácidos

## Medir la calidad de las enzimas generadas

### Generar usando *AlphaFold* la estructura de la proteína

### Analizar la calidad del *bolsillo* generado en torno al sustrato
