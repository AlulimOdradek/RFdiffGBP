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
    
En primera instancia definimos y aplicamos el flujo de trabajo al monómero [6VDZ](https://www.rcsb.org/structure/6vdz), para después tratar de extender el procedimiento al díımero [7KQU](https://www.rcsb.org/structure/7kqu). En ambos casos se trata de oxidoreductasas con sustrato un grupo HEMO.

## Generar el backbone de la nueva enzima 
Para realizar el trabajo hemos utilizado el modelo de difusión descrito en [RF*diffusion*](https://github.com/RosettaCommons/RFdiffusion), por lo que  que hemos construido un entorno local tal como se describe en el GitHub.
Sobre este entorno hemos añadido el código ```mytools/utils.py``` y hemos modificado ```rfdiffusion/inference/utils.py```.

### Definir el motivo a partir del sustrato.

Para definir el motivo en torno al que queremos generar la nueva proteína, consideraremos los residuos cuya distancia al sustrato sea inferior a una longitud propuesta. Definimos la distancia entre residuo y sustrato como la mínima distancia de los átomos del residuo (capturados en el PDB) y los átomos del sustrato.

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
Para una distancia de 5&#x212b;, tras diversas pruebas solo hemos obtenido resultados satisfactorios para esta distancia, el motivo estaría definido por los residuos [191, 194, 195, 198, 231, 234, 238, 266, 268, 269, 271, 274, 277, 278, 281, 282, 306, 309, 310, 313, 316, 317, 318, 319, 320].



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
- ```substrate_contacts```, que trata de guiar la estructura diseñada para que se ajuste a un bolsillo complementario al substrato.
- ```monomer_ROG```, que favorece estructuras compactas, es decir, globulares.

A continuación filtramos los resultados obtenidos utilizando la información resumen que nos proporciona la función ```summ_pdbs()```, definida en ```mytools/utils.py```: 
- ```maxd``` < 65 5&#x212b;: Máxima distancia entre los Cα de los residuos de la enzima. Nos va a dar una idea de la globularidad de la proteína. Para el wild-type tenemos un valor de maxd = 63.70 5&#x212b, por lo que nos parece que un límite de 65 5&#x212b nos puede garantizar hasta cierto punto la globularidad de la enzima.
- ```rmsd``` < 0.3 5&#x212b;: RMSD (Root Mean Square Deviation) entre el motivo en el wild-type y el mismo motivo en el PDB diseñado, una vez que el segundo se ha superpuesto al primero. Para lograr la superposición se utiliza el algoritmo de Kabsch para calcular la matriz de rotación óptima de RMSD mínimo entre dos conjuntos de puntos pareados.

Obtenemos dos estructuras que cumplen con nuestros requisitos.

### Aplicar *inverse folding* para obtener la secuencia de aminoácidos

RF*diffusion* nos proporciona el backbone de la enzima y necesitammos predecir las cadenas laterales de la proteína. Para ello hacemos uso de técnicas de *Inverse Folding*, en concreto ProteinMPNN, Para implementar el procedimiento en local hemos seguido las instrucciones propuesteas por los autores en su [GitHub](https://github.com/dauparas/ProteinMPNN).

En nuestro caso hemos fijado el motivo y hemos generado 10 secuencias por cada estructura
seleccionada usando los comandos:

```
python ./helper_scripts/parse_multiple_chains.py \
--input_path=../TFM/ProteinMPNN/inputs/6vdz/motif_5.0A_4 \
--output_path=../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/parsed_pdbs.jsonl


python ./helper_scripts/assign_fixed_chains.py \
--input_path=../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/parsed_pdbs.jsonl \
--output_path=../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/assigned_pdbs.jsonl \
--chain_list="A"


python ./helper_scripts/make_fixed_positions_dict.py \
--input_path=../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/parsed_pdbs.jsonl \
--output_path=../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/fixed_pdbs.jsonl \
--chain_list "A" \
--position_list "174 177 178 181 214 217 221 249 251 252 254 257 260 261 264 265 289 292 293 296 299 300 301 302 303"


python ./protein_mpnn_run.py \
--jsonl_path ../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/parsed_pdbs.jsonl \
--chain_id_jsonl ../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/assigned_pdbs.jsonl \
--fixed_positions_jsonl ../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/fixed_pdbs.jsonl \
--out_folder ../TFM/ProteinMPNN/outputs/6vdz/motif_5.0A_4/ \
--num_seq_per_target 10 \
--sampling_temp "0.1" \
--seed 1966 \
--batch_size 1
```
Los tres primeros comandos preparan los ficheros de configuración ```parsed_pdbs.jsonl```, ```assigned_pdbs.jsonl``` y ```fixed_pdbs.jsonl```, este último define las posiciones que deben de quedar fijas, que después son utilizados por ```protein_mpnn_run.py```, que nos proporciona 10 ficheros FASTA para cada una de las estructuras almacenadas en el directorio de entrada proporcionado en el primer comando. 

Agrupamos los ficheros FASTA obtenidos en un único fichero utilizando la función ```read_fa()```contenida en ```mytools/utils.py```.

## Medir la calidad de las enzimas generadas

Tratamos de medir la calidad de la enzima en cuanto a su capacidad de reaccionar con el sustrato. Tratar de emular la reacción mediante dinámica molecular es muy complejo, por lo que optamos por trabajar con la complementariedad entre el sustrato y el bolsillo generado. Para ello, en primer lugar generamos mediante *AlphaFold2* las estructuras asociadas a las secuencias de aminoácidos obtenidas mediante inverse folding, para después analizar la calidad del bolsillo generado en estas estructuras.

### Generar usando *AlphaFold* la estructura de la proteína
Para analizar la calidad de las secuencias obtenidas, construimos, mediante *AlphaFold2* las estructuras asociadas a dichas secuencias. Para ello hemos utilizado el método de instalacón [localcolabfold](https://github.com/YoshitakaMo/localcolabfold).

Usamos el parámetro ```--amber``` para el refinamiento de la estructura y ```--templates``` para que durante el proceso de predicción use plantillas de pdb y aplicamos el proceso al fichero fasta construido con los resultados de ProteinMPNN, utilizando un comando de la forma:

```
colabfold_batch --templates --amber ./AF2/inputs/6vdz/motif_5.0A_4/6vdz_Px_5.0A.fa ./AF2/outputs/6vdz/motif_5.0A_4_templates/
```

### Analizar la calidad del *bolsillo* generado en torno al sustrato

Filtramos 
