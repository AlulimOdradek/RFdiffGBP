import torch
import os
import numpy as np

# Utilidades definidas en el algoritmo de RFdiffusion
from rfdiffusion.inference import utils as iu
import rfdiffusion.chemical as ch

def center(m):
    """
    
    Parameters
    ----------
    m : numpy.matrix
        Matriz de coordenadas de un conjunto de puntos.

    Returns
    -------
    list
        Coordenadas del centro de masas del conjunto de puntos
        y matriz de coordenadas de los puntos con respecto al CdM.

    """
    # Devuelve la matriz centrada en el origen de coordenadas
    # Calcula el centro de masas de la matriz m
    center_of_mass_m = m.sum(1)/m.shape[1]
    # Mueve m al origen de coordenadas
    centered_m = m - center_of_mass_m
    return [center_of_mass_m, centered_m]


def sup(x, y):
    """
    Superimpose de las coordenadas 'x' e 'y'. 
    tiliza el algoritmo de Kabsch para calcular la matriz de rotación óptima 
    de RMSD mínimo entre dos conjuntos de puntos pareados.
    
    Construido a partir de https://www.biostars.org/p/486865/
    
    Parameters
    ----------
    x : numpy.matrix
        Matriz de las coordenadas de un conjunto de puntos.
    y : numpy.matrix
        Matriz de las coordenadas de un conjunto de puntos.

    Returns
    -------
    rmsd : numpy.float64
        RMSD medido entre los puntos de ambas matrices.
    desp: numpy.matrix
        Vector entre los centros de masa.
    u : numpy.matrix
        Matriz de rotación.

    """
    # Número de átomos
    N=x.shape[1]

    ########################
    
    # Centra x e y
    x=center(x)
    y=center(y)

    # Matriz de correlación
    r=y[1]*x[1].T

    # SVD((Singular Value Decompposition) de la matriz de correlación
    v, s, wt=np.linalg.svd(r)

    w=wt.T
    vt=v.T

    # Matriz de rotación
    u=w*vt

    # Chequea si hay reflexión
    if np.linalg.det(u)<0:
        z=np.diag((1,1,-1))
        u=w*z*vt
        s[2]=-s[2]

    # Calcula RMSD entre las coordenadas antes y después de rotar
    d=x[1]-u*y[1]
    d=d.A*d.A
    rmsd=np.sqrt(np.sum(d)/N)
    
    desp = x[0] - u*y[0]
    
    # Devuelve RMSD, el vector de traslación y la matriz de rotación
    return rmsd, desp, u


def motif_substr(pdb, substrateName, substrateChain, dist, CA=True):
    """
    Identificación del motivo.
    
    Parameters
    ----------
    pdb : dict
        Estructura generada por RFdiffusion a partir del pdb.
    nameSubstrate : string
        Nombre del sustrato capturado en el pdb.
    dist : floating point
        Distancia mínima al substrato para definir el motivo.
    CA : bool
        Con valor TRUE, calculamos las distancias desde los Ca.
        A FALSE, calculamos las distancias desde todos los átomos del residuo en el PDB

    Returns
    -------
    Lista de ids del motivo junto con la mínima distancia al sustrato.

    """
    idSubstrate = [pdb['info_het'].index(d) for d in pdb['info_het'] if d['name'] in substrateName 
                   and d['chain'] in substrateChain] 
    substrate = torch.tensor(pdb['xyz_het'][idSubstrate], dtype=torch.float64)    
    def minDist(k, CA):
        if CA:
            res = torch.tensor(pdb['xyz'][k,1], dtype=torch.float64)
            pdist = torch.nn.PairwiseDistance(p=2)
            d = pdist(res, substrate)             
        else:
            res = torch.tensor(pdb['xyz'][k,pdb['mask'][k]], dtype=torch.float64)
            d = torch.cdist(res, substrate, p=2)              
        return torch.min(d).tolist()
    
    return [[pdb['pdb_idx'][k], minDist(k, CA)]  for k in range(pdb['xyz'].shape[0]) if minDist(k, CA) < dist]


def maxDist(pdb):
    """
    Parameters
    ----------
    pdb : dict
        Estructura generada por RFdiffusion a partir del pdb..

    Returns
    -------
    float
        Máxima distancia entre CA de los residuos del PDB.

    """
    res = torch.tensor(pdb['xyz'][:,1], dtype=torch.float64)
    d = torch.cdist(res, res, p=2)    
    return torch.max(d).tolist()

def minDist(pdb, idx1, idx2):
    """
    Proporciona la mínima distancia entre los átomos de dos residuos

    Parameters
    ----------
    pdb : dict
        Estructura del PDB proporcionada por RFdiffusion.
    idx1 : int
        Id del residuo 1.
    idx2 : int
        id del residuo 2.

    Returns
    -------
    md : float
        Mínima distancia entre los átomos de los residuos.

    """
    R1 = torch.tensor(pdb['xyz'][idx1,pdb['mask'][idx1]])
    R2 = torch.tensor(pdb['xyz'][idx2,pdb['mask'][idx2]])
    
    md = torch.min(torch.cdist(R1,R2,p=2)).tolist()
    
    return(md)

def identity(pdb1, pdb2):
    """Devuelve el porcentaje de identidad de las dos secuencias, 
    para nuestro caso simple"""
    if pdb1['seq'].size == pdb2['seq'].size:
        return sum(pdb1['seq']==pdb2['seq'])/pdb1['seq'].size
    elif pdb1['seq'].size > pdb2['seq'].size:
        return sum(pdb1['seq'][0:pdb2['seq'].size]==pdb2['seq'])/pdb2['seq'].size
    else:
        return sum(pdb1['seq']==pdb2['seq'][0:pdb1['seq'].size])/pdb1['seq'].size
        

def read_fa(path, pathw, filew):
    """
    Construye un fichero FASTA resumen a partir de un listado de ficheros FASTA independientes.
    
    Parameters
    ----------
    path : string
        Path en el que se encuentra los ficheros FASTA que hemos de leer.
    pathw : string
        Path en donde desaemos escribir el fichero FASTA resumen de todos los leídos.
    filew : string
        Nombre del fichero FASTA resumen.

    Returns
    -------
    None.

    """
    
    # Modificar f.startwith
    pathw = os.path.join(pathw, filew)
    files = os.listdir(path)
    with open(pathw, "w") as fw:
        for f in files:
            ff = os.path.join(path, f)
            if os.path.isfile(ff) and ff.endswith('.fa'):
                with open(ff, 'r') as fff:
                    lines = fff.readlines()
                    for index, l in enumerate(lines):
                        if index > 1:
                            if l[0] != '>':
                                fw.write(l)
                            else:
                                sample = l.split(', ')[1]
                                eq = sample.find('=')
                                sample = int(sample[eq + 1:]) - 1
                                fw.write('>' + f[:-3] + '_s' + str(sample) + '\n')                                    


def summ_pdbs(path, idx_origin, idx_target, pdb_target, substrateName):
    """
    Nos ofrece un resumen de ciertos datos de los PDBs contenidos en un directorio.

    Parameters
    ----------
    path : string
        Path en el que están contenidos los PDBs.
    idx_origin : list
        Id de los residuos del motivo en el pdb origen (en el movimiento)
    idx_target : list
        Id de los residuos del motivo en el pdb destino
    pdb_target : dict
        Estructura generada por RFdiffusion a partir del pdb (del wild-type).
    nameSubstrate : string
        Nombre del sustrato en los registros HETATM del PDB.

    Returns
    -------
    pdbs : dict
        Diccionario en el que a cada fichero PDB (la key es el nombre del fichero) 
        asigna una lista de valores:
            - rmsd: RMSD entre dicho pdb y el pdb_target.
            - mindCA: Mínima distancia entre los CA del motivo y el sustrato.
            - mind: mínima distancia de los residuos al sustrato
            . md: distancia entre dos motivos concretos. Usado para 6vdz
            - maxd: Máxima distancia entre los CA del PDB. Para estimar la globularidad.
            - plddt_mean: Valor medio de los valores de plddt de todos los residuos del pdb.
            - ident: índice de identidad entre el wild-type y la proteina diseñada.

    """
    a = np.matrix(pdb_target['xyz'][idx_target,1]).T
    idSubstrate = [pdb_target['info_het'].index(d) for d in pdb_target['info_het'] if d['name'] in substrateName] 
    substrate = torch.tensor(pdb_target['xyz_het'][idSubstrate])
    files = os.listdir(path)
    pdbs = dict()
    for f in files:
        fp = os.path.join(path, f)
        if os.path.isfile(fp) and f.endswith('.pdb'):
            pdb_diff = iu.parse_pdb(fp, parse_hetatom=False)
            plddt = np.array(pdb_diff['plddt'])
            plddt_mean = np.mean(plddt, dtype=float)
            b = np.matrix(pdb_diff['xyz'][idx_origin,1]).T
            rmsd, v, u = sup(a,b)
            motif_diff_mm = (u*b + v).T
            motif_diff_mm = torch.tensor(motif_diff_mm, dtype=torch.float64)
            mindCA = torch.min(torch.cdist(motif_diff_mm, substrate, p=2)).tolist()
            maxd = maxDist(pdb_diff)
            pdb_diff_all = []
            for k in range(pdb_diff['xyz'].shape[0]):
                for i in range(pdb_diff['xyz'][k,pdb_diff['mask'][k]].shape[0]):
                    pdb_diff_all.append(pdb_diff['xyz'][k,pdb_diff['mask'][k]][i])
            pdb_diff_all = np.matrix(pdb_diff_all)
            pdb_diff_all_mm = torch.tensor((u*pdb_diff_all.T + v).T, dtype=torch.float64)
            mind = torch.min(torch.cdist(pdb_diff_all_mm, substrate, p=2)).tolist()
            md = minDist(pdb_diff, 256, 295)
            ident = identity(pdb_target, pdb_diff)
            pdbs[f] = [rmsd, mindCA, mind, md, maxd, plddt_mean, ident]
    return pdbs 


def pdb2aanum(pdbFile):
    resi = dict()
    with open(pdbFile, 'r') as f:
        for line in f:
            if line.startswith("ATOM"): #¿Sólo ATOM?
                campos = line.split()
                numRes = campos[5]
                nameRes = campos[3]
                chain = campos[4]
                try:
                    resi[chain+numRes] = ch.aa2num[nameRes]
                except KeyError:
                    None
                
    return torch.Tensor([i for i in resi.values()]).to(int)




