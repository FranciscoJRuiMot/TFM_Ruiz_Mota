import argparse
import sys
from pathlib import Path

import pandas as pd
import string
import subprocess
import re

import yaml
import os

import networkx as nx

### Funciones auxiliares
def create_env_inflow(inflows, q_organisms, organisms, dilution, resource_map): #hecho? Aquí ahora los inflows son diccionarios
    organisms_copy = organisms.copy()
    added_resources = []
    with open('head_environment.cfg', 'w') as new_file:
        new_file.write('# set external resource:\n')
        for organism in q_organisms:
            res = resource_map[organism]
            if res not in added_resources:
                #initial =  round(inflows[organism]/dilution) #ESTO CAMBIADO, EL INITIAL LO MEJOR ES QUE SEA 0
                initial = 0
                line = f'RESOURCE {res}:initial={initial}:inflow={inflows[organism]}:outflow={dilution}\n'
                new_file.write(line)
                added_resources.append(res)

        new_file.write('# set resources as by-products:\n')
        for organism in organisms_copy:
            res = resource_map[organism]
            if res not in added_resources:
                # outflow applied only to species acquiring external resources:
                #line = f'RESOURCE res{organism}:initial=0:inflow=0:outflow=0\n'
                # outflow applied to all species:
                line = f'RESOURCE {res}:initial=0:inflow=0:outflow={dilution}\n'
                new_file.write(line)
                added_resources.append(res)

        new_file.write('# set interproducts:\n')
        organisms_copy.extend(q_organisms)
        first_tasks = []
        for organism in organisms_copy:
            line = f'RESOURCE inter{organism}:initial=0\n'
            new_file.write(line)
            first_tasks.append(organism.split('_')[0])

        check_tasks = list(set(first_tasks))
        final_binary = []
        for task in check_tasks:
            line = f'REACTION {task} {task.lower()} process:value=0:type=pow requisite:max_count=1:max_tot_count=1\n'
            new_file.write(line)
            final_binary.append(0)

        return(final_binary)

def create_reactions_from_network_dp(df, intake_map, resource_map, organisms, conversion_rate):
    binary_ids = []
    with open('environment_reactions.cfg', 'w') as new_file:
        for x in organisms:
            new_file.write(f"\n# {x} reactions\n")
            first_tasK = x.split('_')[0]
            second_task = x.split('_')[1]
            
            intake = intake_map[x]
            intake_resource = resource_map[x]
            salidas = df.loc[x]
            n_salidas = sum(salidas != 0)

            if n_salidas == 1:
                salida_value = salidas[salidas != 0].values[0]
                salida_name = salidas[salidas != 0].index[0]
                string_template = f'REACTION {x}_to_{salida_name} {second_task.lower()} process:resource={intake_resource}:value=0:type=pow:min={intake}:max={intake}:product={salida_name}:conversion={conversion_rate} requisite:max_count=1:reaction={first_tasK}'
                new_file.write(f'{string_template}\n')
                binary_ids.append(1)
            
            elif n_salidas > 1:
                salida_value = salidas[salidas != 0].values
                salida_name = salidas[salidas != 0].index
                string_inter = f'REACTION {x}_to_INTER {second_task.lower()} process:resource={intake_resource}:value=0:type=pow:min={intake}:max={intake}:product=inter{x}:conversion={conversion_rate} requisite:max_count=1:reaction={first_tasK}'
                new_file.write(f'{string_inter}\n')
                binary_ids.append(1)
                for i in range(n_salidas):
                    string_template = f'REACTION INTER_{x}_to_{salida_name[i]} {second_task.lower()} process:resource=inter{x}:value=0:type=pow:min={salida_value[i]}:max={salida_value[i]}:product={salida_name[i]}:conversion=1 requisite:reaction={x}_to_INTER:max_count=1'
                    new_file.write(f'{string_template}\n')
                    binary_ids.append(0)

            elif n_salidas == 0:
                string_template = f'REACTION {x} {second_task.lower()} process:resource={intake_resource}:value=0:type=pow:min={intake}:max={intake} requisite:max_count=1:reaction={first_tasK}'
                new_file.write(f'{string_template}\n')
                binary_ids.append(1)
    
    binary_ids.reverse()
    return binary_ids

def actualizar_required_reactions(ruta_script, nuevo_binario):
    with open(ruta_script, 'r', encoding='utf-8') as file:
        lineas = file.readlines()

    with open(ruta_script, 'w', encoding='utf-8') as file:
        for linea in lineas:
            if "-set REQUIRED_REACTIONS" in linea:
                # Reemplaza únicamente el valor binario
                linea = f"    -set REQUIRED_REACTIONS {nuevo_binario} \\\n"
            file.write(linea)

def create_events_from_organisms(organisms, id_organisms, chemo_id, binary_list, end_time, invasion_time):
    #organisms es un diccionario de fenotipo:genoma
    #id_organisms es un diccionario de id:fenotipo
    with open('events.cfg', 'w') as file:
        pos = 0
        file.write(f'u begin SetEnvironmentInputs {binary_list[0]} {binary_list[1]} {binary_list[2]}\n')
        for id, phenotype in id_organisms.items():
            if id in chemo_id:
                file.write(f'#{phenotype}\n')
                file.write(f'u begin InjectSequence {organisms[phenotype]} {pos} {pos+1}\n')
                pos +=1
            else:
                file.write(f'#{phenotype}\n')
                for update in range(100, invasion_time+1, 100): #lo de sumar +1 mirarlo tmb para otro bucle que tiene +100
                    file.write(f'u {update} InjectSequence {organisms[phenotype]} {pos} {pos+1}\n')
                pos +=1
        
        file.write('u 0:100:end SavePopulation\n')
        
        file.write(f'u {end_time} Exit')

def save_result_abundance(out_path,
                   phenotype_to_sequence,
                   data_dir="data",
                   end_point=10100):

    timepoints=range(100, end_point, 100)
    abundance_by_id = pd.DataFrame()
    id_to_sequence = {}

    for t in timepoints:
        try:
            column_names = []
            rows_values = []

            with open(f"{data_dir}/detail-{t}.spop", "r") as f:
                for line in f:
                    # Nombres de columnas
                    m = re.match(r"^#\s+\d+:\s+(.*)$", line)
                    if m:
                        column_names.append(m.group(1))
                        continue

                    # Datos
                    if not line.startswith("#") and line.strip() != "":
                        expected = len(column_names)
                        parts = line.strip().split(maxsplit=expected - 1)
                        if len(parts) < expected:
                            parts += [""] * (expected - len(parts))
                        rows_values.append(parts)

            if not rows_values:
                continue

            detail_df = pd.DataFrame(rows_values, columns=column_names)

            # Guardar secuencias
            for _, row in detail_df.iterrows():
                oid = int(row["ID"])
                seq = row["Genome Sequence"]
                id_to_sequence[oid] = seq

            # Abundancia por ID
            for oid, seq in id_to_sequence.items():
                mask = detail_df["ID"].astype(int) == oid
                abundance = detail_df.loc[mask, "Number of currently living organisms"]
                abundance_by_id.loc[t, oid] = int(abundance.iloc[0]) if len(abundance) else 0

        except FileNotFoundError:
            pass

    abundance_by_id = abundance_by_id.fillna(0).astype(int)

    id_to_phenotype = {}
    for oid, seq in id_to_sequence.items():
        phen = "UNKNOWN"
        for phenotype, pheno_seq in phenotype_to_sequence.items():
            if seq == pheno_seq:
                phen = phenotype
                break
        id_to_phenotype[oid] = phen

    abundance_pheno = abundance_by_id.rename(columns=id_to_phenotype)
    abundance_pheno = abundance_pheno.groupby(level=0, axis=1).sum()

    abundance_pheno.index.name = "update" #ESTO ES NUEVO, AÑADIR AL GIT
    abundance_pheno.to_csv(out_path)
    print(abundance_pheno)
    print(f"[OK] Abundancia por fenotipos guardada en: {out_path}")

    return abundance_pheno
###

def resources_to_species_matrix(pheno_df, df_organisms, args):
    resources_table = pd.read_csv(args.resources_table)

    # Mapeo species_id (A,B,C) → timing (NAND_ORN, EQU_AND, ...)
    id_to_timing = df_organisms.set_index("species_id")["timing"]

    # Mapeo timing → resource_id que consume
    timing_to_resource = {
        id_to_timing[row["species_id"]]: row["resource_id"]
        for _, row in resources_table.iterrows()
    }

    # Mapeo resource_id → columna del pheno_df (1-indexado)
    resource_to_col = {
        i + 1: col for i, col in enumerate(pheno_df.columns)
    }

    species = pheno_df.index.tolist()
    matrix = pd.DataFrame(0.0, index=species, columns=species)

    for i in species:          # especie productora (fila)
        for j in species:      # especie consumidora (columna)
            resource_j = timing_to_resource[j]
            col = resource_to_col[resource_j]
            matrix.loc[i, j] = pheno_df.loc[i, col]

    return matrix

def network_node_properties(network_df, original_df, df_organisms, abundance_df, args):
    G = nx.from_pandas_adjacency(network_df, create_using=nx.DiGraph)
    timing_to_id  = df_organisms.set_index("timing")["species_id"]
    timing_to_gt  = df_organisms.set_index("timing")["gest_time"]
    
    # Última fila del CSV de abundancia (last 10 timepoint)
    last_abundance = abundance_df.iloc[-10:].mean()

    rows = []
    for node in G.nodes():
        sid      = timing_to_id[node]
        in_deg   = G.in_degree(node)
        out_deg  = G.out_degree(node)

        s_out = original_df.loc[node].sum()
        s_in  = network_df[node].sum()

        syntrophy_links = sum(1 for n in G.successors(node) if G.has_edge(n, node))

        rows.append({
            "species":              sid,
            "outgoing_resources":   s_out,
            "incoming_resources":   s_in,
            "rep_rate":             30 / timing_to_gt[node],
            "outgoing_degree":      out_deg,
            "incoming_degree":      in_deg,
            "abundance":            last_abundance.get(node, 0),
            "syntrophy_links":      syntrophy_links,
            "incoming_external_res": 1 if sid in args.chemo_resources else 0,
        })

    result = pd.DataFrame(rows)
    result["network_id"] = Path(args.adj_matrix).stem
    result.to_csv("node_properties.csv", index=False)

    return result

def save_config(args): #REFACTORIZAR
    if args.save_config is None:
        return

    config = vars(args)  # convierte el Namespace a diccionario
    config.pop("save_config")  # no tiene sentido guardarse a sí mismo

    # Convierte Path a string para que yaml pueda serializarlo
    config = {k: str(v) if isinstance(v, Path) else v for k, v in config.items()}

    with open(args.save_config, "w") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)

    print(f"[INFO] Configuración guardada en {args.save_config}")

def select_organisms(df, args): #REVISAR ESTO, ES COPIADO, PARA QUE COINCIDAN LAS COLUMNAS
    n = args.n_organisms

    # ── Filtros que no dependen de env_id ─────────────────────────────────────
    base_mask = pd.Series(True, index=df.index)

    if args.gt_range is not None:
        gt_min, gt_max = map(float, args.gt_range.split("_"))
        base_mask &= df["gest_time"].between(gt_min, gt_max)

    if args.delay_timing is not None:
        base_mask &= df["delay_timing"] <= args.delay_timing #ESTA COLUMNA HAY QUE AÑADIRLA A LA DB

    # ── env_ids a intentar ────────────────────────────────────────────────────
    # Si el usuario lo fija → solo ese. Si no → todos, ordenados por frecuencia
    if args.env_id is not None:
        env_ids_to_try = [args.env_id]
    else:
        env_ids_to_try = df["env_id"].value_counts().index.tolist()

    # ── Búsqueda por env_id ───────────────────────────────────────────────────
    for env_id in env_ids_to_try:
        candidates = df[base_mask & (df["env_id"] == env_id)]
        unique_phens = candidates["phen_label"].unique()

        if len(unique_phens) < n:
            continue  # no hay suficientes fenotipos, prueba el siguiente

        # Escoge n fenotipos al azar y saca un organismo de cada uno
        selected_phens = pd.Series(unique_phens).sample(n)
        result = pd.concat([
            candidates[candidates["phen_label"] == phen].sample(1)
            for phen in selected_phens
        ])

        result = result.reset_index(drop=True)
        result["species_id"] = list(string.ascii_uppercase[:n])
        result["timing"] = result["phen_label"].str.replace("+", "_", regex=False)

        result.to_csv("founders.csv", index=False)

        #Nuevo verbose
        if args.verbose:
            print("="*50)
            print("PARÁMETROS GENERALES DE LOS ORGANISMOS:")
            #print(f"Media del gestation time: {}")

        return result

    raise ValueError(f"Ningún env_id tiene {n} fenotipos distintos con los filtros aplicados")

def build_inflow_map(args):
    resources = args.chemo_resources
    inflows = args.inflow

    if len(inflows) == 1:
        # Un único valor → se asigna a todos los recursos
        inflow_val = float(inflows[0])
        return {r: inflow_val for r in resources}
    elif len(inflows) == len(resources):
        # Un valor por recurso → empareja en orden
        return {r: float(v) for r, v in zip(resources, inflows)}
    else:
        raise ValueError(
            f"--inflow tiene {len(inflows)} valores pero --chemo_resources tiene {len(resources)}. "
            f"Debe ser 1 valor o exactamente {len(resources)}.") #ESTO CAMBIAR Y TRADUCIR

def run_chemostat(args):

    args.inflow_map = build_inflow_map(args)

    network_df = pd.read_csv(args.adj_matrix, index_col=0)
    organisms = network_df.index.tolist()
    #args.prod_value = network_df.iloc[0].sum()
    args.n_organisms = len(organisms)

    if args.organisms_path is not None:
        df_organisms = pd.read_csv(args.organisms_path)
    else:
        db_organisms = pd.read_csv(args.db_file)
        df_organisms = select_organisms(db_organisms, args)

    env_data = pd.read_csv(args.envs_id_file)

    ## A PARTIR DE AQUÍ ES UNA COPIA DEL RUN_NETWORK
    organisms_sequences = {}
    id_organisms = {}
    env_id = df_organisms['env_id'].iloc[0]
    for index, org in df_organisms.iterrows():
        organisms_sequences[org['timing']] = org['sequence']
        id_organisms[org['species_id']] = org['timing']

    #Renombrar los csv
    pheno_df = network_df.copy()
    pheno_df.columns = pheno_df.columns.map(id_organisms)
    pheno_df.index = pheno_df.index.map(id_organisms)

    #Crear reacciones de la red y actualizar IDs binarios
    args.inflow_map = {id_organisms[k]: v for k, v in args.inflow_map.items()} #Nuevo del quimiostato
    q_organisms = [id_organisms[k] for k in args.chemo_resources]
    noq_organisms = [id_organisms[k] for k in organisms if k not in args.chemo_resources] 
    t_organisms = q_organisms + noq_organisms

    intake_map = {organism: pheno_df.loc[organism].sum() / args.conversion_rate for organism in pheno_df.index}
    
    #Estas dos lineas siguientes se añaden para evitar que haya organismos con 0 intake al no producir
    min_intake = min(v for v in intake_map.values() if v > 0)
    intake_map = {k: (1 - args.conversion_rate) * min_intake if v == 0 else v for k, v in intake_map.items()}

    if args.resources_table is None:
        resource_map = {pheno: f'res{pheno}' for pheno in t_organisms}
        resource_df = pheno_df.copy()
        resource_df.columns = resource_df.columns.map(lambda x: f'res{x}')
    else:
        table = pd.read_csv(args.resources_table)
        resource_map = {}
        for _, row in table.iterrows():
            pheno = id_organisms[row['species_id']]
            resource_letter = string.ascii_uppercase[row['resource_id'] - 1]
            resource_map[pheno] = f'res{resource_letter}'
        resource_df = pheno_df.copy()
        resource_df.columns = [f"res{letter}" for letter in string.ascii_uppercase[:len(resource_df.columns)]]
        print(resource_df)

    print(resource_map)
    final_binary = create_env_inflow(args.inflow_map, q_organisms, noq_organisms, args.dilution_rate, resource_map) #Aquí con este basta
    binary_ids = create_reactions_from_network_dp(resource_df, intake_map, resource_map, t_organisms, args.conversion_rate)
    binary_ids.extend(final_binary)
    binary_str = ''.join(map(str, binary_ids))
    actualizar_required_reactions('run-avida.sh', binary_str)

    #Crear events file
    create_events_from_organisms(organisms_sequences, id_organisms, args.chemo_resources, 
                                     env_data[env_data['seed']==env_id][['input_1_hex', 'input_2_hex', 'input_3_hex']].values.tolist()[0], 
                                     args.exp_time, args.invasion_time)
    
    with open('environment.cfg', 'w') as env_file, open('head_environment.cfg', 'r') as head_file, open('environment_reactions.cfg', 'r') as react_file:
        env_file.write(head_file.read())
        env_file.write(react_file.read())

    os.remove('head_environment.cfg')
    os.remove('environment_reactions.cfg')

    #Ejecutar Avida y guardar abundancias
    subprocess.run("./run-avida.sh", shell=True)
    abundance_df = save_result_abundance(out_path=args.output, phenotype_to_sequence=organisms_sequences, end_point=(args.exp_time)+100) #REVISAR QUE ESTO DEVUELVE EL RANGE BIEN
    if pheno_df.shape[0] != pheno_df.shape[1]:
        network_df_square = resources_to_species_matrix(pheno_df, df_organisms, args)
    else:
        network_df_square = pheno_df
    network_node_properties(network_df_square, pheno_df, df_organisms, abundance_df, args)
    
def parse_list(value):
    return [v.strip() for v in value.split(",")]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Descripción del programa",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # ── Argumentos obligatorios ────────────────────────────────────────────────
    parser.add_argument(
        "--adj_matrix",
        type=Path, default=None,
        help="Adjacency Matrix path") #aquí he quitado que sea obligatorio

    # ── Argumentos opcionales ──────────────────────────────────────────────────
    parser.add_argument(
        "--resources_table",
        type=Path, default=None,
        help="A table indicating which resource each species consumes")

    parser.add_argument(
        "--inflow",
        type=parse_list, default=["5"],
        help="Inflow rates for chemostat resources")
    
    parser.add_argument(
        "--chemo_resources",
        type=parse_list, default=["A"],
        help="Chemostat resources ID")
    
    parser.add_argument(
        "--dilution_rate",
        type=float, default=0.01,
        help="Dilution rate for chemostat (mortality rate)")
    
    parser.add_argument(
        "--env_id",
        type=int, default=None,
        help="Environment ID for organism selection")
    
    parser.add_argument(
        "--exp_time",
        type=int, default=10000,
        help="Duration of the experiment")
    
    parser.add_argument(
        "--conversion_rate",
        type=float, default=0.625,
        help="Ratio between production and consumption of an organism.")
    
    parser.add_argument(
        "--invasion_time",
        type=int, default=500,
        help="Time until the invading organisms are injected (in updates of 100 at a time)")
    
    parser.add_argument(
        "--gt_range",
        type=str, default=None,
        help="Range of selection gestation time")
    
    parser.add_argument(
        "--delay_timing",
        type=float, default=None,
        help="Maximum percentage of delay (0-1)")
    
    parser.add_argument(
        "--db_file",
        type=Path, default="./organisms.csv",
        help="Organism Database file path")
    
    parser.add_argument(
        "--envs_id_file",
        type=Path, default="envs_id.csv",
        help="Env description file path")
    
    parser.add_argument(
        "--organisms_path",
        type=Path, default=None,
        help="Selected organisms path (CSV)")
    
    parser.add_argument(
        "--output",
        type=str, default="abundances.csv",
        help="Ruta del fichero de salida") #ESTO HAY QUE CAMBIARLO A ABUNDANCES PATH
    
    parser.add_argument(
        "--save_config",
        type=str, default=None, metavar="PATH",
        help="Guarda los argumentos usados en un fichero YAML (ej: config.yaml)")
    
    parser.add_argument(
        "--config",
        type=Path, default=None, metavar="PATH",
        help="Carga argumentos desde un fichero YAML")

    parser.add_argument(
        "--verbose",
        action="store_true",           # flag booleano: presente = True
        help="Mostrar información detallada")

    # ── Primer parseo ──────────────────────────────────────
    args = parser.parse_args()

    # ── Si hay config, cargarla ────────────────────────────
    if args.config:
        with open(args.config) as f:
            config = yaml.safe_load(f)

        parser.set_defaults(**config)

        # volver a parsear para que CLI tenga prioridad
        args = parser.parse_args()

    if args.adj_matrix is None:
        parser.error("--adj_matrix es obligatorio (por CLI o por --config)")
    return args

def main(args):
    save_config(args)

    if args.verbose:
        print("=" * 50)
        print("CONFIGURACIÓN")
        print("=" * 50)
        print(f"  {'Adj. Matrix:':<20} {args.adj_matrix}")
        print(f"  {'DB File:':<20} {args.db_file}")
        print(f"  {'Output:':<20} {args.output}")
        print("-" * 50)
        print(f"  {'Inflow:':<20} {args.inflow}")
        print(f"  {'Chemo Resources:':<20} {args.chemo_resources}")
        print(f"  {'Dilution Rate:':<20} {args.dilution_rate}")
        print(f"  {'Exp. Time:':<20} {args.exp_time}")
        print("-" * 50)
        print(f"  {'Env ID:':<20} {args.env_id if args.env_id is not None else 'auto'}")
        print(f"  {'GT Range:':<20} {args.gt_range if args.gt_range is not None else 'sin filtro'}")
        print(f"  {'Delay Timing:':<20} {args.delay_timing if args.delay_timing is not None else 'sin filtro'}")
        print("=" * 50)

    run_chemostat(args)

    if args.verbose:
        print(f"[INFO] Guardado en {args.output}")


if __name__ == "__main__":
    args = parse_args()
    main(args)